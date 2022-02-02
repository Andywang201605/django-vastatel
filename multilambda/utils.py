### django relavant
from django.http import request
from django.utils import timezone

from VASTTransients.source import VASTSource
from .models import SourceWeb

### package development
import pkg_resources
from .config import *

### other packages
# normal packages
from datetime import datetime
import numpy as np
import threading
import json
import re
import os
# astropy modules
from astropy.coordinates import SkyCoord
import astropy.units as u

# export PDF
from weasyprint import HTML, CSS

######## tools for parsing input ########
def _parsecoord_(coordinput, galactic=False):
    '''
    function for parsing coordinate input. 
    Possible input:
        1) degree units
        2) sexagesimal units

    Params:
    ----------
    coordinput: str
        coordinate input from users
    galactic: boolean
        choose to use galactic frame if this True, False by default

    Returns:
    ----------
    coordinate: tuple (float, float)
        ra and dec for the given cooridnate
    '''
    frame = 'galactic' if galactic else 'icrs'
    
    coordinput = coordinput.lower().strip() # change to lower case

    ######## parse sexagesimal units ########
    rasex = r'''(?P<raSex>(
        (?P<raHrs>\d+)(:\s*|\s|h\s*)
        (?P<raMin>\d+)(:\s*|\s|m\s*)
        (?P<raSec>\d+\.?\d*)(s?)
    ))'''
    decsex = r'''(?P<decSex>(
        (?P<decDeg>[\+|-]?\d+)(:\s*|\s|d\s*)
        (?P<decMin>\d+)(:\s*|\s|m\s*|'\s*)
        (?P<decSec>\d+\.?\d*)(s?|"\s*)
    ))'''
    restr = '^{},?\s+{}$'.format(rasex, decsex)
    repat = re.compile(restr, re.X | re.S)
    rematch = re.search(repat, coordinput) # None if there is no match
    if rematch is not None:
        if galactic: raise ValueError('Sexagesimal units are not supported for Galactic frame...')
        matchdict = rematch.groupdict()
        coordsex = '{}h{}m{}s {}d{}m{}s'.format(
            matchdict['raHrs'], matchdict['raMin'], matchdict['raSec'],
            matchdict['decDeg'], matchdict['decMin'], matchdict['decSec'],
        )
        coord = SkyCoord(coordsex)
        return coord.ra.value, coord.dec.value

    ######## parse degree units ########
    restr = r'^(?P<rafloat>\d+\.?\d*),?\s+(?P<decfloat>[\+|-]?\d+\.?\d*)$'
    repat = re.compile(restr)
    rematch = re.search(repat, coordinput) # None if there is no match
    if rematch is not None:
        matchdict = rematch.groupdict()
        coord = SkyCoord(
            float(matchdict['rafloat']),
            float(matchdict['decfloat']),
            unit = u.degree,
            frame = frame,
        )
        if galactic: coord = coord.icrs # convert to icrs frame if the original frame is galactic
        return coord.ra.value, coord.dec.value

    raise ValueError('coordinate {} cannot parsed successfully...'.format(coordinput))


######## loading relavant files ########
def _loadjsonsurveys_():
    surveyjsonpath = pkg_resources.resource_filename(
            __name__, "./setups/instrument_survey.json"
        )
    with open(surveyjsonpath) as fp:
        surveyjson = json.load(fp)
    return surveyjson

def _loadjsonradius_():
    radiusjsonpath = pkg_resources.resource_filename(
            __name__, "./setups/instrument_radius.json"
        )
    with open(radiusjsonpath) as fp:
        radiusjson = json.load(fp)
    return radiusjson

def _loadinstruments_():
    surveyjson = _loadjsonsurveys_()
    radiusjson = _loadjsonradius_()

    ### combine these two together ###
    instrument_form = []
    for instrument, surveys in surveyjson.items():
        radius = radiusjson.get(instrument)
        if radius is None: radius = 600
        instrument_form.append((instrument, radius, surveys))
    return instrument_form

def _parserequest_(request):
    '''
    function for parsing post from previous form web

    Params:
    ----------
    request: django HttpRequest object

    Retunrs:
    ----------
    request_info: dict
        a dictionary that contains all information parsed from request.POST
        `coord`: tuple, coordinate of the source in R.A. and Decl.
        `todownload`: list (for each element, (survey[str], radius[float]))
    '''
    ### get coordinate and parse it to ra, dec ###
    galactic = False if request.POST['coord_box'] == '1' else True
    srccoord = _parsecoord_(request.POST['coord'], galactic=galactic)
    ### get survey radius paris ###
    surveyjson = _loadjsonsurveys_()
    download_list = [] # whole list for downloading based on POST
    for band in surveyjson:
        surveys = surveyjson.get(band)
        radius = float(request.POST.get(f'{band}_r'))
        selected = request.POST.getlist(band)
        if selected is None: continue
        ## get information ##
        banddown = [(surveys[int(index) - 1], radius) for index in selected]
        download_list.extend(banddown)
    request_info = {'coord': srccoord, 'todownload': download_list}

    return request_info

######## tools for creating parameters ########
def _sourcename_(coord):
    '''
    get source name based on coordinate - Jxxxx+/-xxxx

    Params:
    ----------
    coord: tuple-like (float, float)

    Returns:
    coordname: str
    '''
    sign = '+' if coord[1] > 0 else '-' # check sign of decl first
    coord = SkyCoord(*coord, unit=u.degree)
    namera = '{:0>2.0f}{:0>2.0f}'.format(coord.ra.hms.h, coord.ra.hms.m)
    namede = '{:0>2.0f}{:0>2.0f}'.format(abs(coord.dec.dms.d), abs(coord.dec.dms.m))
    return f'J{namera}{sign}{namede}'

def _gettaskname_(coord):
    '''
    function for creating task name for a given coordinate. 
    The task name is given by <time_stamp>_<sourcename>

    Params:
    ----------
    coord: tuple-like (float, float)

    Returns:
    ----------
    taskname: str
    '''
    tsnow = datetime.now().timestamp()
    srcname = _sourcename_(coord)
    return '{:.0f}_{}'.format(tsnow, srcname)

# assign task name to a specific source #
def _assigntask_(coord):
    '''
    function for assigning the source to a specific source. 
    Check if there is any available file already there, if not, update the database

    Params:
    ----------
    coord: tuple (float, float)

    Returns:
    ----------
    taskname: str
        task name for this source
    '''
    # assign task if the source is close to a file already existed...
    srccoord = SkyCoord(*coord, unit = u.degree)
    dbsrc = np.array(SourceWeb.objects.values_list('id', 'srcra', 'srcdec'))
    if len(dbsrc) > 0:
        dbcoord = SkyCoord(ra = dbsrc[:,1], dec = dbsrc[:, 2], unit = u.degree)
        seps = srccoord.separation(dbcoord).value
        minseq = seps.min(); minloc = seps.argmin()
        if minseq < 5./3600.:
            pkid = int(dbsrc[minloc,0])
            srcweb = SourceWeb.objects.get(pk = pkid)
            srcweb.taskdate = timezone.now()
            taskname = srcweb.taskname
            srcweb.save()
            return taskname

    ### if a source nearby does not exists... 
    taskname = _gettaskname_(coord)
    srcweb = SourceWeb(
        taskname = taskname,
        taskdate = timezone.now(),
        srcra = coord[0],
        srcdec = coord[1],
        fileExist = True,
    )

    srcweb.save()

    return taskname

def _sourceanalysis_(coord, taskname, download_list):
    vastsource = VASTSource(
        coord,
        f'{LAMBDADIR}/{taskname}/',
        download_list,
    )
    vastsource.sourceAnalysis()
    ## TODO: add export PDF function ## - done!
    html = HTML(f"http://0.0.0.0:8021/multilambda/{taskname}/webpage/")
    html.write_pdf(
        f'{LAMBDADIR}/{taskname}/multiweb_{taskname}.pdf',
        stylesheets=[
            CSS(string='''@page {
                size: A4; 
                margin: 0mm; 
                }''')
        ],
        presentational_hints=True,
    )
    return

def _savejson_(json_data, taskname, jsonname):
    if not os.path.exists(f'{LAMBDADIR}/{taskname}/'):
        os.makedirs(f'{LAMBDADIR}/{taskname}/')
    jsonfile = f'{LAMBDADIR}/{taskname}/{jsonname}.json'
    with open(jsonfile, 'w') as fp:
        json.dump(json_data, fp)

### open a thread ###
def _openthread_(func, args=(), kwargs={}):
    '''
    function for opening a new thread for a given function with specific args and kwargs
    Warning: this function may be changed in the future
    There are possibilities that DoS attack!!!

    Params:
    ----------
    func: function
    args: tuple
        arguments passed to the `func` function
    kwargs: dict
        keyword arguments passed to the `func` function

    Returns:
    ----------
    Nonetype
    '''
    thread = threading.Thread(
        target = func,
        args = args,
        kwargs = kwargs,
        daemon = True,
    )
    thread.start()
    return

### functions for status page ###
def _checkfileexsits_(filepath):
    if os.path.exists(filepath):
        return True
    return False

def _loadrequestjson_(taskname, filenameconvert = False):
    taskdir = f'{LAMBDADIR}/{taskname}/'
    with open(f'{taskdir}/request.info.json') as fp:
        requestjson = json.load(fp)
    ### convert filename if True ###
    if filenameconvert:
        _todownload = [(i.replace(' ', '_'), int(j)) for i, j in requestjson['todownload']]
        requestjson['todownload'] = _todownload
    return requestjson

def _checkdownloadstatus_(taskname):
    taskdir = f'{LAMBDADIR}/{taskname}/'
    status = {}
    ## check if measurements pickle has been created ##
    status['vasttools.query'] = _checkfileexsits_(f'{taskdir}/source_measurement.pickle')
    status['vast.cutouts'] = _checkfileexsits_(f'{taskdir}/img/StokesI_300.jpg')
    ## check downloading status ##
    requestjson = _loadrequestjson_(taskname)
    for survey, radius in requestjson['todownload']:
        survey = survey.replace(' ', '_'); radius = int(radius)
        status[f'{survey}_r{radius}.download'] = _checkfileexsits_(f'{taskdir}/img/{survey}_{radius}.fits')
    return status

### functions for webpage ###
def _getmulticutoutpath_(taskname):
    taskdir = f'{LAMBDADIR}/{taskname}/'
    requestinfo = _loadrequestjson_(taskname, filenameconvert=True)
    overlaypath = []
    for survey, radius in requestinfo['todownload']:
        if _checkfileexsits_(f'{taskdir}/img/{survey}_{radius}.png'):
            overlaypath.append([survey, f'{taskname}/img/{survey}_{radius}.png'])
        else:
            overlaypath.append([survey, "multilambda/noimage.jpg"])
    return overlaypath

def _checkmeasure_(taskname):
    taskdir = f'{LAMBDADIR}/{taskname}/'
    if os.path.exists(f'{taskdir}/source_measurement.pickle'):
        return True
    return False