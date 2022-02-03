### function for source analyzation etc.

from weasyprint import HTML, CSS
import pandas as pd
import numpy as np
import pkg_resources
import json
import os

from .config import *
from . import crawlerUtils

from astropy.coordinates import SkyCoord
import astropy.units as u

from VASTTransients.source import VASTSource

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

def _constructDownlst_():
    surveyjson = _loadjsonsurveys_()
    radiusjson = _loadjsonradius_()

    ### construct download list ###
    downloadLst = []
    for band in radiusjson:
        radius = radiusjson[band]
        banddown = [[survey, radius] for survey in surveyjson.get(band)]
        downloadLst.extend(banddown)
    return downloadLst

def _runVASTTransient_(coord, taskname):
    '''run VASTTransient analysis for a given source'''
    download_list = _constructDownlst_()
    ### save request information for django view usage ###
    request_info = {'coord': coord, 'todownload': download_list}
    workdir = f'{LAMBDADIR}/{taskname}/'; crawlerUtils.makedirs(workdir)
    with open(f'{workdir}/request.info.json', 'w') as fp:
        json.dump(request_info, fp)
    ### analyze source
    vastsource = VASTSource(coord, workdir, download_list)
    vastsource.sourceAnalysis()
    ### search for nearby VAST sources 
    _findVASTrun_(taskname, radius=60.)
    ### export HTMLs
    html = HTML(f"http://0.0.0.0:8021/multilambda/{taskname}/webpage/")
    html.write_pdf(
        f'{LAMBDADIR}/{taskname}/multiweb_{taskname}.pdf',
        stylesheets=[
            CSS(string='''@page {
                size: A4 landscape; 
                margin: 0mm; 
                }''')
        ],
        presentational_hints=True,
    )
    return

##### filter VAST run and find associated sources #####
def _loadrequestjson_(taskname, filenameconvert = False):
    taskdir = f'{LAMBDADIR}/{taskname}/'
    with open(f'{taskdir}/request.info.json') as fp:
        requestjson = json.load(fp)
    ### convert filename if True ###
    if filenameconvert:
        _todownload = [(i.replace(' ', '_'), int(j)) for i, j in requestjson['todownload']]
        requestjson['todownload'] = _todownload
    return requestjson
##-- box search -- ##
# use skycoord will be extremely slow when using astropy.coordinates.SkyCoord.separation
def _calRAangle_(ra1, ra2):
    '''
    calculate the RA angle between two points on the sky with RA = ra1 and ra2 respectively
    '''
    rasep = abs(ra1 - ra2) % 360.
    return 360. - rasep if rasep > 180. else rasep
    
def _boxSearch_(ra, dec, srcdf, radius=60., racol='wavg_ra', deccol='wavg_dec'):
    radius = radius / 3600. # convert to degree
    ### filter dec
    declo = max(-90, dec - radius); dechi = min(90, dec + radius)
    findf = srcdf[(srcdf[deccol] >= declo) & (srcdf[deccol] <= dechi)].copy(deep=True)
    ### filter ra
    rasep = findf[racol].apply(lambda radf:_calRAangle_(ra, radf))
    findf = findf[rasep < min(radius/np.cos(np.deg2rad(dec)), 360)].copy(deep=True)
    return findf
###### --------- ######
def _pipelink_(srcid):
    return f'<a href="https://dev.pipeline.vast-survey.org/sources/{srcid}/">{srcid}</a>'

def _findVASTrun_(taskname, radius=60., overwrite = False):
    if os.path.exists(f'{LAMBDADIR}/{taskname}/VASTrun.source.json') and not overwrite: return
    srcdf = pd.read_pickle(VASTRUNPICKLE)
    ### load coord from taskname
    requestinfo = _loadrequestjson_(taskname, filenameconvert=True)
    coord = requestinfo['coord']
    ### start filtering ###
    findf = _boxSearch_(coord[0], coord[1], srcdf, radius = radius)
    findf['separation'] = SkyCoord(*coord, unit=u.degree).separation(
        SkyCoord(findf['wavg_ra'], findf['wavg_dec'], unit=u.degree)
        ).value*3600. # add separation column, in unit of arcseconds
    findf = findf.sort_values('sep')
    ### if findf is too long - select the first 10 ###
    if len(findf) > 10: findf = findf.iloc[:10]
    findf['pipelink'] = pd.Series(findf.index).apply(_pipelink_).to_numpy()
    ### dump information to json file ###
    findf.to_json(f'{LAMBDADIR}/{taskname}/VASTrun.source.json')

########################
