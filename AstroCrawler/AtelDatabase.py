# update django database for atel

######################################################################
##### set up django configuration #####
from django.conf import settings
import django

from webinterface.settings import DATABASES, INSTALLED_APPS
try:
    settings.configure(DATABASES=DATABASES, INSTALLED_APPS=INSTALLED_APPS)
    django.setup()
except RuntimeError:
    print('Settings already configured... continue...')
# import atel models #
from atel.models import AtelMain, AtelSource, AtelSub
from multilambda.models import SourceWeb
######################################################################

from datetime import datetime
import logging
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u

from django.utils import timezone

######## update database for a given atelpage ########
def updataDatabase(atelpage):
    '''
    Update django database for a given atelpage

    Params:
    ----------
    atelpage: AstroCrawler.AtelParser.AtelPage
    '''
    ### check if coords have been updated ###
    logger = logging.getLogger('atelparser.updatedb')
    try: atelpage.coords
    except AttributeError: 
        logger.critical('no coordinate information found, please run `findCoord`')
        raise AttributeError(
        '''Have you ran `findCoord` function already?'''
        )
    ### update atelmain
    mainModel = _findMain_(int(atelpage.atelid), clear=True)

    ### update subjects
    for subject in atelpage.subjects:
        subModel = _findSub_(subject)
        mainModel.atelSub.add(subModel)

    ### update source
    for srccoord in atelpage.coords:
        srcModel = _findSrc_(srccoord)
        mainModel.atelSrc.add(srcModel)

def _findMain_(newID, clear=True):
    '''find if `newID` existed in the dataframe'''
    try: mainModel = AtelMain.objects.get(pk=int(newID))
    except: 
        mainModel = AtelMain(atelID=int(newID)); mainModel.save()
        return mainModel
    
    ### clear subModel and srcModel
    if not clear: return mainModel
    mainModel.atelSub.clear()
    mainModel.atelSrc.clear()
    return mainModel

def _findSub_(newSub):
    '''find a `newSub` in AtelSub database'''
    subList = list(AtelSub.objects.values_list())
    for pk, sub in subList:
        if sub.lower() == newSub.lower():
            return AtelSub.objects.get(pk=int(pk))
    subModel = AtelSub(subject=newSub); subModel.save()
    return subModel

def _findSrc_(coord):
    srccoord = SkyCoord(*coord, unit=u.degree)
    dbsrc = np.array(AtelSource.objects.values_list('id', 'srcra', 'srcdec'))
    if len(dbsrc) > 0:
        dbcoord = SkyCoord(ra = dbsrc[:,1], dec = dbsrc[:, 2], unit = u.degree)
        seps = srccoord.separation(dbcoord).value
        minseq = seps.min(); minloc = seps.argmin()
        if minseq < 5./3600:
            pkid = int(dbsrc[minloc,0])
            srcModel = AtelSource.objects.get(pk=pkid)
            return srcModel
    ### either nothing in db or minseq is larger than 5 arcsec
    taskname = _findWebSrc_(coord)
    srcModel = AtelSource(
        srcra = coord[0],
        srcdec = coord[1],
        taskname = taskname,
    )
    srcModel.save()
    return srcModel

def _findWebSrc_(coord):
    '''check if the source is in SourceWeb database'''
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
    ### either nothing in db or minseq is larger than 5 arcsec
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

## tools for creating parameters ##
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

#######################################################

######## get information from database ########
class AtelDataBase:
    '''
    object for querying information from django atel database
    '''
    def __init__(
        self,
        atelid,
        ):
        self.atelid = int(atelid)
        self.logger = logging.getLogger('atelparser.ateldb')
        self._getMain_()

    def _getMain_(self):
        '''
        fetch main model from database

        Error:
        -----------
        raise ValueError if `atelid` is not in the database
        '''
        try: mainModel = AtelMain.objects.get(pk=self.atelid); self.mainModel = mainModel
        except: 
            self.logger.critical(f'{self.atelid} does not exist in the database...')
            raise ValueError(f'{self.atelid} does not exist in the database...')

    def _getSub_(self):
        '''get subject from database'''
        self.subject = []
        for subModel in self.mainModel.atelSub.all():
            self.subject.append(subModel.subject)
    
    def _getSrc_(self):
        '''get source and taskname...'''
        self.srctask = []
        for srcModel in self.mainModel.atelSrc.all():
            coord = (srcModel.srcra, srcModel.srcdec)
            srcname = _sourcename_(coord)
            self.srctask.append((srcname, srcModel.taskname))

