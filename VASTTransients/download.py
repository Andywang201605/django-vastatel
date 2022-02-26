# ztwang201605@gmail.com

from astropy.table import Table, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.utils.data import clear_download_cache
from astropy.stats import SigmaClip

from astroquery.ipac.ned import Ned
from astroquery.skyview import SkyView
from astroquery.cadc import Cadc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.image_cutouts.first import First

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pkg_resources
import requests
import threading
import io
import os

### Functions for downloading various multiwavelengths data
def _save_hdulist(hdulist, fname, overwrite=True):
    hdulist.writeto(fname, overwrite=overwrite)

### PanSTARRS download handling ###
def getimages_PanSTARRS(ra,dec,size=240,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

def geturl_PanSTARRS(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages_PanSTARRS(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url

### Skymapper download handling ###
def geturl_skymapper(ra, dec, radius):
        """Fetch cutout data via Skymapper API."""

        radius_degree = radius / 3600
        linka = 'http://api.skymapper.nci.org.au/aus/siap/dr2/'
        linkb = f'query?POS={ra:.5f},{dec:.5f}&SIZE={radius_degree:.3f}&BAND=all&RESPONSEFORMAT=CSV'
        linkc = '&VERB=3&INTERSECT=covers'
        sm_query = linka + linkb + linkc

        link = linka + 'get_image?IMAGE={}&SIZE={}&POS={},{}&FORMAT=fits'

        table = requests.get(sm_query.format(ra, dec, radius))
        df = pd.read_csv(io.StringIO(table.text))
        impos = f'{ra:.2f}, {dec:.2f}'
        assert len(df) > 0, f'No Skymapper image at {impos}'
        assert 'Gateway Time-out' not in df.iloc[0], f'Skymapper Gateway Time-out for image at {impos}'

        df = df[df.band == 'z']
        link = df.iloc[0].get_image
        return link

### DECam LS download handling
def geturl_decam(ra, dec, radius, band='g'):
    """Fetch cutout data via DECam LS API.
    credit: Joshua Pritchard

    Params:
    ----------
    radius: float - cutout radius in arcsec
    """

    # Ensure requested image size is less than DECam maximum of 512x512
    size = int(radius / 0.262)
    if size > 512:
       size = 512
       radius = size * 0.262 / 3600

    link = f"http://legacysurvey.org/viewer/fits-cutout?ra={ra}&dec={dec}"
    link += f"&size={size}&layer=dr8&pixscale=0.262&bands={band}"

    return link

### VLASS downloading - not just get link!
def getimage_vlass(ra, dec, radius=60.):
    coord = SkyCoord(ra, dec, unit=u.degree)
    radius = radius * u.arcsec

    cadc = Cadc()
    try: hdulists = cadc.get_images(coord, radius, collection='VLASS')
    except: return -1 # likely source not in the VLASS footprint
    return hdulists

### get vlass epoch
def _getVLASS_epoch(hdulist):
    header = hdulist[0].header
    return '{}.{}'.format(header['FILNAM01'], header['FILNAM02'])

    
def download_archival(ra,dec,radius,survey,savedir, cache=True):
    '''
    Function for downloading archival fits data

    Params:
    ----------
    ra, dec: float
        coordinate for the position of interest
    radius: float
        radius of the fits file in arcsec
    survey: str
        name of survey, values accepted are those for SkyView and `Skymapper`, `PanSTARRS`
    
    '''
    _survey = survey.replace(' ', '_')
    fits_fname = '{}_{}.fits'.format(_survey, int(radius))
    fitspath = os.path.join(savedir, fits_fname)

    if os.path.exists(fitspath):
        # print(fitspath, 'already exists!')
        return 0
    
    if survey == 'VLASS':
        # this one will be different as there will be different versions of VLASS
        hdulists = getimage_vlass(ra, dec, radius)
        if isinstance(hdulists, int): return -1
        for hdulist in hdulists:
            _survey = _getVLASS_epoch(hdulist)
            fits_fname = '{}_{}.fits'.format(_survey, int(radius))
            fitspath = os.path.join(savedir, fits_fname)
            _save_hdulist(hdulist, fitspath)
        return 0

    if survey == 'PanSTARRS':
        urls = geturl_PanSTARRS(ra, dec, size=radius*4,filters="g",format='fits')
        if len(urls) == 0: return -1
        hdulist = fits.open(urls[0]); hdulists = [hdulist]

    elif survey == 'SkyMapper':
        try: url = geturl_skymapper(ra, dec, radius)
        except: return -1
        hdulist = fits.open(url); hdulists = [hdulist]

    elif survey == 'DECam':
        try:
            url = geturl_decam(ra, dec, radius)
            hdulist = fits.open(url)
            hdulists = [hdulist]
        except: return -1
    elif survey == 'FIRST':
        try:
            coord = SkyCoord(ra, dec, unit=u.degree)
            radius = radius * u.arcsec
            hdulists = First.get_images(coord, image_size=2*radius)
        except: return -1
    else:
        try:
            hdulists = SkyView.get_images(position=f'{ra} {dec}',survey=[survey], radius=radius*u.arcsec,cache=cache)
        except:
            return -1
    ### format survey name
    
    _save_hdulist(hdulists[0], fitspath)
    return 0

# these functions for multiple downloading
# def _parse_downloadlist(survey_radius):
#     '''
#     Parsing survey_radius dictionary into a list for downloading

#     Params:
#     ----------
#     survey_radius: dict - keys: str, values: list
#         a dictionary contains survey name and cutout radius (raidus can be list)
    
#     Returns:
#     ----------
#     download_list: list
#         a list containing downloading parameters (element: [survey, radius])
#     '''

#     download_list = []
#     for survey in survey_radius:
#         for radius in survey_radius[survey]:
#             download_list.append([survey, radius])
#     return download_list

def download_archival_multithreading(ra, dec, download_list, savedir, maxthreads=8, cache=False):
    '''
    Download fits image with multiple threads

    Params:
    ----------
    ra, dec: float
        coordinate for the source of interest
    survey_radius: dict - keys: str, values: list
        a dictionary contains survey name and cutout radius 
    savedir: str
        directory for saving fits files
    maxthreads: int, 8 by default
        maximum threads used for downloading simutaneously
    '''

    process_args = []
    thread_args = []
    for download_pairs in download_list:
        ### create a list containing all args for threading
        thread_args.append([ra, dec, download_pairs[1], download_pairs[0], savedir, cache])
        if len(thread_args) >= maxthreads:
            # print(thread_args)
            # print('***********')
            process_args.append(thread_args)
            thread_args = []
    if len(thread_args) != 0: process_args.append(thread_args)

    # print(process_args)

    ### Start threadings
    for thread_args in process_args:
        download_threads = []
        for thread_arg in thread_args:
            th = threading.Thread(target=download_archival, args=thread_arg)
            th.start()
            download_threads.append(th)

        for r in download_threads:
            r.join()

    ### clear cache
    clear_download_cache()
    
#%% 
### Functions for wise color-color plot based on Ned Query
def _getNedTable_(coord, radius = 5.0*u.arcsec):
    if isinstance(coord, list) or isinstance(coord, tuple):
        coord = SkyCoord(*coord, unit=u.degree)
    nedTable = Ned.query_region(coord, radius=radius).to_pandas()
    
    if len(nedTable) == 0: return None
    return nedTable.sort_values('Separation')

def _getNearWISE_(nedTable):
    '''
    get nearest wise source given a ned query result (already ordered by separation)
    '''
    for i, row in nedTable.iterrows():
        rowName = row['Object Name']
        if len(rowName) >= 4:
            if rowName[:4] == 'WISE': 
                separation = row['Separation']*60
                return rowName, separation
    return None, None

def _getWISEVizier_(coord, radius = 5.0*u.arcsec):
    '''Get WISEA from Vizier (there might be some name convention)'''
    if isinstance(coord, list) or isinstance(coord, tuple):
        coord = SkyCoord(*coord, unit=u.degree)
    v = Vizier(columns=['*', '+_r'])
    tablelist = v.query_region(coord, radius=radius, catalog='II/328/allwise')
    if len(tablelist) == 0: return None, None # for nothing
    wisename = 'WISEA {}'.format(tablelist[0][0]['AllWISE'])
    separation = tablelist[0][0]['_r']
    return wisename, separation

def _addNedTag_(nedPhotDf):
    '''
    add tags for a ned WISE photometry dataframe based on `Spatial Mode` and `Qualifiers` columns
    
    Params:
    ----------
    nedPhotDf: pandas.core.frame.DataFrame
        photometry dataframe from ned for WISE (can be the whole dataframe)
    '''
    tagdict = {
        'From fitting to map, Profile-fit': 1,
        'From fitting to map, Profile-fit;variable': 1,
        'Flux in fixed aperture, r=8.25" COG-corrected': 2,
        'Flux in fixed aperture, r=8.25" COG-corrected;variable': 2,
        'Flux in fixed aperture, r=22.0" aperture': 3,
        'Flux in fixed aperture, r=22.0" aperture;variable': 3,
    }
    
    tagLst = []
    for i, row in nedPhotDf.iterrows():
        spatialMode = row['Spatial Mode']; qualifier = row['Qualifiers']
        checkmsg = f'{spatialMode}, {qualifier}'
        tag = tagdict.get(checkmsg)
        if tag is None: tag = 0
        tagLst.append(tag)
    nedPhotDf['PhotoTag'] = tagLst
    return nedPhotDf.copy(deep=True)
        
def _getWISE3_(nedPhotW3):
    '''
    Find out photometry data for WISE3
    Iterate through all available W3 photometry data, order for checking
    - From fitting to map, profile-fit
    - Flux in fixed aperture, r=8.25"
    - Flux in fixed aperture, r=22.0"
    If all available photometry data are upper limits
    Choose one from fitting to map as the one upper limits
    
    Params:
    ----------
    nedPhotW3: pandas.core.frame.DataFrame
        photometry dataframe from ned for WISE3 only
        
    Returns:
    ----------
    detection: Boolean/None
        whether this is a detection or not, NoneType for there is no value
    photoMeasure: float
        photometry measurement for WISE3
    photoMeasureErr: float
        photometry measurement error for WISE3
    photoType: integer
        type of photometry (see `_addNedTag_` function), -1 for no data
    '''
    ### check if phototag column exists, sort the value then
    if 'PhotoTag' not in nedPhotW3: nedPhotW3 = _addNedTag_(nedPhotW3)
    nedPhotW3 = nedPhotW3.sort_values('PhotoTag')
    
    for i, row in nedPhotW3.iterrows():
        photoMeasure = row['Photometry Measurement']
        photoMeasureErr = row['Uncertainty']
        photoType = row['PhotoTag']
        if photoType == 0: continue # continue if photoType is 0 (i.e., not 1,2,3)
        if not np.isnan(photoMeasure): # check if it is detection
            photoMeasureErr = float(photoMeasureErr[3:])
            return True, photoMeasure, photoMeasureErr, photoType
    # if code runs till here -  that means there is no any detection
    # find the one with tag1 as the upper limit. if not, return nothing
    subrow = nedPhotW3[nedPhotW3['PhotoTag'] == 1]
    if len(subrow) > 0: # there is a row
        row = subrow.iloc[0]
        photoMeasureErr = row['Uncertainty'] # in the format of >xxx.xx
        return False, float(photoMeasureErr[1:]), 0.5, 1
    return None, 0, 0, -1 # last, there is no data

def _getWISE12_(nedPhotSub, photoTag=1):
    '''
    Find out photometry data for WISE1 or WISE2 for a given `photoTag`
    
    '''
    ### check if phototag column exists, sort the value then
    if photoTag == -1: return None, 0, 0, -1
    if 'PhotoTag' not in nedPhotSub: nedPhotWSub = _addNedTag_(nedPhotSub)
    subrow = nedPhotSub[nedPhotSub['PhotoTag'] == photoTag]
    if len(subrow) > 0:
        row = subrow.iloc[0]
        photoMeasure = row['Photometry Measurement']
        photoMeasureErr = row['Uncertainty']
        photoType = row['PhotoTag']
        ### check if it is a detection
        if not np.isnan(photoMeasure):
            photoMeasureErr = float(photoMeasureErr[3:])
            return True, photoMeasure, photoMeasureErr, photoType
        else: # if it is not a detection
            photoMeasureErr = row['Uncertainty'] # in the format of >xxx.xx
            return False, float(photoMeasureErr[1:]), 0.5, photoType
    ### if there is no corresponding value
    return None, 0, 0, -1

def _ParseWISEPhot_(nedPhot):
    '''
    Parse a ned WISE photometry table for a given source
    
    Params:
    ----------
    nedPhot: pandas.core.frame.DataFrame
        photometry dataframe from ned
        
    Returns:
    ----------
    '''
    nedPhot = _addNedTag_(nedPhot)
    nedPhotW3 = nedPhot[nedPhot['Observed Passband'] == 'W3 (WISE)']
    nedPhotW2 = nedPhot[nedPhot['Observed Passband'] == 'W2 (WISE)']
    nedPhotW1 = nedPhot[nedPhot['Observed Passband'] == 'W1 (WISE)']
    
    ### parse WISE3 first
    wise3tuple = _getWISE3_(nedPhotW3)
    ### parse WISE2 and WISE1
    wise2tuple = _getWISE12_(nedPhotW2, photoTag=wise3tuple[-1])
    wise1tuple = _getWISE12_(nedPhotW1, photoTag=wise3tuple[-1])
    
    return wise1tuple, wise2tuple, wise3tuple
    
def getWISEMeasure(coord, radius = 5.0*u.arcsec):
    nedTable = _getNedTable_(coord, radius=radius)
    if nedTable is None: return None, None
    wisename, separation = _getNearWISE_(nedTable)
    if wisename is None: 
        wisename, separation = _getWISEVizier_(coord, radius)
        if wisename is None: return None, None
    photoTab = Ned.get_table(wisename, table='photometry').to_pandas()
    return _ParseWISEPhot_(photoTab), separation

def _ColorCalculate_(band1, band2):
    '''
    calculate `colorBand1 - colorBand2`
    
    Params:
    ----------
    band1, band2: tuple from `getWISEMeasure`
        The content within each tuple is:
        - whether it is a detection or not (None/True/False)
        - flux measurement
        - flux measurement error
        - photometry type
        
    Returns:
    ----------
    color, colorErr, colorlimit
    Note: colorlimit take value from [-1, 0, 1, None]
        -1 for lower limit, 0 for an errorbar, 1 for upper limit, None for no data
    '''
    if band1[0] is None or band2[0] is None:
        return None, None, None
    if band1[0] == False and band2[0] is False: # both bands are lower limit
        return None, None, None
    ### at least one band is a detection
    color = band1[1] - band2[1]
    if band1[0] == False:
        colorLimit = -1
        colorErr = 0.5
        return color, colorErr, colorLimit
    if band2[0] == False:
        colorLimit = 1
        colorErr = 0.5
        return color, colorErr, colorLimit
    colorLimit = 0
    colorErr = np.sqrt(band1[2]**2 + band2[2]**2)
    return color, colorErr, colorLimit  

def _getWISEbandcomment_(bandparse):
    '''get comment for a single band'''
    if bandparse[0] is None: return 'nan'
    if bandparse[0] == True: return '{:.3f}'.format(bandparse[1])
    if bandparse[0] == False: return '>{:.3f}'.format(bandparse[1])

def _getWISEcomment_(photoparse):
    
    tagdict = {
        -1: 'no WISE data within 5.0 arcsec',
        1: 'From fitting to map, Profile-fit',
        2: 'Flux in fixed aperture, r=8.25" COG-corrected',
        3: 'Flux in fixed aperture, r=22.0" aperture',
    }
    
    commentmsg = '{}\n'.format(tagdict[photoparse[0][-1]])
    for i in range(3):
        commentmsg += 'W{}:{} '.format(i+1, _getWISEbandcomment_(photoparse[i]))
    return commentmsg
    
    

def CalcWISEColor(photoparse):
    '''
    calculate W1-W2 and W2-W3 for color-color plot
    
    Returns:
    ----------
    
    '''
    
    if photoparse is None: return None, None, 'no WISE data within 5.0 arcsec'
    
    ### form comments on WISE detections
    commentmsg = _getWISEcomment_(photoparse)
    
    ### calculate colors
    Xcolor, XcolorErr, XcolorLim = _ColorCalculate_(photoparse[1], photoparse[2]) # W2-W3
    Ycolor, YcolorErr, YcolorLim = _ColorCalculate_(photoparse[0], photoparse[1]) # W2-W3
    
    return (
        (Xcolor, XcolorErr, XcolorLim),
        (Ycolor, YcolorErr, YcolorLim),
        commentmsg,
    )

def _addWISEcolor_(ax, photoparse, color='darkred'):
    colorX, colorY, comment = CalcWISEColor(photoparse)
    
    ax.set_title(comment)
    
    spanStyle = {'color':color, 'alpha':0.1,}
    markStyle = {'marker':'*', 'markersize':15, 'color': color}
    if colorX is None or colorY is None: return ax
    if colorX[-1] is None and colorY[-1] is None:
        return ax # do nothing
    if colorX[-1] is None: # colorX is none
        if colorY[-1] == -1:
            ax.axhspan(ymin=colorY[0], ymax=5, **spanStyle)
            ax.errorbar(
                x = 3, y = colorY[0], yerr=colorY[1], xerr=6, lolims=True, 
                **markStyle,
            )
            return ax
        if colorY[-1] == 0:
            ax.axhspan(ymin=colorY[0]-colorY[1], ymax=colorY[0]+colorY[1], **spanStyle)
            ax.errorbar(
                x = 3, y = colorY[0], yerr = colorY[1], xerr=6, 
                **markStyle,
            )
            return ax
        if colorY[-1] == 1:
            ax.axhspan(ymin=-1, ymax=colorY[0], **spanStyle)
            ax.errorbar(
                x = 3, y = colorY[0], yerr=colorY[1], xerr=6, uplims=True, 
                **markStyle,
            )
            return ax
        return ax
    if colorY[-1] is None: # colorY is none
        if colorX[-1] == -1:
            ax.axvspan(xmin=colorX[0], xmax=8, **spanStyle)
            ax.errorbar(
                x = colorX[0], y = 1.75, xerr=colorX[1], yerr=6, xlolims=True, 
                **markStyle,
            )
            return ax
        if colorX[-1] == 0:
            ax.axvspan(xmin=colorX[0]-colorX[1], xmax=colorX[0]+colorX[1], **spanStyle)
            ax.errorbar(
                x = colorX[0], y = 1.75, xerr = colorX[1], yerr=6, 
                **markStyle,
            )
            return ax
        if colorX[-1] == 1:
            ax.axvspan(xmin=-1, xmax=colorX[0], **spanStyle)
            ax.errorbar(
                x = colorX[0], y = 1.75, xerr=colorX[1], yerr=6, xuplims=True, 
                **markStyle,
            )
            return ax
        return ax
    ### for both colors are not none
    uplim, lolim, xuplim, xlolim = False, False, False, False
    if colorX[-1] == -1: xlolim = True
    if colorX[-1] == 1: xuplim = True
    if colorY[-1] == -1: lolim = True
    if colorY[-1] == 1: uplim = True
    
    ax.errorbar(
        x = colorX[0], y = colorY[0], xerr=colorX[1], yerr=colorY[1],
        lolims=lolim, uplims=uplim, xlolims=xlolim, xuplims=xuplim,
        **markStyle,
    )
    return ax

def plot_wise_cc(coord, radius=5.0*u.arcsec):
    ### query NED ###
    photoparse, separation = getWISEMeasure(coord, radius = radius)
    labelmsg = f'no data' if separation is None else f'{separation:.2f} arcsec away'

    ccplot_path = pkg_resources.resource_filename(
        __name__, "./setup/wise_cc.png"
    )
    
    im = plt.imread(ccplot_path)
    fig = plt.figure(figsize=(6, 6), facecolor='white')
    ax = fig.add_subplot(111)

    ax.imshow(im, extent=[-1, 7, -0.5, 4], aspect=2)
    ax = _addWISEcolor_(ax, photoparse, color='darkred')
    
    ax.set_xticks([0, 2, 4, 6])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlabel(r'[$4.6\mu m] - [12\mu m$] mag')
    ax.set_ylabel(r'[$3.4\mu m] - [4.6\mu m$] mag')
    ax.set_xlim(-1, 7)
    ax.set_ylim(-0.5, 4)
    
    ax.scatter([], [], marker='*', color='darkred', label=labelmsg)
    ax.legend(fontsize=12)
    
    return fig, ax
#%%

### Functions for getting archival crossmatch
def get_archival_crossmatch(coord, catalog, radius):
    '''
    Get matches from archival catalogue

    Params:
    ----------
    coord: SkyCoord, tuple or list:
        position of interest
    catalog: str
        catalog reference code from Vizier
    radius: float
        crossmatch radius in arcsec
    
    Returns:
    ----------
    result: int or astropy.Table
    '''
    if isinstance(coord,tuple) or isinstance(coord,list):
        coord = SkyCoord(*coord, unit=u.deg)

    v = Vizier(columns=['*', '+_r'])
    tablelist = v.query_region(coord, radius=radius*u.arcsec, catalog=catalog)
    if len(tablelist) == 0: return -1
    return tablelist[0][0]

def _parse_Vizier_result(r, survey):
    '''
    Parse result from Vizier and output a list contains frequency and its coresponding flux

    Parmas:
    ----------
    r: astropy.Table
        Vizier result from `get_archival_crossmatch` function
    survey: str
        name of the survey, values accepted: `AT20G`, `GLEAM`, `SUMSS`, `NVSS`, `TGSS`

    Returns:
    ----------
    values : list (frequency and flux pairs)
    '''
    if survey == 'AT20G':
        return [[5,r['S5']],[8,r['S8']],[20,r['S20']]]
    if survey == 'GLEAM':
        return [[0.076,r['Fp076']],[0.084,r['Fp084']],[0.092,r['Fp092']],[0.099,r['Fp099']],[0.107,r['Fp107']],
                [0.115,r['Fp115']],[0.122,r['Fp122']],[0.13,r['Fp130']],[0.143,r['Fp143']],[0.151,r['Fp151']],
                [0.158,r['Fp158']],[0.166,r['Fp166']],[0.174,r['Fp174']],[0.181,r['Fp181']],[0.189,r['Fp189']],
                [0.197,r['Fp197']],[0.204,r['Fp204']],[0.212,r['Fp212']],[0.22,r['Fp220']],[0.227,r['Fp227']]]
    if survey == 'SUMSS':
        return [[0.843,r['Sp']]]
    if survey == 'NVSS':
        return [[1.4,r['S1.4']]]
    if survey == 'TGSS':
        return [[0.15,r['Speak']]]

def _get_noise_image(fitspath, index=0):
    '''
    Get noise of one fits image (as a whole)

    Params:
    ----------
    fitspath: str
        path for the fits file
    index: int, 0 by default
        index of the data

    Returns:
    ----------
    noise: float
    '''
    with fits.open(fitspath) as hdulist:
        data = hdulist[index].data

    sigmaclip = SigmaClip()
    noise = sigmaclip(data.ravel()).std()
    return noise

def _get_nondetection_refimage(fitspath, survey):
    '''
    get reference image for a non-detection in archival

    Params:
    ----------
    fitspath: str
        path for downloading all multiwavelength images
    survey: str
        name of the survey, value accepted: `AT20G`, `GLEAM`, `SUMSS`, `NVSS`, `TGSS`

    Returns:
    ----------
    imagepath: str
        full path of a fits file
    '''
    assert survey in ['SUMSS', 'NVSS', 'TGSS', 'GLEAM']

    if survey == 'SUMSS':
        fits_fname = f'{fitspath}/SUMSS_843_MHz_600.fits'
    elif survey == 'NVSS':
        fits_fname = f'{fitspath}/NVSS_600.fits'
    elif survey == 'TGSS':
        fits_fname = f'{fitspath}/TGSS_ADR1_600.fits'
    elif survey == 'GLEAM':
        fits_fname = f'{fitspath}/GLEAM_170-231_MHz_600.fits'

    return fits_fname

def _get_survey_reffreq(survey):
    '''
    Get reference frequency for a survey

    Params:
    ----------
    survey: str
        name of the survey

    Returns:
    ----------
    freq: float
        frequency in GHz
    '''
    freqdict = {'SUMSS': 0.843, 'NVSS':1.4, 'TGSS':0.15, 'GLEAM':0.2}
    return freqdict.get(survey)

def get_archival_data(coord, catalogs, fitspath=None, sigma=5):
    '''
    get all archival data from a list of catalogs you provided

    Params:
    ----------
    coord: tuple, list or SkyCoord
        position of the source
    catalogs: dict, keys accepted (`AT20G`, `GLEAM`, `SUMSS`, `NVSS`, `TGSS`)
        values are a list, which are (in order) `catalog reference code in Vizier`, `year of the observation`, `search radius`
    fitspath: str
        for a non-match, search for any fits file already downloaded and estimate the upperlimit
    sigma: float or int
        upper limit in the final data for a non-detection

    Returns:
    ----------
    archival_data: list of strings
        can be written to a csv file - header: survey, frequency, date, flux, type
    '''
    archival_data = ['survey,freq,date,flux,type\n']
    for survey in catalogs:
        reference_code = catalogs[survey][0]
        yearobs = catalogs[survey][1]
        searchradius = catalogs[survey][2]

        viziertable = get_archival_crossmatch(coord, reference_code, searchradius)

        if isinstance(viziertable, int): # no match found
            if survey == 'AT20G':
                continue
            ### for other survey
            imagepath = _get_nondetection_refimage(fitspath, survey)
            if os.path.exists(imagepath):
                try:
                    noise = _get_noise_image(imagepath)*1e3
                    freq = _get_survey_reffreq(survey)
                    archival_data.append(f'{survey},{freq},{yearobs},{sigma*noise},0\n')
                except:
                    continue
            else:
                continue

        ## for a match
        else:
            vizierlist = _parse_Vizier_result(viziertable, survey)
            for freq, flux in vizierlist:
                archival_data.append(f'{survey},{freq},{yearobs},{flux},1\n')
    return archival_data

def query_simbad(coord, radius=40):
    '''
    Query result from simbad

    Params:
    ----------
    coord: tuple, list or SkyCoord
        coordinate of the source
    radius: int or float
        search radius in arcsecond

    Returns
    ----------
    result: astropy.Table or NoneType
    '''
    if isinstance(coord, tuple) or isinstance(coord, list):
        coord = SkyCoord(*coord, unit=u.deg)
    query = Simbad.query_region(coord,radius*u.arcsec)
    if isinstance(query,Table):
        sep_col = Column(coord.separation(SkyCoord(query['RA'],query['DEC'],unit=(u.hourangle, u.deg))).arcsec,name='separation')
        query.add_column(sep_col,0)
        return query

def get_simbad_url(coord, radius=40):
    '''
    Get url for SIMBAD searching

    Params:
    ----------
    coord: tuple, list or SkyCoord
        coordinate of the source
    radius: int or float
        search radius
    
    Returns:
    ----------
    simbadurl: str
    '''
    if isinstance(coord, tuple) or isinstance(coord, list):
        ra, dec = coord
    else:
        ra = coord.ra.value; dec = coord.dec.value
    return f"http://simbad.u-strasbg.fr/simbad/sim-coo?Coord={ra}d+{dec}d&Radius={radius}&Radius.unit=arcsec&submit=submit+query"



