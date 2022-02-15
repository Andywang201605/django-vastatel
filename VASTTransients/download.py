# ztwang201605@gmail.com

from astropy.table import Table, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.utils.data import clear_download_cache
from astropy.stats import SigmaClip

from astroquery.skyview import SkyView
from astroquery.cadc import Cadc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

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
    hdulists = cadc.get_images(coord, radius, collection='VLASS')
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
    
### Function for wise color-color plot
def plot_wise_cc(position, radius=5):
    """
    Plot a WISE color-color diagram with source at position location overlaid. (same as that in download module)

    Params:
    ----------
    position: SkyCoord, list or tuple
        position of interest
    radius: float or int, 5 by default
        crossmatch radius with WISE catalogue

    Returns:
        fig, ax
    """
    if isinstance(position,tuple) or isinstance(position,list):
        position = SkyCoord(*position, unit=u.deg)

    v = Vizier(columns=['*', '+_r'])
    tablelist = v.query_region(position, radius=radius * u.arcsec, catalog='II/328/allwise')
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    
    ccplot_path = pkg_resources.resource_filename(
        __name__, "./setup/wise_cc.png"
    )

    
    if len(tablelist) > 0:
        result = tablelist[0][0]
        resultcoord = SkyCoord(ra=result['RAJ2000'], dec=result['DEJ2000'], unit=u.deg)
        sep = resultcoord.separation(position).arcsec
        source = pd.Series({'3.4 micron': result['W1mag'],
                            '4.6 micron': result['W2mag'],
                            '12 micron': result['W3mag']})

        w4_12 = source['4.6 micron'] - source['12 micron']
        w3_4 = source['3.4 micron'] - source['4.6 micron']

        ax.scatter(w4_12, w3_4, marker='*', color='magenta', s=100, label=result['AllWISE'])
        im = plt.imread(ccplot_path)
        ax.imshow(im, extent=[-1, 7, -0.5, 4], aspect=2)

        ax.set_xticks([0, 2, 4, 6])
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_xlabel(r'[$4.6\mu m] - [12\mu m$] mag')
        ax.set_ylabel(r'[$3.4\mu m] - [4.6\mu m$] mag')
        ax.set_xlim(-1, 7)
        ax.set_ylim(-0.5, 4)
        ax.set_title('{} - {:.1f} arcsec separation'.format(result["AllWISE"], sep))
        ax.legend()
    else:
        im = plt.imread(ccplot_path)
        ax.imshow(im, extent=[-1, 7, -0.5, 4], aspect=2)

        ax.set_xticks([0, 2, 4, 6])
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_xlabel(r'[$4.6\mu m] - [12\mu m$] mag')
        ax.set_ylabel(r'[$3.4\mu m] - [4.6\mu m$] mag')
        ax.set_title(f"No WISE crossmatch")
        ax.legend()

    return fig, ax

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



