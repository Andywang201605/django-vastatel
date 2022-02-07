# ztwang2016@gmail.com
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip
from astropy import units as u
from astropy.visualization import ZScaleInterval, ImageNormalize

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astroquery.vizier import Vizier

import json
import pkg_resources

### Handle the fits file, Perform the cutout
class FITSIMAGE:
    '''
    Class for fits files
    '''

    def __init__(self, fitspath, index = 0):
        '''
        Initiate function for FITSIMAGE class

        Params:
        ----------
        fitspath: str
            Path for the fits file
        index: int, 0 by default
            Index of the data to extract from
        '''

        with fits.open(fitspath) as hdu:
            self.data = hdu[index].data
            self.header = hdu[index].header

        self.wcs = WCS(self.header).celestial


    def cutout(self, ra, dec, radius=20.):
        '''
        Make a cutout with a radius of `radius` arcsec centered at coordinate (ra, dec)

        Params:
        ----------
        ra, dec: float
            coordinate for the position of interests
        radius: float, 20.0 by default
            radius of the cutout

        Returns:
        ----------
        cutout: Cutout2D object
        '''
        coord = SkyCoord(ra, dec, unit=u.deg)
        size = (radius*u.deg / 1800., radius*u.deg / 1800.)
        return Cutout2D(np.squeeze(self.data), coord, size, wcs=self.wcs)

### Functions for plotting
def plot_fits(data, ax, norms=None, **kwargs):
    '''
    Plot fits data on a given axis.

    Params:
    ----------
    data: numpy.ndarray
        Data to be plotted
    ax: matplotlib.axis
        Axis to be plotted on
    norms: astropy.visualization.ImageNormalize or Nonetype
    **kwargs:
        arguments passed to plt.imshow function

    Returns:
    ---------
    ax, im
    '''
    ### set default values
    kwargs.setdefault('cmap', 'gray_r')

    ### set Zscale normalization
    if not isinstance(norms,ImageNormalize):
        norms = ImageNormalize(data,interval=ZScaleInterval())

    im = ax.imshow(data,norm=norms,**kwargs)

    return ax, im

def source_crosshair(coord, ax, plots='udrl', sep=3, length=5, **kwargs):
    '''
    Add crosshair in a plot (axis with given wcs in)

    Params:
    ---------
    coord: tuple, list, astropy.coordinate.SkyCoord
        Coordinate to put the crosshair
    ax: matplotlib.axis
        axis with wcs projection
    plots: str, 'urdl' by default
        position for plotting detection crosshair. u - up, d - down, l - left, r - right
    sep: float or int
        seperation between crosshairs
    length: float or int
        length of a single part of crosshair
    **kwargs: argument passed to plt.plot function

    Returns:
    ---------
    ax
    '''
    if isinstance(coord,tuple) or isinstance(coord,list):
        coord = SkyCoord(*coord, unit=u.deg)

    ### transfer skycoord to frame position
    src_x, src_y = ax.wcs.world_to_array_index(coord)

    ### plot crosshairs
    if 'u' in plots:
        ax.plot([src_x,src_x],[src_y+sep,src_y+sep+length],**kwargs)
    if 'd' in plots:
        ax.plot([src_x,src_x],[src_y-sep,src_y-sep-length],**kwargs)
    if 'r' in plots:
        ax.plot([src_x+sep+length,src_x+sep],[src_y,src_y],**kwargs)
    if 'l' in plots:
        ax.plot([src_x-sep-length,src_x-sep],[src_y,src_y],**kwargs)
    return ax

def plot_contour(data, wcs, ax, numlevels=5, **kwargs):
    '''
    Plot contour on an axis with wcs projection

    Params:
    ----------
    data: numpy.ndarray
        data for contour
    wcs: astropy.wcs.WCS
        wcs for contour
    ax: matplotlib.axis
        axis with a wcs projection
    numlevels: int
        number of levels for the contour (if levels are presented in **kwargs, ignore this param)
    **kwargs: arguments passed to plt.contour function

    Returns:
    ----------
    ax, ct
    '''
    ### make levels
    if not kwargs.get('levels'):
        sigmaclip = SigmaClip() # set a SigmaClip instance
        filtered_data = sigmaclip(data)

        ### calculate some metrics
        maximum_value = data.max()
        std_value = filtered_data.std()
        maximum_snr = maximum_value / std_value

        ### get levels
        snr_seperation = (maximum_snr - 3) / numlevels
        levels = (np.arange(numlevels) + 3)*snr_seperation*std_value

        kwargs['levels'] = levels

    ct = ax.contour(data, transform=ax.get_transform(wcs), **kwargs)

    return ax, ct

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
        ax.set_xlim(-1, 7)
        ax.set_ylim(-0.5, 4)
        ax.set_xlabel(r'[$4.6\mu m] - [12\mu m$] mag')
        ax.set_ylabel(r'[$3.4\mu m] - [4.6\mu m$] mag')
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

### Functions for plotting multi-epochs images
def _get_figure_layout(numaxes, ncols, colwidth=5, rowwidth=5):
    '''
    Work out the figure size for multi-axes plot

    Params:
    ----------
    numaxes: int
        number of axes are there in the figure
    ncols: int
        number of columns in the figure
    colwidth, rowwidth: float or int
        width of column and row for a single axis
    
    Returns:
    ----------
    subplot_layout, figsize: tuple
    '''
    nrows = int((numaxes - 1) // ncols) + 1
    subplot_layout = (nrows, ncols)
    figsize = (ncols*colwidth, nrows*rowwidth)
    return subplot_layout, figsize

def plot_multiepoch_cutout(ra, dec, measurements, images, radius=300, imagepath_col='path', imageid_col='image_id', time_col='time', samescale=True):
    '''
    Plot multi-epoch cutout in one image 
    in order to plot stokesV plot, you can add a row in images dataframe and pass the corresponding stokesV path in imagepath_col

    Params:
    ----------
    ra, dec: float
        position of the source
    measurements: pandas.DataFrame
        dataframe contains all measurements, image_id, time of observation etc.
    images: pandas.DataFrame or NoneType #todo
        dataframe contains all images information
    radius: int or float, 300 by default
        radius of the cutout, in arcsec
    imagepath_col: str
        column that saves image path in images dataframe
    imageid_col: str
        column that saves imageid in measurements dataframe
    time_col: str
        column that saves time of the observation in measurements daraframe
    samescale: bool, True by default
        If the image use the same scale and limit for all subplots or not

    Returns:
    ----------
    fig
    '''
    ### sort the measurements based on time
    sorted_measurements = measurements.sort_values('time')
    nepochs = len(sorted_measurements)
    subplot_layout, figsize = _get_figure_layout(nepochs, 4)

    ### start plotting
    fig = plt.figure(figsize=figsize, facecolor='w')

    norms = None; plotcount = 1
    for i, row in sorted_measurements.iterrows():
        imageid = row[imageid_col]
        imagepath = images.loc[imageid][imagepath_col]

        try: # handle fits file doesnot exist
            fitsimage = FITSIMAGE(imagepath)
            cutout = fitsimage.cutout(ra, dec, radius)
        except:
            plotcount += 1
            continue

        if (norms == None) or (samescale == False):
            norms = ImageNormalize(cutout.data,interval=ZScaleInterval())

        ### add subplot
        ax = fig.add_subplot(*subplot_layout, plotcount, projection=cutout.wcs)
        ax, im = plot_fits(cutout.data, ax, norms)
        ax = source_crosshair((ra, dec), ax, color='red', sep=radius/12., length=radius/12.)

        ax.set_title(row['time'])

        # set axislable invisible
        ax.coords[0].set_axislabel(' ')
        ax.coords[1].set_axislabel(' ')

        plotcount += 1


    return fig

### Functions for lightcurves
def plot_VAST_lightcurve(measurements):
    '''
    Plot lightcurve from VAST measurements

    Params:
    ----------
    measurements: pandas.DataFrame
        dataframe contains all measurements

    Returns:
    ----------
    fig, ax
    '''
    ### read lightcurve setup from jsonfile
    lightcurve_setup_path = pkg_resources.resource_filename(
        __name__, './setup/lightcurve_setup.json'
    )

    with open(lightcurve_setup_path) as fp:
        lightcurve_setting = json.load(fp)

    uplim = lightcurve_setting['uplim'] # upper limit threshold to plot
    lightcurve_setting['uplim_style'].setdefault('label', f'{uplim} sigma')

    ### plot lightcurve
    fig = plt.figure(figsize=(12, 6), facecolor='w')
    ax = fig.add_subplot(111)

    ### forced measurements
    forced_measure = measurements[measurements['forced']]
    ax.errorbar(
        forced_measure['time'],
        uplim*forced_measure['local_rms'],
        yerr = forced_measure['local_rms'],
        **lightcurve_setting['uplim_style'],
    )
    ax.errorbar(
        forced_measure['time'],
        forced_measure['flux_peak'],
        yerr = forced_measure['local_rms'],
        **lightcurve_setting['forced_style']
    )

    ### selavy measurements
    selavy_measure = measurements[~measurements['forced']]
    ax.errorbar(
        selavy_measure['time'],
        selavy_measure['flux_peak'],
        yerr = selavy_measure['local_rms'],
        **lightcurve_setting['selavy_style']
    )

    ax.legend()

    return fig, ax

def plot_archival_lightcurve(archival_measures, measurements=None):
    '''
    Plot archival lightcurve for the source

    Params:
    ----------
    archival_measures: pandas.DataFrame
        dataframe contains the archival data
    measurements: pandas.DataFrame or NoneType
        dataframe from VAST pipeline (nonetype for not plotting VAST)

    Returns:
    ----------
    fig
    '''
    fig = plt.figure(figsize=(12, 4), facecolor='w')
    ax1 = fig.add_subplot(121) # flux-time plot
    ax2 = fig.add_subplot(122) # flux-freq plot

    colors = {'AT20G':'C0','GLEAM':'C1','SUMSS':'C2','NVSS':'C3','TGSS':'C4','ASKAP':'C5','ASKAP_forced':'C6'}

    for i, row in archival_measures.iterrows():
        if row['type'] == 1:
            ax1.scatter(row['date'], row['flux'], color=colors[row['survey']])
            ax2.scatter(row['freq'], row['flux'], color=colors[row['survey']])
        else:
            ax1.errorbar(row['date'], row['flux'], yerr=row['flux']*0.3, color=colors[row['survey']], uplims=True, marker='o')
            ax2.errorbar(row['freq'], row['flux'], yerr=row['flux']*0.3, color=colors[row['survey']], uplims=True, marker='o')

    if isinstance(measurements, pd.DataFrame):
        for i, row in measurements.iterrows():
            if row['forced'] == True:
                ax1.errorbar(row['time'].year, 5*row['local_rms'], yerr=row['local_rms'], color=colors['ASKAP_forced'], uplims=True, marker='o')
                ax2.errorbar(0.887, 5*row['local_rms'], yerr=row['local_rms'], color=colors['ASKAP_forced'], uplims=True, marker='o') # need to change if there is mid-band
            else:
                ax1.scatter(row['time'].year, row['flux_peak'], color=colors['ASKAP'])
                ax2.scatter(0.887, row['flux_peak'], color=colors['ASKAP']) # need to change if there is mid-band

    ### add legend
    for survey in colors:
        ax2.scatter([], [], color=colors[survey], label=survey)

    ax1.set_yscale('log'); ax2.set_yscale('log')
    ax1.set_xlabel('year'); ax2.set_xlabel('frequency [GHz]')
    ax1.set_ylabel('flux [mJy/beam]'); ax2.set_ylabel('flux [mJy/beam]')
    ax2.legend()

    return fig

### Functions for overlay
def plot_VAST_overlay(ra, dec, imagepath, measurements, images, index=0):
    '''
    plot VAST contour over the image

    Params:
    ----------
    ra, dec: float
        coordinate of the source
    imagepath: str
        path for the fits file
    measurements: pandas.DataFrame
        dataframe contains all measurements for the source
    images: pandas.DataFrame
        dataframe contains all image information
    index: int, 0 by default
        index for the fits file

    Returns:
    ----------
    fig, ax
    '''
    with fits.open(imagepath) as hdulist:
        data = hdulist[index].data
        header = hdulist[index].header
        wcs = WCS(header).celestial

    fig = plt.figure(figsize=(5, 5), facecolor='w')
    ax = fig.add_subplot(111, projection=wcs)
    ax, im = plot_fits(data, ax)
    ax = source_crosshair((ra, dec), ax, color='red', sep=10, length=10)

    ### work out the highest selavy point
    selavy_measures = measurements[~measurements['forced']]
    maxflux_row = selavy_measures.iloc[selavy_measures['snr'].argmax()]
    imageid = maxflux_row['image_id']
    vastpath = images.loc[imageid]['path']

    try: # handle fits file doesnot exist
        fitsimage = FITSIMAGE(vastpath)
        cutout = fitsimage.cutout(ra, dec, 40)
    except:
        # set axislable invisible
        ax.coords[0].set_axislabel(' ')
        ax.coords[1].set_axislabel(' ')
        return fig, ax

    ax, ct = plot_contour(cutout.data, cutout.wcs, ax, colors='blue')

    # set axislable invisible
    ax.coords[0].set_axislabel(' ')
    ax.coords[1].set_axislabel(' ')
    return fig, ax

def plot_multicut(ra, dec, imagepath, index=0):
    '''
    plot VAST contour over the image

    Params:
    ----------
    ra, dec: float
        coordinate of the source
    imagepath: str
        path for the fits file
    measurements: pandas.DataFrame
        dataframe contains all measurements for the source
    images: pandas.DataFrame
        dataframe contains all image information
    index: int, 0 by default
        index for the fits file

    Returns:
    ----------
    fig, ax
    '''
    with fits.open(imagepath) as hdulist:
        data = hdulist[index].data
        header = hdulist[index].header
        wcs = WCS(header).celestial

    fig = plt.figure(figsize=(5, 5), facecolor='w')
    ax = fig.add_subplot(111, projection=wcs)
    ax, im = plot_fits(data, ax)
    ax = source_crosshair((ra, dec), ax, color='red', sep=10, length=10)

    return fig, ax



