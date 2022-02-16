#ztwang201605@gmail.com

from .image_data import *
from .download import *

from vasttools.query import Query

import glob

class PipelineSource:
    '''
    Run all analysis based on pipeline run result
    '''
    def __init__(self, coord, measurements, images, sourcepath, download_list):
        '''
        Initiate Function and variable checking for PipelineSource Object

        Params:
        ----------
        coord: list, tuple or SkyCoord
            coordinate of the source
        measurements: pandas.DataFrame
            measurement for the source from the pipeline
        images: pandas.DataFrame
            dataframe contains all images detail
        sourcepath: str
            path for storing all stuff for the source
        '''
        ### source position
        if isinstance(coord, SkyCoord):
            ra = coord.ra.value; dec = coord.dec.value
        if isinstance(coord, tuple) or isinstance(coord, list):
            ra, dec = coord
        self.ra = ra; self.dec = dec

        ### pipeline products
        self.measurements = measurements
        self.images = images

        ### directory
        # base path
        self.sourcepath = sourcepath
        self._makepath(self.sourcepath)
        # image path
        self.imagepath = os.path.join(self.sourcepath, 'img/')
        self._makepath(self.imagepath)
        # download list
        self.download_list = download_list

    def _makepath(self, path):
        '''
        check if the path exists, if not make a new directory

        Params:
        ----------
        path: str
        '''
        if not os.path.exists(path):
            os.makedirs(path)

    def _savefig(self, fig, figpath, **kwargs):
        '''
        Save figure to figpath, with **kwargs parameters passed to fig.savefig

        Params:
        ----------
        fig: matplotlib.pyplot.figure
        figpath: str
        '''
        kwargs.setdefault('bbox_inches', 'tight')
        fig.savefig(figpath, **kwargs)

        plt.clf()
        plt.close(fig)

    def plotVASTStokesI(self, radius=300., imagepath_col='path', samescale=True):
        '''
        Plot multiepoch StokesI cutout for the source

        Params:
        ----------
        radius: float or int
            cutout radius
        imagepath_col: str
            column name for imagepath (StokesI)
        samescale: bool
            If use the same scale and limit for all subplots
        '''
        fig = plot_multiepoch_cutout(
            self.ra, self.dec,
            self.measurements, self.images,
            radius = radius, imagepath_col=imagepath_col,
            samescale=samescale,
        )
        if samescale == True:
            self._savefig(fig, os.path.join(self.imagepath, 'StokesI_{}.jpg'.format(int(radius))))
        else:
            self._savefig(fig, os.path.join(self.imagepath, 'StokesI_{}_scale.jpg'.format(int(radius))))

    def plotVASTStokesV(self, radius=300., imagepath_col='Vpath', samescale=True):
        '''
        Plot multiepoch StokesV cutout for the source

        Params:
        ----------
        radius: float or int
            cutout radius
        imagepath_col: str
            column name for imagepath (StokesV)
        samescale: bool
            If use the same scale and limit for all subplots
        '''
        assert imagepath_col in self.images, f'Column {imagepath_col} not exist in Images DataFrame!'

        fig = plot_multiepoch_cutout(
            self.ra, self.dec,
            self.measurements, self.images,
            radius = radius, imagepath_col=imagepath_col,
            samescale=samescale,
        )
        self._savefig(fig, os.path.join(self.imagepath, 'StokesV_{}.jpg'.format(int(radius))))

    def download_archival(self):
        '''
        Download all archival data based on ./setups/multiwavelength_information.json
        '''

        download_archival_multithreading(
            self.ra, self.dec,
            self.download_list,
            self.imagepath,
            maxthreads=32
        )

        clear_download_cache()

    def fetch_archival_data(self):
        '''
        Download archival data and save it to a .dat file
        '''
        archivalcatalog_path = pkg_resources.resource_filename(
            __name__, './setup/archival_catalog.json'
        )
        with open(archivalcatalog_path) as fp:
            archivalcatalog = json.load(fp)

        archivaldata = get_archival_data(
            (self.ra, self.dec),
            archivalcatalog,
            fitspath=self.imagepath
        )

        with open(os.path.join(self.sourcepath, 'archival_flux.dat'), 'w') as fp:
            fp.writelines(archivaldata)

    def plot_archival_lightcurve(self):
        '''
        Plot lightcurve for archival data
        '''
        if not os.path.exists(os.path.join(self.sourcepath, 'archival_flux.dat')):
            return
        archival_measures = pd.read_csv(os.path.join(self.sourcepath, 'archival_flux.dat'))
        fig = plot_archival_lightcurve(archival_measures, self.measurements)
        self._savefig(fig, os.path.join(self.imagepath, 'archival_lightcurve.png'))

    def _getVLASS_download(self, radius):
        '''
        get epochs in VLASS that have been downloaded
        '''
        fileFormat = f'{self.imagepath}/VLASS?.?_{int(radius)}.fits'
        return glob.glob(fileFormat)

    def plot_multiwavelength_overlay(self):
        '''
        plot multiwavelength results overlaid with radio countour
        '''
        for survey, radius in self.download_list:
            if survey == 'VLASS':
                fitspaths = self._getVLASS_download(radius)
            else:
                survey = survey.replace(' ', '_')
                fitspaths = [os.path.join(self.imagepath, f'{survey}_{radius}.fits')]

            for fitspath in fitspaths:
                try:
                    fig, ax = plot_VAST_overlay(
                        self.ra, self.dec,
                        fitspath, self.measurements,
                        self.images
                    )
                except:
                    continue

                self._savefig(fig, fitspath.replace('.fits', '.png'))

    def plot_wise_cc(self):
        '''
        Plot wise color-color plot
        '''
        fig, ax = plot_wise_cc((self.ra, self.dec))
        self._savefig(fig, os.path.join(self.imagepath, 'wise-cc.png'))

    def plotVASTlightcurve(self):
        '''
        Plot lightcurve from VAST only
        '''
        fig, ax = plot_VAST_lightcurve(self.measurements)
        self._savefig(fig, os.path.join(self.imagepath, 'VASTlightcurve.png'))

    # def makewebpage(self):
    #     pipeweb = PipelineWeb(
    #         (self.ra, self.dec),
    #         'source_web.html',
    #         self.sourcepath
    #     )
    #     pipeweb.makefullweb()

    def sourceAnalysis(self):
        ### Plot VAST cutouts
        self.plotVASTStokesI(radius=300., samescale=True)
        self.plotVASTStokesI(radius=300.,samescale=False)
        self.plotVASTStokesI(radius=600., samescale=True)
        self.plotVASTStokesV(radius=300.)

        ### Plot VAST lightcurves
        self.plotVASTlightcurve()

        ### Download archival
        self.download_archival()
        self.plot_multiwavelength_overlay()
        self.plot_wise_cc()
        
        self.fetch_archival_data()
        self.plot_archival_lightcurve()
        
        ###
        #self.makewebpage()

def _getStokesVpath(StokesIpath):
        return StokesIpath.replace('STOKESI', 'STOKESV').replace('.I.', '.V.')

class VASTSource:
    '''
    Use vasttool to get source measurements
    '''
    def __init__(self, coord, sourcepath, download_list, pickleoverwrite=False, pilotbasefolder='/import/ada1/askap/PILOT/release/', ncpu=1):
        '''
        Initiate Function for VASTSource object

        Params:
        ----------
        coord: list, tuple or SkyCoord
        sourcepath: str
            The folder where all data put
        pilotbasefolder: str
            The base folder for pilot survey data, used in vasttools
        '''
        ### source position
        if isinstance(coord, SkyCoord):
            ra = coord.ra.value; dec = coord.dec.value
        if isinstance(coord, tuple) or isinstance(coord, list):
            ra, dec = coord
        self.ra = ra; self.dec = dec

        self.ncpu = ncpu
        self.pilotbasefolder = pilotbasefolder
        self.download_list = download_list

        ### directory
        # base path
        self.sourcepath = sourcepath
        self._makepath(self.sourcepath)
        # image path
        self.imagepath = os.path.join(self.sourcepath, 'img/')
        self._makepath(self.imagepath)

        ### add measurements
        try:
            self._vasttoolquery(pickleoverwrite=pickleoverwrite)
        except:
            print('No source with vasttools...')

    def _makepath(self, path):
        '''
        check if the path exists, if not make a new directory

        Params:
        ----------
        path: str
        '''
        if not os.path.exists(path):
            os.makedirs(path)

    def _vasttoolquery(self, pickleoverwrite=True):
        '''
        Query measurements with the help from vasttools
        '''
        measurement_picklepath = os.path.join(self.sourcepath, 'source_measurement.pickle')
        if os.path.exists(measurement_picklepath) and pickleoverwrite == False:
            self.measurements = pd.read_pickle(measurement_picklepath)
            return

        coord = SkyCoord([self.ra], [self.dec], unit=u.deg)
        # TODO: vasttools doesn't support epoch12, 13 and so on, replaced with 0-11x 
        # TODO: vasttools will raise errors if ncpu > 1.
        query = Query(
            coords = coord,
            crossmatch_radius = 10.,
            epochs = '0,1,2,3x,4x,5x,6x,7x,8,9,10x,11x',
            forced_fits = True,
            ncpu=self.ncpu,
            base_folder=self.pilotbasefolder
        )
        query.find_sources()
        self.measurements = query.results.iloc[0].measurements
        self._formatmeasurement()

        self.measurements.to_pickle(measurement_picklepath)

    def _formatmeasurement(self):
        '''
        Format self.measurements - add a column for StokesV image, forced, localrms
        '''
        self.measurements['path'] = self.measurements['image']
        self.measurements['Vpath'] = self.measurements['path'].apply(_getStokesVpath)
        self.measurements['forced'] = ~self.measurements['detection']
        self.measurements['local_rms'] = self.measurements['rms_image']
        self.measurements['time'] = self.measurements['dateobs']
        self.measurements['image_id'] = self.measurements.index

        self.measurements['flux_peak'][self.measurements['forced']] = self.measurements['f_flux_peak'][self.measurements['forced']].to_list()
        self.measurements['snr'] = self.measurements['flux_peak'] / self.measurements['local_rms']

    def _savefig(self, fig, figpath, **kwargs):
        '''
        Save figure to figpath, with **kwargs parameters passed to fig.savefig

        Params:
        ----------
        fig: matplotlib.pyplot.figure
        figpath: str
        '''
        kwargs.setdefault('bbox_inches', 'tight')
        fig.savefig(figpath, **kwargs)

        plt.clf()
        plt.close(fig)

    def plotVASTStokesI(self, radius=300., imagepath_col='path', samescale=True):
        '''
        Plot multiepoch StokesI cutout for the source

        Params:
        ----------
        radius: float or int
            cutout radius
        imagepath_col: str
            column name for imagepath (StokesI)
        samescale: bool
            If use the same scale and limit for all subplots
        '''
        fig = plot_multiepoch_cutout(
            self.ra, self.dec,
            self.measurements, self.measurements,
            radius = radius, imagepath_col=imagepath_col,
            samescale=samescale,
        )
        if samescale == True:
            self._savefig(fig, os.path.join(self.imagepath, 'StokesI_{}.jpg'.format(int(radius))))
        else:
            self._savefig(fig, os.path.join(self.imagepath, 'StokesI_{}_scale.jpg'.format(int(radius))))

    def plotVASTStokesV(self, radius=300., imagepath_col='Vpath', samescale=True):
        '''
        Plot multiepoch StokesV cutout for the source

        Params:
        ----------
        radius: float or int
            cutout radius
        imagepath_col: str
            column name for imagepath (StokesV)
        samescale: bool
            If use the same scale and limit for all subplots
        '''
        assert imagepath_col in self.measurements, f'Column {imagepath_col} not exist in Images DataFrame!'

        fig = plot_multiepoch_cutout(
            self.ra, self.dec,
            self.measurements, self.measurements,
            radius = radius, imagepath_col=imagepath_col,
            samescale=samescale,
        )
        self._savefig(fig, os.path.join(self.imagepath, 'StokesV_{}.jpg'.format(int(radius))))

    def download_archival(self):
        '''
        Download all archival data based on self.download_list
        '''
        # load setup files
        # print(self.download_list)

        download_archival_multithreading(
            self.ra, self.dec,
            self.download_list,
            self.imagepath,
            maxthreads=16,
        )

        clear_download_cache()

    def fetch_archival_data(self):
        '''
        Download archival data and save it to a .dat file
        '''
        archivalcatalog_path = pkg_resources.resource_filename(
            __name__, './setup/archival_catalog.json'
        )
        with open(archivalcatalog_path) as fp:
            archivalcatalog = json.load(fp)

        archivaldata = get_archival_data(
            (self.ra, self.dec),
            archivalcatalog,
            fitspath=self.imagepath
        )

        with open(os.path.join(self.sourcepath, 'archival_flux.dat'), 'w') as fp:
            fp.writelines(archivaldata)

    def plot_archival_lightcurve(self):
        '''
        Plot lightcurve for archival data
        '''
        if not os.path.exists(os.path.join(self.sourcepath, 'archival_flux.dat')):
            return
        archival_measures = pd.read_csv(os.path.join(self.sourcepath, 'archival_flux.dat'))
        fig = plot_archival_lightcurve(archival_measures, self.measurements)
        self._savefig(fig, os.path.join(self.imagepath, 'archival_lightcurve.png'))

    def plot_multiwavelength_overlay(self):
        '''
        plot multiwavelength results overlaid with radio countour
        '''
        # load setup files
        
        for survey, radius in self.download_list:
            survey = survey.replace(' ', '_')
            fitspath = os.path.join(
                self.imagepath,
                f'{survey}_{int(radius)}.fits'
            )

            try:
                fig, ax = plot_VAST_overlay(
                    self.ra, self.dec,
                    fitspath, self.measurements,
                    self.measurements
                )
            except:
                try:
                    fig, ax = plot_multicut(
                        self.ra, self.dec,
                        fitspath,
                    )
                except:
                    continue

            self._savefig(fig, os.path.join(self.imagepath, f'{survey}_{int(radius)}.png'))

    def plot_wise_cc(self):
        '''
        Plot wise color-color plot
        '''
        fig, ax = plot_wise_cc((self.ra, self.dec))
        self._savefig(fig, os.path.join(self.imagepath, 'wise-cc.png'))

    def plotVASTlightcurve(self):
        '''
        Plot lightcurve from VAST only
        '''
        fig, ax = plot_VAST_lightcurve(self.measurements)
        self._savefig(fig, os.path.join(self.imagepath, 'VASTlightcurve.png'))

    def _measurementexists_(self):
        if os.path.exists(os.path.join(self.sourcepath, 'source_measurement.pickle')):
            return True
        return False

    def sourceAnalysis(self):
        ### Plot VAST cutouts
        if not self._measurementexists_():
            self.download_archival()
            self.plot_wise_cc()
            self.plot_multiwavelength_overlay()

            return

        self.plotVASTStokesI(radius=300., samescale=True)
        self.plotVASTStokesI(radius=300.,samescale=False)
        self.plotVASTStokesI(radius=600., samescale=True)
        self.plotVASTStokesV(radius=300.)

        ### Plot VAST lightcurves
        self.plotVASTlightcurve()

        ### Download archival
        self.download_archival()
        self.fetch_archival_data()
        self.plot_archival_lightcurve()
        self.plot_wise_cc()

        self.plot_multiwavelength_overlay()

        return

        ###
        # self.makewebpage()

        
