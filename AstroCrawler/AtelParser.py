# ztwang201605@gmail.com

# this code is used for parsing local atel webpages
# check astronomical sources within the webpage, update django database, write relavent files

### html parsing
from bs4 import BeautifulSoup
from lark import logger
from numpy import rec
import requests
import re

### astropy packages
from astropy import coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u

### package development
import re
import os
import pkg_resources
import logging
import crawlerUtils
import coordMatch

class AtelHome:
    '''
    Class for parsing ATel home page

    Key Attributes:
    ----------
    maxcount: int
        saves the total Telegram count appearing on the home page
    
    Initiate Params:
    -----------
    htmlpath: str, atel home page by default
    url: boolean, True by default
        if the htmlpath is a url or a filepath
    '''
    def __init__(
        self, 
        htmlpath = 'https://www.astronomerstelegram.org/',
        url = True,
        ):
        self.logger = logging.getLogger('atelparser.atelhome')
        if not url:
            with open(htmlpath, encoding='utf-8') as fp:
                self.html = fp.read()
        else:
            headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36"}
            response = requests.get(htmlpath, headers=headers)
            assert response.status_code == 200, f'Unsuccessful attempt for requesting {url}...'
            self.html = response.text
        self._get_count() # get atel count

    def _get_count(self):
        reCount = re.compile('selected of ([0-9]*?) telegrams', re.IGNORECASE)
        CountFind = re.findall(reCount, self.html)
        assert len(CountFind) > 0, 'No pattern found for count in atel home...'

        self.maxcount = int(CountFind[0])


########## for atel page ############
class AtelPage:
    '''
    Class for parsing a single telegram page
    '''
    def __init__(
        self,
        pagepath,
        pageurl = False,
        datapath = None,
        searchMethod = 'main,tns,title',
        ):
        '''
        Initiate function for `AtelPage` class

        Params:
        ----------
        pagepath: str
            page url or page directory (if saved locally)
        pageurl: boolean, False by default
            whether `pagepath` is a url or not
        datapath: string or Nonetype
            if not Nonetype, save all relavent data in that place
        searchMethod: string
            coordinates matching methods in `findCoord`.
        '''
        self.logger = logging.getLogger('atelparser.atelpage') # set logging
        ### check and make directory ###
        if datapath is None:
            self.logger.critical('`datapath` cannot be Nonetype...')
            raise ValueError('None for datapath')
        else:
            crawlerUtils.makedirs(datapath)
        ######
        if pageurl: raise NotImplementedError('`pageurl = True` for Atelpage is under development...')
        with open(pagepath, encoding='utf-8') as fp:
            self.html = fp.read()
        self.logger.info(f'loading atel telegram from {pagepath}')
        self.htmlsoup = BeautifulSoup(self.html, 'html.parser')
        self.searchMethod = searchMethod

    def _extractMain_(self):
        '''extract main text from html'''
        self.logger.info('extracting atel main telegram from the html...')
        telegrams = self.htmlsoup.findAll('div', attrs={'id': 'telegram'})
        assert len(telegrams) == 1, 'Error, {} telegram(s) found in this telegram...'.format(len(telegrams))
        self.telegram = telegrams[0].__str__() # convert from `Tag` to `string`

    def _extractAtelID_(self):
        '''extract atel id for this page'''
        ### extract title first
        self.logger.info('extracting atel id from the html...')
        titles = self.htmlsoup.findAll('title')
        assert len(titles) > 0, 'No title found in the html...'
        if len(titles) > 1: self.logger.warning(f'{len(titles)} titles found in the html...')
        titleString = titles[0].string.strip()
        ### use regex to parse the id
        reAtelID = re.compile(r'ATel #(.*?):')
        atelids = re.findall(reAtelID, titleString)
        assert len(atelids) > 0, 'No atel ID found in the html...'
        if len(atelids) > 1: self.logger.warning(f'{len(atelids)} atel ids found in the html...')
        self.atelid = atelids[0]
        ### save atel title to an attribute
        self.ateltitle = ':'.join(titleString.split(':')[1:]).strip() # to avoid titles like xxx: xxxxxx

    def _extractSubjects_(self):
        '''extract subject from the atel'''
        reSubjects = re.compile('<p class="subjects">Subjects:(.*?)<\/p>', re.I)
        subjectStr = reSubjects.findall(self.telegram)
        assert len(subjectStr) > 0, 'No atel subjects found in the html...'
        subjectRaw = subjectStr[0].strip().lower().split(',')
        self.subjects = [subject.strip() for subject in subjectRaw]

    def _findTNS_(self):
        '''find TNS source within the telegram'''
        reTNSurl = r'\"https://wis-tns.weizmann.ac.il/object/(?P<tnsObj>[0-9a-zA-Z]{5,9})\"'
        reTNSEx = re.compile(reTNSurl, re.S | re.X)
        
        tnsObjs = []
        for reMatch in re.finditer(reTNSEx, self.telegram):
            tnsObj = reMatch.groupdict()['tnsObj']
            if tnsObj not in tnsObjs: tnsObjs.append(tnsObj)

        self.tnsObjs = tnsObjs

    def _colorMatch_(self, coordRaws, color='powderblue'):
        ''' replace the matching strings to a different format '''
        self.coloredTelegram = self.telegram
        for coordRaw in coordRaws:
            formatRaw = r'''<font style="background-color:{};"> {} </font>'''.format(color, coordRaw)
            self.coloredTelegram = self.telegram.replace(coordRaw, formatRaw)

    def findCoord(self):
        '''
        search for coordinates in the atel
        ''' 
        coords = []
        ### search for coordinates
        if 'main' in self.searchMethod:
            coordmatcher = coordMatch.CoordMatch(self.telegram, module='atelparser')
            textCoords, coordRaws = coordmatcher.matchAll()
            self._colorMatch_(coordRaws, color='powderblue') # replace the matched string with formatting
            if len(textCoords) > 0: coords.append(SkyCoord(textCoords))
            self.logger.info(f'{len(textCoords)} sources found in the main text...')

        ### search for TNS coordinates
        if 'tns' in self.searchMethod:
            self._findTNS_(); tnsCoords = []
            for tnsName in self.tnsObjs:
                tnsquery = coordMatch.TNSQuery(tnsName, module='atelparser')
                try: tnsCoords.append(tnsquery.checkcoord())
                except: self.logger.warning(f'{tnsName} is not existed in the transient name server...')
            if len(tnsCoords) > 0: coords.append(SkyCoord(tnsCoords, unit=(u.hourangle, u.degree)))
            self.logger.info(f'{len(tnsCoords)} sources found in the TNS server...')
        
        ### search for sources in coordinate
        if 'title' in self.searchMethod:
            titlesearcher = coordMatch.TitleSearch(self.ateltitle, module='atelparser')
            try: titleCoords = titlesearcher.simbadsearch()
            except: self.logger.warning('Simbad Title search failed... Likely Connection Error...'); titleCoords = []
            if len(titleCoords) > 0: coords.append(SkyCoord(titleCoords, unit=(u.hourangle, u.degree)))
            self.logger.info(f'{len(titleCoords)} sources found in the title...')
        
        ### convert these strings to astropy.coordinates.SkyCoord objects
        if len(coords) == 0: self.coords = []; return
        ### use coordfilter to clean the coordinate
        coordfilter = coordMatch.CoordFilter(coords, module='atelparser')
        coordfilter.filterSources(radius = 5., method=None)
        self.coords = coordfilter.uniqueSources; return 



    




