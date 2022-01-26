# ztwang201605@gmail.com

# this code is used for parsing local atel webpages
# check astronomical sources within the webpage, update django database, write relavent files

### html parsing
from bs4 import BeautifulSoup
from numpy import rec
import requests
import re

### astropy packages
from astropy.coordinates import SkyCoord
import astropy.units as u

### package development
import os
import pkg_resources
import logging
import crawlerUtils

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

    def _extractMain_(self):
        '''extract main text from html'''
        self.logger.info('extracting atel main telegram from the html...')
        telegrams = self.htmlsoup.findAll('div', attrs={'id': 'telegram'})
        assert len(telegrams) == 1, 'Error, {} telegram(s) found in this telegram...'.format(len(telegrams))
        self.telegram = telegrams[0]

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

    


