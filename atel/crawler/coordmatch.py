### this is a separate script mainly focus on matching coordinates in html/txt
import logging
import re
from bs4 import BeautifulSoup
import requests

from astropy.coordinates import SkyCoord
import astropy.units as u

import time


class CoordMatch:
    '''
    Class for match coordinate in a given string
    '''
    def __init__(self, rawString, module='main'):
        self.rawString = rawString
        self.logger = logging.getLogger(f'{module}.coordmatch')
        self._buildAux_()

    def _buildAux_(self):
        '''
        building auxiliary re pattern within this function
        adapted from https://github.com/thespacedoctor/atelParser/blob/master/atelParser/mysql.py

        function supports different starts/middles/ends pattern
        P.S.: this looks like a good tool for interpreting regex: https://regexr.com/3f4vo
        '''
        ### format for R.A.
        #           #   RA or R.A.  #          #         (J2000.0)       #          #      date       #
        self.startRA = r'''(((R\.?\s*A\.?\s*|Coord\s*)(\(?\s*J2000(\.0)?\s*\)?\s*)?(=|:|\s)|(\d{4}-\d{2}-\d{2}))\s*)?'''
        self.endRA = r'''\s*(deg)?(d)?(and)?(,)?(;)?(o)?\s+'''
        # end with deg, d and and
        ### format for Decl.
        self.startDec = r'''\s*(Decl?\.?\s*(=|:)?)?\s*'''
        self.endDec = r'''(\s*deg)?(\s*d)?(\s*o)?'''

    def matchSex(self):
        '''
        match sexegesimal coordinates in the text

        Returns:
        ----------
        coordStrings: list
            each element is coordinate string parsed from `self.rawString` - units included
        coordRaws: list
            each element is the raw text for a given coordinate in coordStrings
        '''
        raSex = r"""(?P<raSex>(
                                (?P<raHrs>(\d)|([0-1]\d)|([2][0-3]))(\s*(:|h|\s)\s*)
                                (?P<raMin>[0-5][0-9])(\s*(:|m|\s)\s*)
                                (?P<raSec>[0-5]\d|\d(?!\d))s?(?P<raSubSec>\.\d{1,})?(\s|\s?s)?
                            )
                    )"""

        decSex = r"""(?P<decSex>(
                                (?P<decDeg>(\+|-|–)?[0-8]\d)(\s*(:|\s|d|deg|o)\s*)
                                (?P<decMin>[0-5][0-9])(\s*(:|\s|m|')\s*)
                                (?P<decSec>[0-5]?\d)'?\s?(?P<decSubSec>\.\d{1,})?'?"?s?
                            )
                    )"""

        reSex = r'''(?:(?<=^)|(?<=\s)|(?<=\>))(?P<coordSex> %s %s %s %s %s %s)''' % (self.startRA, raSex, self.endRA, self.startDec, decSex, self.endDec)
        reSexEx = re.compile(reSex, re.S | re.I | re.X)

        ### convert match group to coordinate ###
        coordStrings = []; coordRaws = []
        for reMatch in re.finditer(reSexEx, self.rawString):
            matchDict = reMatch.groupdict()
            coordStr = "{}h{}m{}{}s {}d{}m{}{}s".format(
                matchDict['raHrs'], matchDict['raMin'], matchDict['raSec'], matchDict['raSubSec'],
                matchDict['decDeg'], matchDict['decMin'], matchDict['decSec'], matchDict['decSubSec'],
            )
            coordStr = coordStr.replace('None', '').replace('–', '-') # replace None for raSubSec and decSubSec
            coordStrings.append(coordStr); coordRaws.append(matchDict['coordSex'])
        return coordStrings, coordRaws

    def matchDeg(self):
        '''
        match decimal coordinates in the text - we strongly suggest you do `self.matchSex` first and then do `self.matchDeg`

        Returns:
        ----------
        coordStrings: list
            each element is coordinate string parsed from `self.rawString` - units included
        coordRaws: list
            each element is the raw text for a given coordinate in coordStrings
        '''
        # TODO: remove Sex match first and then perform Deg match
        raDeg = r'''(?P<raDeg>(([0-3]\d{2})|(\d{2})|\d)(\.\d{3,}))'''
        decDeg = r'''(?P<decDeg>[\+\-\–]?([0-8]?\d(\.\d{3,})))''' # given the precision of the instruments, we assert that all decimal coordinates are with decimal places
        reDeg = r'''(?:(?<=^)|(?<=\s)|(?<=\>))(?P<coordDeg> %s %s %s %s %s %s)''' % (self.startRA, raDeg, self.endRA, self.startDec, decDeg, self.endDec)
        reDegEx = re.compile(reDeg, re.S | re.I | re.X)

        ### convert match group to coordinate ###
        coordStrings = []; coordRaws = []
        for reMatch in re.finditer(reDegEx, self.rawString):
            matchDict = reMatch.groupdict()
            coordStr = '{}d {}d'.format(matchDict['raDeg'], matchDict['decDeg'])
            coordStr.replace('–', '-')
            coordStrings.append(coordStr); coordRaws.append(matchDict['coordDeg'])
        return coordStrings, coordRaws


class TNSQuery:
    '''
    Literally we need to use TNS-API to do so. The `TNSQuery` class here uses requests module to achieve that (NOT RECOMMENDED) 
    '''
    def __init__(self, tnsname, module='main'):
        self.logger = logging.getLogger(f'{module}.tnsquery')
        self.tnsname = tnsname
        pass
        
    def _getresponse_(self):
        ### use requests to get response
        headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36"}
        tnsurl = f'https://www.wis-tns.org/object/{self.tnsname}'
        r = requests.get(tnsurl, headers=headers)
        self.logger.info(f'fetching data from TNS server... status code {r.status_code}')
        if r.status_code == 404: raise ValueError("The tns name is invalid or TNS server is not accessable for now...")
        assert r.status_code == 200, "cannot request data from tns for object {} - status_code {}".format(self.tnsname, r.status_code)
        ### convert response to BeautifulSoup Object
        self.soup = BeautifulSoup(r.text, 'html.parser')

    def _findcoord_(self):
        '''find coordinate from html'''
        _field = self.soup.findAll('div', attrs={'class':'field field-radec'})
        assert len(_field) > 0, f"no radec field found... Perhaps you get a wrong object name ({self.tnsname})?"
        return _field[0].findAll('div', attrs={'class':'value'})[0].text

    def checkcoord(self, sleeptime=60, max_attempt=3):
        '''
        check coordinate for a given TNS object

        Params:
        ----------
        sleeptime: int or float
            seconds for the system to get another attempt if failed
        max_attempt: int
            the maximum attempt for a single try
        
        Returns:
        ----------
        coord: string
            coordinate of the object - in the format of HH:MM:SS.SSS DD:MM:SS.SS
        '''
        ### try get the html page
        tnsrespond = False; attempt = 0
        while not tnsrespond:
            try:
                self._getresponse_()
                tnsrespond = True   # get response successfully
            except ValueError:
                self.logger.critical(f'The tns name {self.tnsname} is invalid or TNS server is not accessable for now...')
            except AssertionError:
                if attempt > max_attempt: raise ValueError(f'Reached maximum ({max_attempt}) attempts, Abort!')
                time.sleep(sleeptime) # sleep and attempt again
                attempt = attempt + 1

        ### parse coordinate
        return self._findcoord_()
