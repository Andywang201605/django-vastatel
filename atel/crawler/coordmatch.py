### this is a separate script mainly focus on matching coordinates in html/txt
import logging
import re
import logging

from astropy.coordinates import SkyCoord
import astropy.units as u


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
        self.endRA = r'''\s*(deg)?(d)?(and)?(,)?(;)?(o)?\s*'''
        # end with deg, d and and
        ### format for Decl.
        self.startDec = r'''\s*(Decl?\.?\s*(=|:)?)?\s*'''
        self.endDec = r'''\s*(deg)?(d)?(o)?'''

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

        reSex = r'''(?P<coordSex> %s %s %s %s %s %s)''' % (self.startRA, raSex, self.endRA, self.startDec, decSex, self.endDec)
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
        raDeg = r'''(?P<raDeg>(([0-3]\d{2})|(\d{2})|\d)(\.\d{1,})?)'''
        decDeg = r'''(?P<decDeg>[\+\-\–]?([0-8]?\d(\.\d{1,})?))'''
        reDeg = r'''(?P<coordDeg> %s %s %s %s %s %s)''' % (self.startRA, raDeg, self.endRA, self.startDec, decDeg, self.endDec)
        reDegEx = re.compile(reDeg, re.S | re.I | re.X)

        ### convert match group to coordinate ###
        coordStrings = []; coordRaws = []
        for reMatch in re.finditer(reDegEx, self.rawString):
            matchDict = reMatch.groupdict()
            coordStr = '{}d {}d'.format(matchDict['raDeg'], matchDict['decDeg'])
            coordStr.replace('–', '-')
            coordStrings.append(coordStr); coordRaws.append(matchDict['coordDeg'])
        return coordStrings, coordRaws

