### main script for monitoring/running atel parsing ###

from operator import mod
import os
import logging
import shutil
import time

from . import AtelParser, crawlerUtils, AtelDatabase, analysisUtils, SlackBot
from .config import *

### execute bot configuration file ###
with open(BOTCONFIG) as fp:
    botConfigCode = fp.read()
exec(botConfigCode)

### running all analysis for a given atel

class AtelRun:
    '''
    Analysis a given atel

    Usage:
    ----------
    >>> run = AtelRun(atelid=10000)
    >>> run.RunAll()

    Note: you need to either provide `ATELDIR` in `config.py`
    or give `atelpath` when initiating AtelRun
    '''
    def __init__(
        self,
        atelid,
        atelpath = None,
        ):
        self.logger = logging.getLogger('atelparser.run')
        ### initialize ###
        if atelpath is None: atelpath = f'{ATELDIR}/ATEL{atelid}'
        crawlerUtils.makedirs(atelpath) # make directory
        self.atelid = int(atelid); self.atelpath = atelpath
        


    def _downloadATel_(self, maxAttempt=5, sleepTime=60, overwrite=False):
        htmlfile = f'{self.atelpath}/ATEL{self.atelid}.html'
        if os.path.exists(htmlfile) and not overwrite: 
            self.logger.info(f'{htmlfile} already exists... aborted!')
            return
        self.logger.info(f'downloading atel page from `https://astronomerstelegram.org/?read={self.atelid}`')
        crawlerUtils.downloadURL(
            url = f'https://astronomerstelegram.org/?read={self.atelid}',
            filepath = htmlfile,
            maxAttempt = maxAttempt, sleepTime = sleepTime
        )

    def _parseATel_(self):
        '''parse ATel and run'''
        self.logger.info('starting parsing information from atel page...')
        atelpage = AtelParser.AtelPage(
            pagepath = f'{self.atelpath}/ATEL{self.atelid}.html',
            datapath = self.atelpath
        )
        self.logger.info('finding coordinates...')
        atelpage.findCoord()
        self.logger.info('writing information to local files...')
        atelpage.dumpinfo()
        self.logger.info('updating django database...')
        atelpage.updateDatabase()

    def _runVASTTransient_(self):
        ateldbRecord = AtelDatabase.AtelDataBase(self.atelid)
        ### donot run if it is a comet ###
        if 'comet' in ateldbRecord.subject: return 
        for srcname, coord, taskname in ateldbRecord.srctask:
            self.logger.info(f'analysis source {srcname}...')
            analysisUtils._runVASTTransient_(coord, taskname)
            ### fetch pdf files from LAMBDADIR
            shutil.copy(f'{LAMBDADIR}/{taskname}/multiweb_{taskname}.pdf', f'{self.atelpath}/')


    def RunParse(self):
        ### download atel ###
        self._downloadATel_(overwrite=False)
        ### parse atel ###
        self._parseATel_()
        ### running VAST match analysis ###

    def RunVAST(self):
        ### run VAST analysis on the source
        self._runVASTTransient_()

    def RunAll(self):
        self.RunParse(); self.RunVAST()

class AtelMonitor:
    '''
    xxxxx
    '''
    def __init__(
        self,
        sleepTime = 600,
        ):
        ### initiate parameters ###
        self.sleepTime = sleepTime
        self.slacktoken = SLACKTOKEN
        self.botuid = BOTUID 
        self.channelname = CHANNELNAME 
        self.channelid = CHANNELID 
        self.errorchannel = ERRORCHANNEL 
        self.errorchannelid = ERRORCHANNELID 
        self.notifyuid = ERRORUID 
        ### initiate other things ###
        self.logger = logging.getLogger('atelparser.monitor')
        self._connectSlack_()

    def _connectSlack_(self):
        self.postBot = SlackBot.ATelVASTBot(
            token = self.slacktoken,
            channel = self.channelname,
            channel_id = self.channelid,
            bot_uid = self.botuid,
            module = 'atelparser'
        )
        self.errorBot = SlackBot.ATelVASTBot(
            token = self.slacktoken,
            channel = self.errorchannel,
            channel_id = self.errorchannelid,
            bot_uid = self.botuid,
            module = 'atelparser'
        )

    ### loop content ###
    def _checkATelNew_(self):
        home = AtelParser.AtelHome()
        return home.maxcount

    def _getRunIDs_(self):
        maxATel = self._checkATelNew_()
        maxDB = AtelDatabase._getLastest_()
        self.logger.info(f'Latest ATel online is {maxATel}, latest ATel in our database is {maxDB}')
        return list(range(maxDB+1, maxATel+1))

    def _runSingleATel_(self, atelid):
        self.logger.info(f'Parse and Analyse sources in ATel{atelid}')
        run = AtelRun(atelid)
        run.RunAll()
    ##########################

    def _sendStatus_(self):
        self.errorBot.postMessage(
            f'''Hi! <@{self.botuid}> is still running!'''
        )

    def _sendError_(self, message):
        self.errorBot.postMessage(
            f'''*CRITICAL* Bot has met an error <@{self.notifyuid}> - \n {message}'''
        )

    ##########################
    def run(self):
        ### setting logger ###
        crawlerUtils.setLogger(
            module = 'atelparser',
            streamlevel = logging.INFO,
        )
        running = True
        while running:
            self._sendStatus_()
            ######################
            try: runList = self._getRunIDs_() # atelids for further running...
            except Exception as error: 
                self._sendError_(f'{error} in `_getRunIDs_`')
                self.logger.error(f'{error} in `_getRunIDs_`')
            ######################
            if len(runList) != 0:
                self.logger.info(f'number of new atels to be updated: {len(runList)}')
                for atelid in runList:
                    self.logger.info(f'sleeping for 20 seconds to download atel{atelid} page')
                    time.sleep(20) # stop for 20 seconds for next one
                    ### download and parse atel ###
                    try:
                        run = AtelRun(atelid)
                        run.RunAll()
                    except Exception as error:
                        self._sendError_(f'{error} in `AtelRun` where atelid = {atelid}')
                        self.logger.error(f'{error} in `AtelRun` where atelid = {atelid}')
                    ### send information ###
                    try:
                        self.postBot.postAtel(atelid)
                    except Exception as error:
                        self._sendError_(f'{error} in `postAtel` where atelid = {atelid}')
                        self.logger.error(f'{error} in `postAtel` where atelid = {atelid}')
            self.logger.info(f'waiting for {self.sleepTime} seconds to check for new atels...')
            time.sleep(self.sleepTime)
            