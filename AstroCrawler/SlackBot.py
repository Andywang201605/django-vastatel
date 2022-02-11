
from slack import WebClient
import pandas as pd
import logging
import json
import os

from .config import *
from . import AtelDatabase

from astropy.coordinates import SkyCoord
import astropy.units as u

class SlackBot:
    '''
    vast-atelbot for sending messages and upload files
    '''
    def __init__(self, token, channel, channel_id, bot_uid, module='main'):
        ### read token
        self.logger = logging.getLogger(f'{module}.slackbot')
        self._init_client(token)
        self.channel = channel
        self.channel_id = channel_id
        self.bot_uid = bot_uid

    def _init_client(self, token):
        self.slack_client = WebClient(token)

    def postMessage(self, text, unfurl_links=False):
        response = self.slack_client.chat_postMessage(
            channel = self.channel,
            text = text,
            unfurl_links = unfurl_links,
        )
        self.logger.info('postMessage - {}'.format(response['ok']))
        return response

    def replyMessage(self, text, message, reply_broadcast=False, unfurl_links = False):
        '''
        Reply to a message

        Params:
        ----------
        text: str
            text/message you want the robot to reply
        message: dict
            json format message, get directly from api.WebClient.conversations_history function
        reply_broadcast: bool, True by default
            whether to post in the main channel too
        '''
        thread_ts = message['ts']
        response = self.slack_client.chat_postMessage(
            channel = self.channel_id,
            text = text,
            thread_ts = thread_ts,
            reply_broadcast = reply_broadcast,
            unfurl_links = unfurl_links,
        )
        self.logger.info('replyMessage - {}'.format(response['ok']))
        return response

    def replyFile(self, file, comment, message):
        thread_ts = message['ts']
        response = self.slack_client.files_upload(
            channels = self.channel,
            file = file,
            initial_comment = comment,
            thread_ts = thread_ts,
        )
        self.logger.info('replyFile - {} - {}'.format(file, response['ok']))
        return response

    def uploadFile(self, file, comment):
        response = self.slack_client.files_upload(
            channels = self.channel,
            file = file,
            initial_comment = comment,
        )
        self.logger.info('uploadFile - {} - {}'.format(file, response['ok']))
        return response

    def checkHistory(self, tsstart, limit=200):
        response = self.slack_client.conversations_history(
            channel = self.channel_id,
            oldest = '{:.6f}'.format(tsstart),
            limit = limit,
        )

        self.logger.debug('{} messages found from {}'.format(len(response['messages']), tsstart))
        return response['messages']

################################################################################
### ATelBot ###

class ATelVASTBot(SlackBot):
    '''
    object inherited from SlackBot
    ATelVASTBot is used for posting ATel information to Slack
    '''
    ### send information ###
    def _loadAtel_(self, atelid):
        with open(f'{ATELDIR}/ATEL{atelid}/atelinfo.json') as fp:
            atelinfo = json.load(fp)
        return atelinfo
        
    def _postNewAtel_(self, atelid, manual=False, dryrun=False):
        '''Post Atel in the main conversation'''
        atelinfo = self._loadAtel_(atelid)
        ### start a message ###
        messageSubject = ','.join(atelinfo['subject'])
        atelLink = f'https://www.astronomerstelegram.org/?read={atelid}'
        numSrc = len(atelinfo['source'])
        # message content #
        if manual: line1 = f'''*[Update ATel {atelid}]* - {numSrc} sources parsed'''
        else: line1 = f'''*[ATel {atelid}]* - {numSrc} sources parsed'''
        if dryrun: 
            print(f'''{line1}\n*Title:* {atelinfo['title']}\n*Subject:* {messageSubject}\n*Check original telegram* <{atelLink}|here>!''')
            return
        atelmessage = f'''{line1}\n*Title:* {atelinfo['title']}\n*Subject:* {messageSubject}\n*Check original telegram* <{atelLink}|here>!'''
        response = self.postMessage(atelmessage, unfurl_links=False); return response

    ### for replying a source ###
    def _findVASTLink_(self, taskname):
        vastSrc = pd.read_json(f'{LAMBDADIR}/{taskname}/VASTrun.source.json')
        if len(vastSrc) == 0: return '*PipeRun:* No VAST Source within 1 arcmin'
        nearSrc = vastSrc.iloc[0]; srcID = nearSrc.name; separation = nearSrc['separation']
        srcLink = f'https://dev.pipeline.vast-survey.org/sources/{srcID}'
        return f'''*PipeRun:* VAST source ({separation:.2f} arcsec away): <{srcLink}|{srcID}>\ncheck *PDF* for more information'''

    def _checkVASTmeasure_(self, taskname):
        meaPath = f'{LAMBDADIR}/{taskname}/source_measurement.pickle'
        if not os.path.exists(meaPath): return '*VAST/RACS footprint*: no'
        vastMea = pd.read_pickle(meaPath)
        if len(vastMea) == 0: return '*VAST/RACS footprint*: no'
        if len(vastMea) == 1: return '*VAST/RACS footprint*: RACS'
        return '*VAST/RACS footprint*: VAST/RACS'
        

    def _replySource_(self, response, atelid, replyFormatHTML = True, dryrun=False):
        ateldbRecord = AtelDatabase.AtelDataBase(atelid)
        if response is None: dryrun = True # if there is no response - run a dry run
        if 'comet' in ateldbRecord.subject:
            if dryrun: print('This source is probably a comet :comet:')
            else:
                self.replyMessage(
                    'This source is probably a comet :comet:',
                    response['message'],
                )
            return

        for _, coord, taskname in ateldbRecord.srctask:
            ## source information related
            srcname = self._longName_(coord)
            footLine = self._checkVASTmeasure_(taskname)
            linkLine = self._findVASTLink_(taskname)
            ## file related
            pdfPath = f'{ATELDIR}/ATEL{atelid}/multiweb_{taskname}.pdf'
            if dryrun:
                print(f"*Source:* {srcname}\n{footLine}\n{linkLine}")
                print(f'Send File - {pdfPath}')
            else:
                self.replyFile(
                    file = pdfPath,
                    comment = f"*Source:* {srcname}\n{footLine}\n{linkLine}",
                    message = response['message'],
                )
        ### post parsed atel to the thread ###
        if replyFormatHTML:
            if dryrun: 
                print('You can refer to a formatted telegram here (objects parsed are highlighted)')
                print(f'Send File - {ATELDIR}/ATEL{atelid}/coloredATel.png')
            else:
                self.replyFile(
                    file = f'{ATELDIR}/ATEL{atelid}/coloredATel.png',
                    comment = 'You can refer to a formatted telegram here (objects parsed are highlighted)',
                    message = response['message']
                )

    def postAtel(self, atelid, manual=False, dryrun=False):
        response = self._postNewAtel_(atelid, manual=manual, dryrun=dryrun)
        self._replySource_(response=response, atelid=atelid, dryrun=dryrun)

        


    ### utiles
    def _shortName_(self, coord):
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

    def _longName_(self, coord):
        sign = '+' if coord[1] > 0 else '-' # check sign of decl first
        coord = SkyCoord(*coord, unit=u.degree)
        namera = '{:0>2.0f}{:0>2.0f}{:0>4.1f}'.format(coord.ra.hms.h, coord.ra.hms.m, coord.ra.hms.s)
        namede = '{:0>2.0f}{:0>2.0f}{:0>4.1f}'.format(abs(coord.dec.dms.d), abs(coord.dec.dms.m), abs(coord.dec.dms.s))
        return f'J{namera}{sign}{namede}'



