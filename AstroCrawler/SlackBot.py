
from slack import WebClient
import logging

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

    def postMessage(self, text):
        response = self.slack_client.chat_postMessage(
            channel = self.channel,
            text = text,
        )
        self.logger.info('postMessage - {}'.format(response['ok']))
        return response

    def replyMessage(self, text, message, reply_broadcast=False):
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
    def __init__():
        pass


### dev function ###
def loadDevSlack():
    localDict = locals()
    with open('../../slackbotTokens.py') as fp:
        tokenscript = fp.read()
    exec(tokenscript, globals(), localDict)
    return SlackBot(
        token = localDict["SLACKTOKEN"],
        channel = localDict["CHANNELNAME"],
        channel_id = localDict["CHANNELID"],
        bot_uid = localDict["BOTUID"],
    )

