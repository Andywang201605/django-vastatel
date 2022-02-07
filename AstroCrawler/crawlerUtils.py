import logging

def setLogger(
    module = 'main',
    loglevel = logging.DEBUG,
    streamlevel = logging.DEBUG,
    logfname = None,
    filelevel = logging.INFO,
    ):
    '''
    setting a `logging.Logger` instance for a given `model`

    Params:
    ----------
    model: str
        your main model name, 'main' by default
    loglevel: logging.LEVEL instance
        set the whole level for the logger
    streamlevel, filelevel: logging.LEVEL instance
        set level for the logger with regard to stream and file respectively
    logfname: str or Nonetype
        the log file name for your logger, for stream logger only use None
    '''
    logger = logging.getLogger(module)
    logger.setLevel(loglevel)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # stream logger handler
    console = logging.StreamHandler()
    console.setLevel(streamlevel)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # file logger handler
    if logfname is not None:
        file = logging.FileHandler(logfname, 'w')
        file.setLevel(filelevel)
        file.setFormatter(formatter)
        logger.addHandler(file)

import os
### make new directory
def makedirs(path):
    '''
    making new directory if the path is not existed
    '''
    if not os.path.exists(path):
        os.makedirs(path)

### response downloading related ###
import requests
import time

def _requestURL_(url):
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36"}
    response = requests.get(url, headers=headers)
    assert response.status_code == 200, f'Unsuccessful attempt for request "{url}"'
    return response

def _responseText_(response, filepath):
    '''save response.text to a file'''
    with open(filepath, 'w', encoding='utf-8') as fp:
        fp.write(response.text)

def downloadURL(url, filepath, recursive=True, maxAttempt=5, sleepTime=10):
    '''
    download response.text from a given url
    '''
    if recursive: r = _recursiveRequest_(url, maxAttempt=maxAttempt, sleepTime=sleepTime)
    else: r = _requestURL_(url)
    _responseText_(r, filepath)

def _recursiveRequest_(url, maxAttempt=5, sleepTime=10):
    '''
    recursively request a given url, with maximum attempt given in `maxAttempt`
    '''
    successRequest = False; attemptCount = 0
    while not successRequest:
        try: r = _requestURL_(url); return r
        except:
            attemptCount += 1
            if attemptCount > maxAttempt: raise ValueError('Request failed...')
            time.sleep(sleepTime)
