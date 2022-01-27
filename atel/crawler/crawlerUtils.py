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
