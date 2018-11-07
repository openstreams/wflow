"""
wflow_logging -- Logging method

Deltares 2013

"""

import logging
import logging.handlers


def setuplog(logfilename, loggername):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    
    Input:
        
        - logfilename - filename to log to (console is also used)
        - loggername - name of this logger
    """
    try:
        # create logger
        logger = logging.getLogger(loggername)
        logger.setLevel(logging.DEBUG)
        ch = logging.FileHandler(logfilename)
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        # create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
        )
        # add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        # add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
        return None
