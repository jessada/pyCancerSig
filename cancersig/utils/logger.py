import logging
import sys

lg = logging.getLogger('cancersig')
lg.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('## [%(levelname)s]%(asctime)25s - %(name)-12s - %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
lg.addHandler(ch)

def set_log_file(file_name, level=logging.DEBUG):
    # create file handler which logs even debug messages
    fh = logging.FileHandler(file_name)
    fh.setLevel(level)
    # create formatter and add it to the handlers
    fh.setFormatter(formatter)
    # add the handlers to the logger
    lg.addHandler(fh)


def getLogger(name):
    global lg
    lg = logging.getLogger(name)

def info(msg=""):
    lg.info(msg)

def warning(msg=""):
    lg.warning(msg)

def dbg(msg=""):
    lg.debug(msg)

def throw(err_msg):
    lg.exception(err_msg)
