"""Logging file
"""
import os, logging
from datetime import datetime

## Version
__version__ = '1.0.0'

## create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
#logging.basicConfig(datefmt='%Y-%m-%d %H:%M:%S')

## set file handler and formatter
abs_path = os.path.dirname(os.path.abspath(__file__))
abs_root_path = os.path.dirname(abs_path)
log_dir_path = abs_root_path+os.sep+'log'
if not os.path.exists(log_dir_path):
    os.makedirs(log_dir_path)

now = datetime.now()
log_filename = now.strftime("%Y-%m-%d-%H%M%S")

fh = logging.FileHandler(log_dir_path + os.sep + log_filename + '.log')
fh.setLevel(logging.WARNING)
f_format = logging.Formatter('%(asctime)s %(module)s.%(funcName)s: %(levelname)s: %(message)s')
fh.setFormatter(f_format)
logger.addHandler(fh)

## set console hander and formatter
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
c_format = logging.Formatter('%(module)s.%(funcName)s: %(levelname)s: %(message)s')
ch.setFormatter(c_format)
logger.addHandler(ch)

logger.debug(f"Logger set successfully in {__name__}")
