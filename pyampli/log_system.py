#!/usr/bin/python

import logging, os, sys


def initiate_log(debug_modes):
    os_version_check()
    if not debug_modes:
        logging.basicConfig(filename=os.getcwd() + '/log.txt',
                            filemode='w',
                            format='%(asctime)s \t-\t %(name)s \t-\t %(levelname)s \t-\t %(message)s',
                            datefmt='%d/%m/%Y %I:%M:%S',
                            level=logging.INFO)
    else:
        logging.basicConfig(filename=os.getcwd() + '/log.txt',
                            filemode='w',
                            format='%(asctime)s \t-\t %(name)s \t-\t %(levelname)s \t-\t %(message)s',
                            datefmt='%d/%m/%Y %I:%M:%S',
                            level=logging.NOTSET)
    stderrLogger = logging.StreamHandler()
    stderrLogger.setFormatter(logging.Formatter('%(asctime)s \t-\t %(name)s \t-\t %(levelname)s \t-\t %(message)s'))
    logging.getLogger().addHandler(stderrLogger)

    return logging


def os_version_check():
    os_check = 1
    os_platform = sys.platform.lower()
    # linux
    if os_platform == 'linux' or os_platform == 'linux2':
        os_check = 0
    # OS X
    elif os_platform == 'darwin':
        os_check = 0

    if not os_check:
        logging.error('Your operating system is not supported! Please contact matthias.beyens@uantwerpen.be')
        sys.exit(0)
    else:
        logging.info('Your operating system is supported (%s)', os_platform)