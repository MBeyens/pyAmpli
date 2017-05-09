#!/usr/bin/python

import logging, os


def initiate_log(debug_modes):
    if not debug_modes:
        logging.basicConfig(filename=os.path.dirname(os.path.realpath(__file__)) + '/log.txt',
                            filemode='w',
                            format='%(asctime)s \t-\t %(name)s \t-\t %(levelname)s \t-\t %(message)s',
                            datefmt='%d/%m/%Y %I:%M:%S',
                            level=logging.INFO)
    else:
        logging.basicConfig(filename=os.path.dirname(os.path.realpath(__file__)) + '/log.txt',
                            filemode='w',
                            format='%(asctime)s \t-\t %(name)s \t-\t %(levelname)s \t-\t %(message)s',
                            datefmt='%d/%m/%Y %I:%M:%S',
                            level=logging.NOTSET)
    stderrLogger = logging.StreamHandler()
    stderrLogger.setFormatter(logging.Formatter('%(asctime)s \t-\t %(name)s \t-\t %(levelname)s \t-\t %(message)s'))
    logging.getLogger().addHandler(stderrLogger)

    return logging
