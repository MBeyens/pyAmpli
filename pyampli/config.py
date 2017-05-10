#!/usr/bin/python

import os, sys, logging, yaml


def config_load(config_filename):
    
    try:
        config_filename = os.path.realpath(config_filename)
        with open(config_filename, 'r') as yml_file:
            cfg = yaml.load(yml_file)
    logging.info('Using config file: %s', config_filename)
        return cfg
    except IOError:
        logging.error('Cannot find/load config file (%s)', config_filename, exc_info=True)
        sys.exit(0)
