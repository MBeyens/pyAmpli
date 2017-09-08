#!/usr/bin/python

import os, sys, logging, yaml


def config_load(config_filename):

    if os.path.isfile(infile):
        config_filename = os.path.realpath(config_filename)
        try:
            with open(config_filename, 'r') as yml_file:
                cfg = yaml.load(yml_file)
                logging.info('Using configuration file: %s', config_filename)
            return cfg
        except:
            logging.error('Cannot load configuration file (%s), although it exists', config_filename, exc_info=True)
            sys.exit(0)
    else:
        logging.error('Cannot find configuration file (%s)', config_filename, exc_info=True)
        sys.exit(0)