#!/usr/bin/python

import os, sys, logging, yaml


def config_load():
    try:
        config_filename = os.path.dirname(os.path.realpath(__file__)) + '/config.yaml'
        with open(config_filename, 'r') as yml_file:
            cfg = yaml.load(yml_file)
        return cfg
    except IOError:
        logging.error('Cannot find/load config file (%s)', config_filename, exc_info=True)
        sys.exit(0)
