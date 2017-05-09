#!/usr/bin/python

import os
from setuptools import setup

PACKAGES = ['pyampli']
NAME = 'pyAmpli'
DESCRIPTION = 'Python package for amplicon filtering (germline and somatic)'
AUTHOR = 'Matthias Beyens'
AUTHOR_EMAIL = 'matthias.beyens@uantwerpen.be'
URL = 'https://github.com/MBeyens/pyAmpli'
VERSION = '0.1.0'
LICENSE = 'GPL'
DATAFILES=[('pyampli',['pyampli/config.yaml'])]

def read(*paths):
    """ read files """
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

setup(
    name=NAME,
    version=VERSION,
    packages=PACKAGES,
    url=URL,
    license=LICENSE,
    author=AUTHOR,
    data_files=DATAFILES,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=(read("README.md")),
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7'
    ],
    scripts = ['pyampli/pyAmpli.py'],
    install_requires = ['pysam','pyvcf','pyyaml'
    ]

)
