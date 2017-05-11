#!/usr/bin/python

import os,sys
from setuptools import setup

PACKAGES = ['pyampli']
NAME = 'pyAmpli'
DESCRIPTION = 'Python package for amplicon-based variant filtering (germline and somatic)'
AUTHOR = 'Matthias Beyens'
AUTHOR_EMAIL = 'matthias.beyens@uantwerpen.be'
URL = 'https://github.com/MBeyens/pyAmpli'
VERSION = '0.1.0'
LICENSE = 'GPL'
DATAFILES=[('pyampli',['pyampli/config.yaml'])]

if sys.version_info.major == 2 and sys.version_info.minor != 7:
    sys.stderr.write("ERROR: pyAmpli is only for python 2.7 or greater you are running %d.%d\n", (sys.version_info.major, sys.version_info.minor))
    sys.exit(1)


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
    install_requires = ['pysam==0.8.4','pyvcf==0.6.8','pyyaml==3.11', 'setuptools == 20.2.2'
    ]

)