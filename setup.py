#!/usr/bin/python

import os
from distutils.core import setup

PACKAGES = ['bin']
NAME = 'pyAmpli'
DESCRIPTION = 'Python package for amplicon filtering (germline and somatic)'
AUTHOR = 'Matthias Beyens'
AUTHOR_EMAIL = 'matthias.beyens@uantwerpen.be'
URL = 'https://github.com/MBeyens/pyAmpli'
VERSION = '0.8.0'
LICENSE = 'GPL'


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
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=(read("README.md")),
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7'
    ],
    scripts = ['bin/pyAmpli'],
    install_requires = ['pysam','pyvcf'
    ]

)
