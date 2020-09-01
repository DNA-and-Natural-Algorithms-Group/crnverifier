#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'crnverifier',
    version = '0.1',
    description = 'Verify the equivalence of chemical reaction networks (CRNs), or the correctness of an implementation CRN with respect to a formal CRN.',
    long_description = LONG_DESCRIPTION,
    author = 'Stefan Badelt, Seung Woo Shin, Robert Johnson, Qing Dong, Erik Winfree',
    author_email = 'winfree@caltech.edu',
    #url = 'http://www.github.com/DNA-and-Natural-Algorithms-Group/crnverifier/',
    license = 'MIT',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        ],
    #install_requires = ['pyparsing'],
    test_suite = 'tests',
    packages = ['crnverifier'],
    entry_points = {
        'console_scripts': [
            'crnverifier=crnverifier.verifier:main'
            ],
        }
)
