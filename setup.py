#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'crnverifier',
    version = '0.3',
    description = 'Verify the equivalence of chemical reaction networks (CRNs), or the correctness of an implementation CRN with respect to a formal CRN.',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = "text/markdown",
    author = 'Stefan Badelt, Seung Woo Shin, Robert Johnson, Qing Dong, Erik Winfree',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'http://www.github.com/DNA-and-Natural-Algorithms-Group/crnverifier/',
    license = 'MIT',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.8',
    install_requires = ['pyparsing'],
    packages = ['crnverifier'],
    test_suite = 'tests',
    entry_points = {
        'console_scripts': [
            'crnverifier=crnverifier.verifier:main'
            ],
        }
)
