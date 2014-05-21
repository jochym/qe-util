# -*- coding: utf-8 -*-

from distutils.core import setup

setup(
    name='qeutil',
    version='0.1.16',
    packages=['qeutil'],
    license='GPLv3',
    description = 'A set of utilities for Quantum Espresso',
    author = 'Paweł T. Jochym',
    author_email = 'Pawel.Jochym@ifj.edu.pl',
    url = 'https://github.com/jochym/qe-util',   
    download_url = 'https://github.com/jochym/qe-util/', 
    keywords = ['science', 'physics', 'quantum-espresso'], 
    requires = ['ase','numpy','scipy','pyspglib'],
    classifiers = [],
)
