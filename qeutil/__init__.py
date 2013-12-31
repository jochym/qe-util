#!/usr/bin/python
# -*- coding: utf-8 -*-

"""This module defines an interface to Quantum-Espresso.
http://www.quantum-espresso.org/
"""



from ase import Atom, Atoms
import os
from glob import glob
from os.path import join, isfile, islink

import numpy as np

from ase.data import atomic_numbers
from ase.units import Bohr, Hartree, Rydberg, fs
from ase.data import chemical_symbols
from ase.calculators.calculator import FileIOCalculator, Parameters, kpts2mp, ReadError

from quantumespresso.io import *


class QuantumEspresso(FileIOCalculator):
    """Class for doing Quantum Espresso calculations.

    The default parameters are very close to those that 
    the QE Fortran code would use.  These are the exceptions::

    """

    implemented_properties = ['energy', 'forces', 'stress']

    pass


