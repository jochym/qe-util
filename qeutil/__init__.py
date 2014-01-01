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

# Doubtful practice should be more selective in future
from readers import *
from writers import *
from qeio import *
from analyzers import *


class QuantumEspresso(FileIOCalculator):
    """Class for doing Quantum Espresso calculations.

    The default parameters are very close to those that 
    the QE Fortran code would use.  These are the exceptions::

    """

    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self,label=None,atoms=None,
                    pw_cmd='pw.x',
                    ph_cmd='ph.x',
                    matdyn_cmd='matdyn.x',
                    q2r_cmd='q2r.x',
                    **kwargs):
        FileIOCalculator.__init__(self,label=label,atoms=atoms,command=None,**kwargs)
        self.atoms=atoms.copy()
        self.label=label
        self.pw_cmd=pw_cmd
        self.ph_cmd=ph_cmd
        self.matdyn_cmd=matdyn_cmd
        self.q2r_cmd=q2r_cmd




