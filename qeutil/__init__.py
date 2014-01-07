#!/usr/bin/python
# -*- coding: utf-8 -*-

"""This module defines an interface to Quantum-Espresso.
http://www.quantum-espresso.org/
"""

import os
from glob import glob
from os.path import join, isfile, islink

import numpy as np

import ase.units
from ase.units import Bohr, Hartree, Rydberg, fs
from ase import Atom, Atoms
from ase.data import atomic_numbers, chemical_symbols
from ase.calculators.calculator import FileIOCalculator, Parameters, kpts2mp, ReadError, all_changes

# Doubtful practice should be more selective in future
from readers import *
from writers import *
from qeio import *
from analyzers import *

# The exception for the calc runnin but not ready.
class CalcNotReadyError(Exception):
    pass


class QuantumEspresso(FileIOCalculator):
    """Class for doing Quantum Espresso calculations.

    The default parameters are very close to those that 
    the QE Fortran code would use.  These are the exceptions::

    """

    implemented_properties = ['energy', 'forces', 'stress']
    
    pw_cmd='pw.x <pw.in >pw.out'
    ph_cmd='ph.x <ph.in >ph.out'
    matdyn_cmd='matdyn.x <matdyn.in >matdyn.out'
    phdos_cmd= 'matdyn.x <phdos.in  >phdos.out'
    q2r_cmd='q2r.x <q2r.in >q2r.out'

    # Command for copying the data out to the computing system
    copy_out_cmd=''
    
    # Command for copying the data in after the calculation
    # If you cannot mount the data directory into your system it is best 
    # to use the rsync command to transfer the results back into the system.
    copy_in_cmd=''
    
    # Queue system submit command
    qsub_cmd=''

    default_parameters = {
                'calc':'scf',
                'outdir': 'tmp',
                'wfcdir': 'tmp',
                'asr':'crystal',
                'pseudo_dir': '../pspot',
                'ndos': 200,
                'points': 100,
                'tstress': True,
                'ecutwfc': 50,
                'ibrav': 0,
                'pp_type': 'nc',
                'pp_format': 'ncpp',
                'use_symmetry': False,
                'kpt_type': 'automatic',
                'kpt_shift': [0,0,0]
    }
    'Default parameters'

    def __init__(self,label=None,atoms=None, wdir='./', **kwargs):
        FileIOCalculator.__init__(self,label=label,atoms=atoms,command=None,**kwargs)
        self.label=label
        self.prefix=label
        self.directory=make_calc_dir(self.prefix,wdir)
        self.submited=False
        self.use_symmetry=self.parameters['use_symmetry']

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if self.submited :
            raise CalcNotReadyError
        try :
            if {'energy','stress'} and set(properties) :
                self.command=self.pw_cmd % {'infile':'pw.in', 'outfile':'pw.out'}
                FileIOCalculator.calculate(self, atoms, 
                    list({'energy','stress'} and set(properties)),
                    system_changes)
            elif 'phonons' in properties :
                self.calculate(self, atoms, ['energy'], system_changes)
                FileIOCalculator.calculate(self, atoms, ['phonons'], system_changes)
        except CalcNotReadyError :
            # The command finished without error but the 
            # results are not ready. 
            # We assume this is a remote job submitted to 
            # some queuing system.
            self.submited=True

            

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input file(s).
            
        """
        self.set(prefix=self.prefix)
        self.atoms2params()
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write_pw_in(self.directory, self.atoms, self.parameters)

    def read_results(self):
        """Read energy, forces, ... from output file(s)."""
        
        fn=os.path.join(self.directory,'pw.out')
        # Read the pan-ultimate line of the output file
        try: 
            ln=open(fn).readlines()[-2]
            if ln.find('JOB DONE.')>-1 :
                # Job is done we can read the output
                r=read_quantumespresso_textoutput(fn)
                self.results['energy']=r['etotal']
                s=array(r['stress'])* 1e-1 * ase.units.GPa
                self.results['stress']=array([s[0, 0], s[1, 1], s[2, 2],
                                       s[1, 2], s[0, 2], s[0, 1]])
                self.submited=False
            else :
                # Job not ready.
                raise CalcNotReadyError
        except IOError :
            # Job not ready.
            raise CalcNotReadyError
        
    def atoms2params(self):
        '''
            Populate parameters dictionary with data extracted from atoms.
            Return the updated dictionary.
        '''
        p=self.parameters
        return p


