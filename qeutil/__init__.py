#!/usr/bin/python
# -*- coding: utf-8 -*-

"""This module defines an interface to Quantum-Espresso.
http://www.quantum-espresso.org/
"""

#    Copyright 1998-2014 by Paweł T. Jochym <pawel.jochym@ifj.edu.pl>
#
#    qe-util is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    qe-util is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Elastic.  If not, see <http://www.gnu.org/licenses/>.


import os
import sys
from glob import glob
from os.path import join, isfile, islink
import subprocess
import time

import numpy as np

import ase.units
from ase.units import Bohr, Hartree, Rydberg, fs
from ase import Atom, Atoms
from ase.data import atomic_numbers, chemical_symbols
from ase.calculators.calculator import Calculator, FileIOCalculator, Parameters, kpts2mp, ReadError, all_changes

# Doubtful practice should be more selective in future
from readers import *
from writers import *
from qeio import *
from analyzers import *

import copy

# The exception for the calc runnin but not ready.
class CalcNotReadyError(Exception):
    pass


class QuantumEspresso(FileIOCalculator):
    """Class for doing Quantum Espresso calculations.

    The default parameters are very close to those that 
    the QE Fortran code would use.  These are the exceptions::

    """

    implemented_properties = ['energy', 'forces', 'stress']
    
    pw_cmd='mpiexec -n %(procs)d  pw.x <pw.in >pw.out'
    ph_cmd='mpiexec -n %(procs)d  ph.x <ph.in >ph.out'
    matdyn_cmd='mpiexec -n %(procs)d matdyn.x <matdyn.in >matdyn.out'
    phdos_cmd= 'mpiexec -n %(procs)d matdyn.x <phdos.in  >phdos.out'
    q2r_cmd='q2r.x <q2r.in >q2r.out'

    cmds={
        'scf': pw_cmd,
        'ph': ph_cmd,
        'matdyn': matdyn_cmd,
        'dos': phdos_cmd,
        'q2r': q2r_cmd
    }

    default_parameters = {
                'calc':'scf',
                'outdir': 'tmp',
                'wfcdir': 'tmp',
                'asr':'crystal',
                'pseudo_dir': '../pspot',
                'ndos': 200,
                'points': 100,
                'tstress': True,
                'tprnfor': True,
                'ecutwfc': 50,
                'ibrav': 0,
                'pp_type': 'nc',
                'pp_format': 'ncpp',
                'use_symmetry': True,
                'kpt_type': 'automatic',
                'kpt_shift': [0,0,0],
                'procs': 1
    }
    'Default parameters'

    def __init__(self, label=None, atoms=None, wdir='./', restart=None, **kwargs):
        FileIOCalculator.__init__(self,label=label,atoms=atoms,command=None,**kwargs)
        self.label=label
        self.prefix=label
        self.wdir=wdir
        if restart :
            self.directory=restart
        else :
            self.directory=make_calc_dir(self.prefix,wdir)
        self.submited=False

    def copy(self):
        c=copy.deepcopy(self)
        c.directory=make_calc_dir(c.prefix,c.wdir)
        return c

    def build_command(self, prop=['energy'], params={}):
        params=params.copy()
        params.update({'prop': prop})
        cmd=''
        if {'energy','stress','forces'} & set(prop):
            cmd+= self.pw_cmd % params
            cmd+='\n'
        if 'd2' in prop :
            cmd+= self.ph_cmd % params
            cmd+='\n'
            cmd+= self.q2r_cmd % params
            cmd+='\n'
        if 'frequencies' in prop :
            cmd+= self.matdyn_cmd % params
            cmd+='\n'
        if 'phdos' in prop :
            cmd+= self.phdos_cmd % params
            cmd+='\n'
        return cmd
            

    def run_calculation(self, atoms, properties, system_changes):
        #print self.command
        FileIOCalculator.calculate(self, atoms, properties, system_changes)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

            self.command=self.build_command(properties,self.parameters)
            self.run_calculation(atoms, properties, system_changes)
        
#        if {'energy','stress'} & set(properties) :
#            self.command=self.build_command(properties,self.parameters)
#            prop=list({'energy','stress'} & set(properties))
#            self.run_calculation(atoms, prop, system_changes)
#        elif 'bands' in properties :
#            raise NotImplementedError
#        elif 'd2' in properties :
#            self.command=self.build_command(['d2'],self.parameters)
#            self.run_calculation(atoms, ['d2'], system_changes)
#        elif 'frequencies' in properties :
#            self.command=self.build_command(['frequencies'],self.parameters)
#            self.run_calculation(atoms, ['frequencies'], system_changes)
#        elif 'phdos' in properties :
#            self.command=self.build_command(['phdos'],self.parameters)
#            self.run_calculation(atoms, ['phdos'], system_changes)
#        else :
#            raise NotImplementedError


    def write_input(self, atoms=None, properties=None, system_changes=None):
        """
        Write input file(s). This is called by FileIOCalculator just before 
        it executes the external command.
        """
        self.set(prefix=self.prefix)
        self.atoms2params()
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        
        if {'energy','stress','forces'} & set(properties) :
            write_pw_in(self.directory, atoms, self.parameters)
        if 'd2' in properties :
            write_ph_in(self.directory, atoms, self.parameters)
            write_q2r_in(self.directory, atoms, self.parameters)
        if 'frequencies' in properties :
            write_matdyn_in(self.directory, atoms, self.parameters)
        if 'phdos' in properties:
            write_phdos_in(self.directory, atoms, self.parameters)



    def read_results(self):
        """Read energy, forces, ... from output file(s)."""
        
        fn=os.path.join(self.directory,'pw.out')
        xml_fn=os.path.join(self.directory,self.parameters['outdir'],self.prefix+'.save','data-file.xml')
        # Read the pan-ultimate line of the output file
        try: 
            ln=open(fn).readlines()[-2]
            if ln.find('JOB DONE.')>-1 :
                # Job is done we can read the output
                r=read_quantumespresso_textoutput(fn)
                self.results.update(r)
                self.results_xml={}
                self.results_xml.update(read_quantumespresso_xmloutput(xml_fn,'all'))
                self.results['energy']=r['etotal']
                s=-array(r['stress'])* 1e-1 * ase.units.GPa
                self.results['stress']=array([s[0, 0], s[1, 1], s[2, 2],
                                       s[1, 2], s[0, 2], s[0, 1]])
                rk=self.results.keys()
                if 'cell' in rk :
                    try :
                        self.results['cell']=array(self.results['cell'])*self.results['alat']*Bohr
                    except TypeError:
                        # cell data not in order just clear the mess
                        del(self.results['cell'])
                if 'atoms_forces' in rk :
                    # Translate from the Ry/au of QE to the eV/A units of ASE
                    self.results['forces']=array(self.results['atoms_forces'])/(ase.units.Rydberg/ase.units.Bohr)

                try :
                    self.results['phdos']=np.loadtxt(os.path.join(self.directory,self.prefix+'.dos')).T
                except IOError:
                    # No frequency data. We assume there is none and silently
                    pass

                try :
                    self.results['frequencies']=np.loadtxt(os.path.join(self.directory,self.prefix+'.freq.gp')).T
                except IOError:
                    # No frequency data. We assume there is none and silently
                    pass
            else :
                # Job not ready.
                raise CalcNotReadyError
        except (IOError, IndexError) :
            # Job not ready.
            raise CalcNotReadyError
        
    def clean(self):
        '''
        Celaning procedure. Non-functional for now.
        We need to decide which files should be left intact.
        '''
        pass
        
    def atoms2params(self):
        '''
            Populate parameters dictionary with data extracted from atoms.
            Return the updated dictionary.
        '''
        p=self.parameters
        return p



class RemoteQE(QuantumEspresso):
    
    # Queue system submit command
    qsub_tool='qsub'
    qstat_tool='qstat'
    
    qsub_cmd='cd %(rdir)s ; %(qsub_tool)s -N %(title)s -l procs=%(procs)d ./run-pw.pbs'

    # Remote execution command
    remote_exec_cmd='ssh %(user)s@%(host)s "%(command)s"'

    # If you cannot mount the data directory into your system it is best 
    # to use the rsync command to transfer the results back into the system.

    # Command for copying the data out to the computing system
    copy_out_cmd='rsync -a "%(ldir)s" "%(user)s@%(host)s:%(rdir)s"'
    # Command for copying the data in after the calculation
    copy_in_cmd='rsync -a "%(user)s@%(host)s:%(rdir)s" "%(ldir)s"'

    # Template for the PBS batch job
    pbs_template=''

    # Command to check the state of the job
    pbs_check_cmd='''%(qstat_tool)s -f %(jobid)s |grep job_state |awk '{print $3}' '''
    
    # Access data
    host=''
    user=''
    
    # Location:
    # local working directory
    wdir='.'
    # Remote working directory relative to the home directory or absolute
    rdir='.'

    # Repetition timer (seconds) for checkin the state of the job.
    job_check_time=15
    
    def __init__(self, block=True, **kwargs):
        QuantumEspresso.__init__(self, user=self.user, host=self.host, 
                                wdir=self.wdir, rdir=self.rdir, **kwargs)
        self.jobid=None
        self.block=block

    def write_pbs_in(self,properties):
        fh=open(os.path.join(self.directory,'run-pw.pbs'),'w')

        fh.write(self.pbs_template % {
            'command': QuantumEspresso.build_command(self,prop=properties,
                                                        params=self.parameters)
            })        
        
        fh.close()

    def build_command(self,prop=['energy'],params={}):
        cmd=self.qsub_cmd % {
            'qsub_tool': self.qsub_tool,
            'qstat_tool': self.qstat_tool,
            'title': self.label,
            'procs': self.parameters['procs'],
            'rdir': os.path.join(self.parameters['rdir'],os.path.split(self.directory)[-1])
        }
        cmd=self.remote_exec_cmd % {
                'command': cmd,
                'user': self.parameters['user'],
                'host': self.parameters['host']
         }
        return cmd

    def write_input(self, atoms=None, properties=['energy'], system_changes=all_changes):
        '''Write input file(s).'''
        QuantumEspresso.write_input(self, atoms, properties, system_changes)
        self.write_pbs_in(properties)
        subprocess.call(self.copy_out_cmd % {
                            'ldir': self.directory,
                            'rdir': self.parameters['rdir'],
                            'user': self.parameters['user'],
                            'host': self.parameters['host']
                        }, shell=True)

    def job_ready(self):
        try :
            cmd=self.remote_exec_cmd % {
                'command': self.pbs_check_cmd % {
                                    'qsub_tool': self.qsub_tool,
                                    'qstat_tool': self.qstat_tool,
                                    'jobid':self.jobid
                                    },
                'user': self.parameters['user'],
                'host': self.parameters['host']
                }
            state=subprocess.check_output(cmd, shell=True).split()[-1]
        except (subprocess.CalledProcessError, IndexError) :
            # Unknown state. We assume it has finished and continue
            state='N'

        return not (state in ['Q','R'])


    def run_calculation(self, atoms=None, properties=['energy'],
                            system_changes=all_changes):
        '''
        Internal calculation executor. We cannot use FileIOCalculator
        directly since we need to support remote execution.
        
        This calculator is different from others. 
        It prepares the directory, launches the remote process and
        raises the exception to signal that we need to come back for results
        when the job is finished.
        '''
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please configure RemoteQE calculator!')
        olddir = os.getcwd()
        errorcode=0
        try:
            os.chdir(self.directory)
            output = subprocess.check_output(self.command, shell=True)
            self.jobid=output.split()[0]
            self.submited=True
            #print "Job %s submitted. Waiting for it." % (self.jobid)
            # Waiting loop. To be removed.
        except subprocess.CalledProcessError as e:
            errorcode=e.returncode
        finally:
            os.chdir(olddir)
        
        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))
        self.read_results()


    def read_results(self):
        """Read energy, forces, ... from output file(s)."""
        
        if self.submited:
            # The job has been submitted. Check the state.
            if not self.job_ready() :
                if self.block :
                    while not self.job_ready() :
                        time.sleep(self.job_check_time)
                else :
                    raise CalcNotReadyError
                    

            # Assume the calc finished. Copy the files back.
            subprocess.call(self.copy_in_cmd % {
                'ldir': self.wdir,
                'rdir': os.path.join(self.parameters['rdir'],os.path.split(self.directory)[-1]),
                'user': self.parameters['user'],
                'host': self.parameters['host']
            }, shell=True)
                
        
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
                self.jobid=None
            else :
                # Job not ready.
                raise CalcNotReadyError
        except (IOError, IndexError) :
            # Job not ready.
            raise CalcNotReadyError
        
        # All is fine - read the results
        QuantumEspresso.read_results(self)

    @classmethod
    def ParallelCalculate(cls,syslst,properties=['energy'],system_changes=all_changes):
        '''
        Run a series of calculations in parallel using (implicitely) some 
        remote machine/cluster. The function returns the list of systems ready
        for the extraction of calculated properties.
        '''
        print 'Launching:',
        sys.stdout.flush()
        for n,s in enumerate(syslst):
            try :
                s.calc.block=False
                s.calc.calculate(atoms=s,properties=properties,system_changes=system_changes)
            except CalcNotReadyError:
                s.calc.block=True
            print n+1, 
            sys.stdout.flush()
        print
        print '     Done:',
        sys.stdout.flush()
        for n,s in enumerate(syslst):
            s.calc.read_results()
            print n+1, 
            sys.stdout.flush()
        print
        return syslst

