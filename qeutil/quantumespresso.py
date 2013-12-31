"""
This module defines an ASE interface to QUANTUM-ESPRESSO

http://www.quantum-espresso.org
"""

# the string
#           ## [missed]
# tags points where I do have uncertainties

import os
import subprocess
from glob import glob
from os.path import join, isfile, islink

import numpy as np
import string

from ase.data import chemical_symbols
from ase.data import atomic_numbers, atomic_masses
from ase.units import Bohr, Hartree, GPa


class QEpw:
    """Class for doing QUANTUM-ESPRESSO calculations.

    The default parameters are very close to those that the QE
    Fortran code would use.  

      calc = QE(label='qe-pw', xc='scf', pulay=5, mix=0.1)

    Use the set_inp_X method to set extra INPUT parameters, where
    X stands for (control | system | electrons | ions). For instance

      calc.set_inp_control('nstep', 30)

    """

    __name__ = "QEpw"
    __version__ = "1.0"
    calculations_available = ['scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md']

    def __init__(self, atoms=None, label='pw', xc='scf', nstep=50, etot_conv_thr=1.0e-4,       #CONTROL
                 restart_mode=None, tstress=None, prefix='pwscf',
                 dt=None, pseudodir=None,
                 ibrav=0, cell_parameters_units=None, celldm=None, crystal = None,nbands=None, #SYSTEM
                 ecutwfc=None, ecutrho=None, occupations=None,
                 input_dft=None, smearing='gauss', degauss=None, nspin=None,
                 electron_maxstep=None, conv_thr=None, adaptive_thr=None, conv_thr_init=None,  #ELECTRONS
                 conv_thr_multi=None, mixing_mode=None, mixing_beta=None, diag=None,
                 spoint=None, swfc=None,
                 cell_dynamics=None, press=None, wmass=None, press_conv_thr=None,              #CELL
                 cell_dofree=None,
                 kpts=None,                                                                    #K_POINTS
                 kpts_type=None,                                                               
                 cell_parameters=None):                                                        #CELL_PARAMETERS
    
        """Construct QuantumEspresso-calculator.

        Parameters
        ==========

        ---------- CONTROL section
        
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'qe'.
        xc: str
            the operation QE is supposed to perform (scf, nscf, bands, relax, md, vc-relax, vc-md)
            Default is 'scf'.
        nstep: integer
            number of ionic + electronic steps
            Default is 1 if calc is 'scf', 'nscf', 'bands'; 50 for the other cases
        etot_conv_thr: float
            convergence threshold on total energy (a.u) for ionic minimization
            Default is 1.0e-4

        prefix: str
            to be written in the input file
            Default is 'pwscf'
        outdir: str
            input, temporary, output files are found in this directory
            Default is value of the ESPRESSO_TMPDIR environment variable if set; current directory ('./') otherwise
        pseudodir: str
            directory containing pseudopotential files
            Default is value of the $ESPRESSO_PSEUDO environment variable if set; '$HOME/espresso/pseudo/' otherwise 

        ---------- SYSTEM section

        ibrav: integer
            Bravais-lattice index (see QE manual for insights).Might be cell parameters, idx of BRAVET grid+ cell dimension, crystal param.
            Default expectation are cell parameters, since you are supposed to inherit cell from atoms object
            Required
        units: string
            either 'bohr' or 'angstrom', it is needed since ase does not handle units of measure when building an atomic system.
            So, if ibrav = 0 and celldm(1) is given, it is given in Bohr and the CELL card is well-defined. Otherwise, we need to
            specify the units (either Bohr or Angstrom) in which the Cell is given.
            This variable is the way we do it.
            Default value here is Bohr, in case alat is not given.
        celldm(i), i=1,6: float
            Crystallographic constants
            must be passed as a list even for a single number. only those needed might be specified.
        crystal (A, B, C, cosAB, cosAC, cosBC) : float
            Traditional crystallographic constants
            must be passed as a list with 6 entries
        nat: integer
            number of atoms in the unit cell
            Default obtained from Atoms
        ntyp: integer
            number of types of atoms in the unit cell
            Default obtained from Atoms
        nbands: integer
            number of electronic states (bands) to be calculated
            Default (see QE documentation)
        ecutwfc: float
            kinetic energy cutoff (Ry) for wavefunctions
            Required
        ecutrho: float
            kinetic energy cutoff (Ry) for charge density and potential
            Default is 4 * ecutwfc
        smearing: str          
            'gaussian'/'gauss', 'methfessel-paxton'/'m-p'/'mp', 'marzari-vanderbilt'/'cold'/'m-v'/'mv', 'fermi-dirac'/'f-d'/'fd'
            Default is 'gaussian'
        input_dft: str
            Exchange-correlation functional: eg 'PBE', 'BLYP', etc. Overrides the value read from pseudopotential files.
            Default read from pseudopotential files

        ---------- ELECTRONS section

        electron_maxstep: int
            maximum number of iterations in a scf step
            default is 100
        conv_thr: float
            conv. threshold for selconsistency
            default is 1.0e-6
        adaptive_thr: logical

        conv_thr_init: float
            if adaptive_thr is true, this is used for first cycle
            default is 1.0e-3
        conv_thr_multi: float
            if adaptive_thr is true, the convergence th for each cycle is given by min(conv_thr, conv_thr_multi * dexx)
            default is 0.1
        mixing_mode: str
            may be 'plain', 'TF', 'local-TF'
            Default is 'plain'
        mixing_beta: float
            mixing factor for self-consistency
            Default is 0.7
        diag: str
            may be 'david', 'cg'
            Default is 'david'
        spoint: str
            'atomic', 'file'
            Default for scf, *relax, *md is 'atomic'; for others the only possibility is 'file'
        swfc: str
            'atomic', 'atomic+random', 'random', 'file'
            Default is 'atomic+random'

        ---------- CELL section
        present only if calculation is vc-relax or vc-md

        ---------- CELL PARAMETERS card

        cell_parameters: float
            Optional card, needed only if ibrav = 0 is specified
            made by 3 vectors

        cell_parameters_units: str
            alat | bohr | angstrom
            
            
        """
        
        # initialize dictionaries for cards in input files
        self.inp_control = {}
        self.inp_control_writeall = {}
        self.inp_system = {}
        self.inp_system_writeall = {}
        self.inp_electrons = {}
        self.inp_electrons_writeall = {}        
        self.inp_ions = {}
        self.inp_ions_writeall = {}      
        self.inp_cell = {}
        self.inp_cell_writeall = {}

        # initialize key variables that must exist
        self.error = False
        self.error_message=[]
        self.if_pos = None
        self.PP = None
        self.kpts =[]
        self.kpts_num= 0
        self.nelect = 0
        self.alat = 0.0
        self.etotal = 0.0
        self.etotal_accuracy = 0.0
        self.dt = 0.0
        self.conv_thr=conv_thr
        self.degauss=degauss
        
        
        # | CONTROL input parameters
        
        self.label = label
        if not xc.lower() in self.calculations_available:
            raise ValueError, "The specified calculation (%s) is not available for this calculator" % xc.lower()

        self.calculation = xc.lower()
        if( self.calculation in ['vc-relax', 'vc-md']) and cell_dynamics is None:
            raise ValueError, 'your calculation requires you specify the cell-dynamic'

        self.nstep = nstep
        self.etot_conv_thr = etot_conv_thr

        if not dt is None:
            self.set_inp_control('dt', dt)
            self.dt = dt

        if not restart_mode is None:
            self.set_inp_control('restart_mode', restart_mode)

        if not tstress is None:
            self.set_inp_control('tstress', tstress)
            self.tstress = tstress

        if not prefix is None:                              # in case it is not given, default will be use by the code.
            self.prefix = prefix
            self.set_inp_control('prefix', prefix)
            
        # | SYSTEM input parameters

        self.set_ibrav_values()
        if (ibrav is None):
            # note that default value is 0, which means to copy the cell matrix from atoms structure
            # so being None here might point to some problem
            raise ValueError('Bravais-lattice index is required and not set')
        if not ibrav in self.ibrav_allowed:
            raise ValueError('ibrav given value is not recognizable as valid')
            
        self.ibrav = ibrav
        self.set_inp_system('ibrav', ibrav)

        if cell_parameters_units is None:
            if ((celldm is None) and (crystal is None)):
                # in this case cell is not specified at all; actually, that might be the most common case,
                # with ibrav=0 a<nd the cell matrix being adopted from Atoms object
                self.celldm = None
                self.cell_parameters_units = "bohr"
                # NOTE: self.alat is 0.0 at this point: then it will become cell[0]
            elif ((not celldm is None) or (not crystal is None)): 
                # check whether the given celldm is consistent with the specified ibrav value
                # and sets it with the given values
                self.set_celldm(celldm, crystal)
                self.cell_parameters_units = "bohr"
                if self.ibrav == 0:
                    self.alat = self.celldm[0]
                    self.cell_parameters_units = "alat"
            else:
                #((not celldm is None ) and (not crystal is None)):
                raise ValueError('celldm and crystallographic constant can not be specified together')
        else:
            if isinstance(cell_parameters_units, str):
                if (cell_parameters_units.lower() in ['alat', 'bohr', 'angstrom']):
                    self.cell_parameters_units = cell_parameters_units
                else:
                    raise ValueError, "call_parameters_units %s is unknown" %cell_parameters_units            

                
        self.nbands = nbands                                             # called nbands in QE [Q] how to define the default?
        if not nbands is None:
            self.set_inp_system('nbnd', nbands)
            self.nbands = nbands

        if ecutwfc is None:
            raise ValueError('kinetic energy cutoff (Ry) fpr wavefunction is required and not set')
        else:
            self.ecutwfc = ecutwfc
            self.set_inp_system('ecutwfc', ecutwfc)

        if not ecutrho is None:
            self.ecutrho = 4.0 * self.ecutwfc
            self.set_inp_system('ecutrho', ecutrho)

        if not smearing is None:
            self.set_inp_system('smearing', smearing)
            self.semaring = smearing

        if not degauss is None:
            self.set_inp_system('degauss', degauss)
            self.degauss=degauss

        if not input_dft is None:
           self.set_inp_system('input_dft', input_dft)
        
        if not mixing_mode is None:
            self.set_inp_system('mixing_mode', mixing_mode)

        if not occupations is None:
            self.set_inp_system('occupations', occupations)

        # | ELECTRONS input parameters
        
        if not electron_maxstep is None:
             self.set_inp_electrons('electrons_maxtep', electron_maxstep)

        if not conv_thr is None:
            self.set_inp_electrons('conv_thr', conv_thr)
            
        if not adaptive_thr is None:
            self.set_inp_electrons('adaptive_thr', adaptive_thr)

        if not conv_thr_init is None:
            if not hasattr(self, 'adaptive_thr'):
                raise ValueError('conv_thr_init can not be used in absence of adaptive_thr')
            if (hasattr(self, 'adaptive_thr') and not self.adaptive_thr):
                raise ValueError('conv_thr_init can not be used if adaptive_thr is false')
            set_inp_electrons('self.conv_thr_init', conv_thr_init)
            
        if not conv_thr_multi is None:
            if not hasattr(self, 'adaptive_thr'):
                raise ValueError('conv_thr_multi can not be used in absence of adaptive_thr')
            if (hasattr(self, 'adaptive_thr') and not self.adaptive_thr):
                raise ValueError('conv_thr_multi can not be used if adaptive_thr is false')
            self.set_inp_electrons('self.conv_thr_multi', conv_thr_multi)

        if not mixing_beta is None:
            self.set_inp_electrons('self.mixing_beta', mixing_beta)

        if not diag is None:
            self.set_inp_electrons('diagonalization', diag)

        if not spoint is None:
            self.set_inp_electrons('startingpoint', spoint)
        
        if not swfc is None:
            self.set_inp_electrons('startingwfc', swfc)

        # | CELL input parameters

        if not cell_dynamics is None:
            if ((self.calculation.lower() is 'vc-realx') and (cell_dynamics.lower() in ['none', 'sd', 'damp-pr', 'damp-w', 'bfgs'])) or \
                    ((self.calculation.lower() is 'vc-md') and (cell_dynamics.lower() in ['none', 'pr', 'w'])):
                self.set_inp_cell('cell_dynamics', cell_dynamics)
            if not press is None:
                self.set_inp_cell('press', press)
            if not wmass is None:
                self.set_inp_cell('wmass', wmass)            

        self.cell_press = press


        # | K_POINTS

        self.set_kpts(kpts, kpts_type)
            
        # | CELL_PARAMETERS

                # done here above, in SYSTEM section
            
        # | ------------------------------------------

        if 'ESPRESSO_OUTDIR' in os.environ:
            self.outdir = os.environ['ESPRESSO_OUTDIR']                                
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)
        else:
            self.outdir = './'
        if not os.access(self.outdir,os.W_OK):
                raise IOError, "you do not have write access to provided outdir"

        if pseudodir != None :
            self.pseudodir = pseudodir
        else :
            if 'ESPRESSO_PSEUDO' in os.environ:
                self.pseudodir = os.environ['ESPRESSO_PSEUDO']
            else:
                self.pseudodir = '$HOME/espresso/pseudo/'
            if not os.access(self.outdir,os.R_OK):
                    raise IOError, "you do not have write access to provided outdir"

        if 'ESPRESSO_ROOT' in os.environ:
            self.rootdir = os.environ['ESPRESSO_ROOT']
            self.bindir = self.rootdir + '/PW/src/'
        else:
            self.rootdir= '/usr/'
            self.bindir = self.rootdir + '/bin/'
            
        self.execfile = "pw.x"
        if not os.access(self.bindir+self.execfile, os.F_OK):
            raise IOError, "binary %s does not exist" % (self.bindir+self.execfile)
        if not os.access(self.bindir+self.execfile, os.X_OK):
            raise IOError, "you do not have execution permission on the binary %s" % self.bindir+self.execfile
            
        self.converged = False


    def initialize(self, atoms, if_pos=None):
        """Prepare the input files required to
        start the program (calculator).  """

        self.positions = atoms.get_positions().copy()
        # note: cell will be copied as it comes in the init file ONLY IF
        #       cell parameters are NOT given in input to the calculator
        self.cell = atoms.get_cell().copy()
        if (len(self.cell.shape) == 1):
            # it is an orthombic cell: make it 3x3 array
            _local= array=np.zeros((3,3), dtype=float)
            for i in range(3):
                _local[i][i] = self.cell[i]
            self.cell = _local.copy()
            if self.alat == 0:
                self.alat = self.cell[0][0]
            
        
        self.pbc = atoms.get_pbc().copy()

        if not if_pos is None:
            self.set_ifpos(if_pos)
        
        self.nat = atoms.get_number_of_atoms()  # how many atoms are in the atoms objext 
        self.set_inp_system('nat', self.nat)
        
        self.numbers = atoms.get_atomic_numbers().copy()
        
        self.species = []
        for a, Z in enumerate(self.numbers):
            if not Z in self.species:
                self.species.append(Z)
        self.ntyp = len(self.species)
        self.set_inp_system('ntyp', self.ntyp) #how many atom types are in the atoms object

        #NOTE: outdir, pseudodir and prefix have been initialized in __init__

        self.converged = False


        
# ============================================================================
# :
# : I/O and info methods
# :

    def get_name(self):
        """Return the name of the calculator (string)"""

        p=subprocess.Popen("echo $'\cc'| ./pw.x | grep Program", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err=p.communicate()
        print out.split(" ")[out.find("Program")+1]

            
    def get_version(self):
        """Return the version of the calculator (string)"""

        p=subprocess.Popen("echo $'\cc'| ./pw.x | grep Program", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err=p.communicate()
        print out.split(" ")[out.find("Program")+2]


    def get_cell(self):
        """Get the three unit cell vectors as a 3x3 ndarray."""
        return self.cell.copy()
    
    def write_infile(self):
        """Write input parameters to .in file"""
        fh = open(self.label + '.in', 'w')
        
        if os.access(self.outdir,os.W_OK):
            import getpass
            username=getpass.getuser()
            scratch = '/scratch/'+username
        elif os.access(os.curdir,os.W_OK):
            scratch = os.curdir
            print "Warning: outdir is set to current dir"
        else:
            raise IOError,"No write access neither to outdir nor to current dir"


        # ----------------------------------------------------------
        # | CONTROL section
        # ----------------------------------------------------------

        fh.write(' &CONTROL\n')
        fh.write("  calculation = '%s',\n" % self.calculation)

        #check for default cases
        if ((self.calculation in ["scf", "nscf", "bands"]) and (self.nstep > 1)) or \
                ((not self.calculation in ["scf", "nscf", "bands"]) and (self.nstep != 50)):
            fh.write('  nstep = %d,\n' % self.nstep)

        if len(self.inp_control.keys()) == 0:
            raise NotImplementedError, "CONTROL section appears to be empty: that is not allowed"
        
        self.write_section_params(self.inp_control, self.inp_control_writeall, fh)
        
        fh.write("  pseudo_dir = '%s',\n" % self.pseudodir)
        fh.write("  outdir = '%s',\n" % self.outdir)

        fh.write(' /\n')
        
        # ----------------------------------------------------------
        # | SYSTEM section
        # ----------------------------------------------------------

        if len(self.inp_system.keys()) == 0:
            raise NotImplementedError, "SYSTEM section appears to be empty: that is not allowed"
        
        fh.write(' &SYSTEM\n')
        self.write_section_params(self.inp_system, self.inp_system_writeall, fh)
        
        fh.write(' /\n')


        # ----------------------------------------------------------
        # | ELECTRONS section
        # ----------------------------------------------------------

        fh.write(' &ELECTRONS\n')
        self.write_section_params(self.inp_electrons, self.inp_electrons_writeall, fh)
                    
        fh.write(' /\n')


        # ----------------------------------------------------------
        # | IONS section
        # ----------------------------------------------------------

        fh.write(' &IONS\n')
        self.write_section_params(self.inp_ions, self.inp_ions_writeall, fh)
        
        fh.write(' /\n')


        # ----------------------------------------------------------
        # | CELL section
        # ----------------------------------------------------------

        if self.calculation in ['vc-relax', 'vc-md']:
            fh.write('&CELL\n')
            self.write_section_params(self.inp_cell, self.inp_cell_writeall, fh)
        
            fh.write('/\n')

        
        # ----------------------------------------------------------
        # | Card: ATOMIC_SPECIES
        # ----------------------------------------------------------
            
        fh.write('ATOMIC_SPECIES\n')

        PP = self.PP[:]
        PP_symbols = self.PP_symbols[:]
        occurrences = {}
        
        for Z in self.species:
            symbol = chemical_symbols[abs(Z)]
            mass = atomic_masses[abs(Z)]

            if not occurrences.has_key(symbol):
                occurrences[symbol] = 1
            else:
                occurrences[symbol] += 1

            if PP_symbols.count(symbol) == 0:
                raise ValueError, "you did not specify PP file for element %s" % symbol
            
            idx = PP_symbols.index(symbol)

            # write the correct entry that is formed by SYM MASS PP-FILE-NAME
            fh.write(' %s %f %s\n' % (symbol, mass, PP[idx][2]))
            
            # PP contains now the no-unique list of PP files
            # make it a unique list
            if(PP[idx][1] == occurrences[symbol]):
                PP.remove(symbol)
                PP_symbols.remove(symbol) 


        # ----------------------------------------------------------
        # | Card: ATOMIC_POSITIONS
        # ----------------------------------------------------------
        
        fh.write('ATOMIC_POSITIONS\n')

        if_pos = None
        if not self.if_pos is None:
            if_pos = self.if_pos[:]
            if_pos_symbols = self.if_pos_symbols[:]

        occurrences = {}
        
        for Z, pos in zip(self.numbers, self.positions):
            symbol = chemical_symbols[abs(Z)]
            if not occurrences.has_key(symbol):
                occurrences[symbol] = 1
            else:
                occurrences[symbol] += 1
            mystring = ' '.join("%f " % f for f in pos)
            fh.write(' %s %s' % (symbol, mystring))

            if not if_pos is None:
                if if_pos_symbols.count(symbol) > 0:
                    idx = if_pos_symbols.index(symbol)
                    if (if_pos[idx][1] == occurrences[symbol]) or (if_pos[idx][1] == 0):
                        mystring = ' '.join("%f " % f for f in if_pos[idx][2])
                        fh.write(' %s' % mystring)
                        if(if_pos[idx][1] == occurrences[symbol]):
                            if_pos.remove(symbol)
                            if_pos_symbols.remove(symbol)

            fh.write('\n')
                        

        # ----------------------------------------------------------
        # | Card: K_POINTS
        # ----------------------------------------------------------

        fh.write('K_POINTS %s\n' % self.kpts_type)
        if self.kpts_type is 'automatic':
            fh.write('   %d %d %d   %d %d %d\n' % tuple(self.kpts))
        elif not self.kpts_type is 'gamma':
            fh.write('  %d\n' % self.kpts_num)
            for array in self.kpts:
                fh.write('   %.14f %.14f %.14f %.14f\n' % tuple(array))
        

        # ----------------------------------------------------------
        # | Card: CELL_PARAMETERS
        # ----------------------------------------------------------

        if self.ibrav == 0:
            if (self.alat is None) and (self.cell_parameters_units.lower() is "alat"):
                raise ValueError, "it seems you want to specify CELL PARAMETERS in alat units, but you did not specify alat"

            if not self.cell_parameters_units.lower() in ["alat", "bohr", "angstrom"]:
                raise ValueError, "you specified un uknown type for CELL_PARAMETERS"
            
            fh.write('CELL_PARAMETERS {%s}\n' % self.cell_parameters_units)

            # Use self.cell ONLY IF ibrav == 0 AND cell_parameters have not been specified
            if self.cell_parameters is None:
                self.cell_parameters = self.cell

            if (self.cell_paramters.shape[0] != 3) or (self.cell_parameters.shape[1] != 3):
                raise ValueError, "cell parameters must be 3x3 "
            for i in range(3):
                fh.write('%f %f %f\n' % tuple(self.cell_parameters[i]))                        


        # ----------------------------------------------------------
        # | Card: OCCUPATIONS
        # ----------------------------------------------------------

                    
        if hasattr(self, 'occupations'):
            if self.occupations.lower() is 'from_input':
                fh.write('OCCUPATIONS\n')

                for i in range(self.nbands):
                    counter = 1
                    while counter < 10:
                        fh.write('%.14f ' % self.occupations1[i])
                        counter += 1

                    fh.write('\n')

                    if self.nspin == 2:
                        for i in range(self.nbands):
                            counter = 1
                            while counter < 10:
                                fh.write('%.14f ' % self.occupations2[i])
                                counter += 1

                            fh.write('\n')


        # ----------------------------------------------------------
        # | CLOSE IN FILE
        # ----------------------------------------------------------
                    
        fh.close()


    def write_section_params(self, inp, inp_writeall, fh):
        """ Automatically writes all the specified parameters for this section.
            Multidimensional parameters are taken into account properly.
        """
        for key, value in inp.items():
            if value is None:
                continue

            if isinstance(value, list):
                if inp_writeall[key] == 1:
                    for element in value:
                        fh.write('    %s(%d) = %s' % (key, value.index(element), str(value)))
                else:
                    for i in range(len(value)):
                        if inp_writeall[key][i] == 1:
                            fh.write('    %s(%d) = %s,\n' % (key, i+1, str(value[i])))   # i+1 instead of i because enumeration must start from 1

            elif isinstance(value, bool):
                fh.write('    %s = .%s.,\n' % (key, str(value).lower()) )
            elif isinstance(value, float):
                fh.write('    %s = %g,\n' % (key, value))
            elif isinstance(value, int):
                fh.write('    %s = %d,\n' % (key, value))
            elif isinstance(value, basestring):  
                fh.write("    %s = '%s',\n" % (key, value))

# END write_infile()


                
# ============================================================================
# :
# : Calculation methods
# :



    def update(self, atoms):                                         # [Q] what are the conditions for performing a new calculation on the system?
        if (not self.converged or                                    # not yet calculated
            len(self.numbers) != len(atoms) or                       # number of atoms has chenged
            (self.numbers != atoms.get_atomic_numbers()).any()):     # some atom has changed
            self.initialize(atoms)                                   
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or     # some atoms moved
              (self.pbc != atoms.get_pbc()).any()) :                 # pbc has changed
             # Not a good test self.cell and atoms.get_cell() are not properly handled yet 
             # (self.cell != atoms.get_cell().any())):                # configuration has changed 
            self.calculate(atoms)

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

#        if force_consistent:
#            return self.efree
#        else:
            # Energy extrapolated to zero Kelvin:
        return  (self.etotal)

    def get_potential_energy_accuracy(self, atoms, force_consistent=False):
        return self.etotal_accuracy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        '''Returns stress tensor (3x3)'''
        if (not self.tstress):
            self.stress = True
        self.update(atoms)
        return -np.array(self.stress)* 1e-1 * GPa

    def calculate(self, atoms):

        self.write_infile()

        infilename = self.label + '.in'
        outfilename = self.label + '.out'
        #rename the output file if it does already exist
        if islink(outfilename) or isfile(outfilename):
            os.rename(outfilename, outfilename+'.bak')

        _commandstring = '%s < %s > %s' % ( self.bindir+self.execfile, infilename, outfilename)
        exitcode = os.system(_commandstring)
        
        # if exitcode != 0:
        #     raise RuntimeError((self.execfile+' exited with exit code: %d.  ' +
        #                         'Check %s.log for more information.') %
        #                        (exitcode, self.label))

        print "calculation done, read output.."
        results = self.read_text_output()             # read in results from output
        #print "Keys obtained:", results.keys()
        self.nelect = results['nelec']
        self.nbnd = results['nbands']
        self.kpts_num = results['kpts_num']
        self.xc = results['exchange_correlation']
        self.niter = results['niter']
        self.pressure = results['pressure']
        self.stress = results['stress']
        self.fermi_energy = results['fermi_energy']
        self.kpts = results['kpts']
        self.kpts_wk = results['kpts_wk']
        self.atoms_forces = results['atoms_forces']
        self.cell = results['cell']
        self.atomic_positions = results['atomic_positions']
        self.etotal = results['etotal']
        self.etotal_accuracy = results['etotal_accuracy']
        self.fermi_energy = results['fermi_energy']
        self.total_magnetization = results['total_magnetization']
        self.absolute_magnetization = results['absolute_magnetization']
        self.alat = results['alat']
        
        
        if results is None:
            print "an error has occurred"
            print "%s" % (''.join(s+'\n' for s in self.error_message))
            self.converged = False
        else:
            self.converged = True


        
# ============================================================================
# :
# : Read methods
# :

    
    def read_text_output(self):
        from ase.io.qe import read_quantumespresso_textoutput
        filename = self.label + '.out'
        return read_quantumespresso_textoutput(filename)

    
    def read_xml_output(self):
        from ase.io.qe import read_quantumespresso_xmloutput
        filename = self.label + '.save/data-file.xml'        
        return read_quantumespresso_xmloutput(filename, "all")


    def get_number_of_iterations(self):
        return self.niter


    def get_number_of_electrons(self):
        return self.nelect

    
    def set_atoms(self, atoms):
        # atoms.pbc = [True, True, True]
        self.atoms = atoms.copy()

        
    def get_atoms(self):
        atoms = self.atoms.copy()
        # atoms.set_calculator(self)
        return atoms


    
# =====================================================================================================
#  routines that set up internal variables, permitted values, relationships among variables and so on
# =====================================================================================================

    
    def set_inp_control(self, key, value, write_all=1):
        """Sets INPUT parameter for CONTROL section.
           If the key is a list or a dictionary AND write_all is zero,
           only the non-empty/non-zero values will be written in .in file
        """
        self.inp_control[key] = value
        if isinstance(value, list) or isinstance(value, dict):
            self.inp_control_writeall[key] = write_all

            
    def set_inp_system(self, key, value, write_all=1):
        """Sets INPUT parameter for SYSTEM section.
           If the key is a list or a dictionary AND write_all is zero,
           only the non-empty/non-zero values will be written in .in file
        """
        self.inp_system[key] = value
        if isinstance(value, list) or isinstance(value, dict):
            self.inp_system_writeall[key] = write_all


    def set_inp_electrons(self, key, value, write_all=1):
        """Sets INPUT parameter for ELECTRONS section.
           If the key is a list or a dictionary AND write_all is zero,
           only the non-empty/non-zero values will be written in .in file
        """
        self.inp_electrons[key] = value
        if isinstance(value, list) or isinstance(value, dict):
            self.inp_electrons_writeall[key] = write_all

        
    def set_inp_ions(self, key, value, write_all=1):
        """Sets INPUT parameter for IONS section.
           If the key is a list or a dictionary AND write_all is zero,
           only the non-empty/non-zero values will be written in .in file        
        """
        self.inp_ions[key] = value
        if isinstance(value, list) or isinstance(value, dict):
            self.inp_ions_writeall[key] = write_all

            
    def set_inp_cell(self, key, value, write_all=1):
        """Sets INPUT parameter for CELL section.
           If the key is a list or a dictionary AND write_all is zero,
           only the non-empty/non-zero values will be written in .in file        
        """
        self.inp_cell[key] = value
        if isinstance(value, list) or isinstance(value, dict):
            self.inp_cell_writeall[key] = write_all


    def set_if_pos(self, if_pos):
        """Sets the values for the if_pos arrays in the ATOMIC_POSITIONS CARD.
           if_pos is meant to be a list of tuples, having an entry for each atoms you
           want to put a force on.
           Each tuple is formed as
                [Symbol, occurrence, force]
           where 'Symbol' is the chemical symbol of the atom you are referring to, 'occurrence'
           is the occurrence of that atoms, 'force' is 3-component list of integers that
           actually represents the force.

           note: the occurrence is meant to be in fortran-style, i.e. to start from 1 not from 0.
                 A 0 value means "all occurrences"
           
        """
        if isinstance(if_pos, list):
            # check that if_pos object is well-formed
            for el in if_pos:
                if (isinstance(el[0], str) and isinstance(el[1], int) and isinstance(el[2], list)):
                    for subel in el[2]:
                        if not isinstance(el[2][subel], int):
                            raise TypeError, "if_pos object is malformed"
                else:
                    raise TypeError, "if_pos object is malformed"
                
            # sort the list by chemical symbol as primary key and occurrence as secondary key, so
            # that the imposed forces may be retrieved in order of appearence of elements
            if_pos.sort(key=lambda l: (l[0], l[1]))
            self.if_pos = if_pos[:]
            self.if_pos_symbols = []
            for el in self.if_pos:
                self.if_pos_symbols.append[if_pos[el][0]]
                
        else:
            raise TypeError, "if_pos object is malformed"

                                           
    def set_PP(self, PP):
        """Sets the values for the PseudoPotential file names in the ATOMIC_SPECIES CARD.
           PP is meant to be a list of tuples, having an entry for each atom (or atom specie,
           see below) is in the system.
           Each tuple is formed as
                [Symbol, occurrence, filename]
           where 'Symbol' is the chemical symbol of the atom you are referring to, 'occurrence'
           is the occurrence of that atoms, 'filename' is a string that contains the file name
           of the pseudopotential file

           note: the occurrence is meant to be in fortran-style, i.e. to start from 1 not from 0
                 A 0 value means "all occurrences"          
        """
        if isinstance(PP, list):
            # check that PP object is well-formed
            for el in PP:
                if (not isinstance(el[0], str)) or (not isinstance(el[1], int)) or (not isinstance(el[2], str)):
                    raise TypeError, "PP object is malformed"
                
            # sort the list by chemical symbol as primary key and occurrence as secondary key, so
            # that the pseudo potential may be retrieved in order of appearence of elements
            PP.sort(key=lambda l: (l[0], l[1]))
            self.PP = PP[:]
            self.PP_symbols = []
            for el in PP:
                self.PP_symbols.append(el[0])
        else:
            raise TypeError, "PP object is malformed"
                                           
            
    def set_ibrav_values(self):
        """build the list of allowed values for ibrav
           build the dictionary of what entries are expected to be specified for the celldm()
           array as a function of ibrav value
        """
        # build the list of permitted values for ibrav
        self.ibrav_allowed = range(15)
        others=[-5, -9, -12]
        self.ibrav_allowed.append(others)

        # build the 'expectation list' of values for celldm given an ibrav value
        self.celldm_expected_values ={0: [1,0,0,0,0,0],
                                      1: [1,0,0,0,0,0],
                                      2: [1,0,0,0,0,0],
                                      3: [1,0,0,0,0,0],
                                      4: [1,0,1,0,0,0],
                                      5: [1,0,0,1,0,0],
                                      -5: [1,0,0,1,0,0],
                                      6: [1,0,1,0,0,0],
                                      7: [1,0,1,0,0,0],
                                      8: [1,1,1,0,0,0],
                                      9: [1,1,1,0,0,0],
                                      -9: [1,1,1,0,0,0],
                                      10: [1,1,1,0,0,0],
                                      11: [1,1,1,0,0,0],
                                      12: [1,1,1,1,0,0],
                                      -12: [1,1,1,0,1,0],
                                      13: [1,1,1,1,0,0],
                                      14: [1,1,1,1,1,1]}
        
    def set_celldm(self, celldm, crystal):
        """check that the given celldm array is consistent with the specified ibrav value and
           assign its the values
           note that celldm is an array that contains only the value you need to assign. For
           instance, let's say ibrav=5: then you need to set only the first and the fourth entries
           and the celldm list that you pass to the calculator might be a 2-entries list and
           does not need to be a 6-entries list with 4 zero entries.
           Of course you may pass as well a full 6-entries list.
           If crystal is given, uses crystal.
        """
        if (celldm is None) and (not crystal is None):
            self.celldm = None
            self.crystal = crystal[:]
            # adds input parameters to the system section for the purpose of automaticly writing them
            self.set_inp_system('A', self.crystal[0])
            self.set_inp_system('B', self.crystal[1])
            self.set_inp_system('C', self.crystal[2])
            self.set_inp_system('cosAB', self.crystal[3])
            self.set_inp_system('cosAC', self.crystal[4])
            self.set_inp_system('cosBC', self.crystal[5])
            self.alat = self.crystal[0]
            
        else:
            # checks whether celldm is a list
            if not isinstance(celldm, list):
                if not isinstance(celldm, float):
                    raise TypeError, "celldm must be a list of floats or at least a single float"
                # if it is a float, makes it a one-entry list
                else:
                    celldm = [celldm]

            # checks whether celldm has as many entries as it is expected accordingly with the ibrav value    
            if len(np.nonzero(celldm)[0]) != self.celldm_expected_values[self.ibrav].count(1) :
                raise ValueError, "the number specified entries for celldm (%d) is different than expected for the given ibrav"

            # now insert the values at the right positions
            self.celldm = [0, 0, 0, 0, 0, 0]
            ctrl = self.celldm_expected_values[self.ibrav]
            idxs = [i for i in range(6) if ctrl[i]==1]
        
            for i in range(len(idxs)):
                self.celldm[idxs[i]]= celldm[i]
        
            self.set_inp_system('celldm', self.celldm, self.celldm_expected_values[self.ibrav])    
            self.alat = self.celldm[0]

            # note: celldm values will be automatically written as "celldm(i) = value" for all, and only, the needed i


    def check_kpts_type(self, kpts_type):
        """Checks whether the given type is among the allowed ones"""

        if kpts_type.lower() in ['gamma', 'automatic', 'tpiba', 'crystal', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c']:
            return True
        else:
            return False

        
    def set_kpts(self, kpts=None, kpts_type=None):
        """Sets the kpts accordingly with kpts_type """
        
        if not kpts_type is None:
            # if kpts_type is given and is valid, it overrides self.kpts_type
            if self.check_kpts_type(kpts_type):
                self.kpts_type = kpts_type
        else:
            # otherwise, check whether it already exists in self
            if self.kpts_type is None:
                raise ValueError, "no kpts_type has been defined"
            else:
                kpts_type = self.kpts_type

        if kpts_type.lower() == 'gamma':
            self.kpts = None

        elif kpts_type.lower() == "automatic":
            if isinstance(kpts, list):
                if (len(kpts) == 6):
                    self.kpts = kpts[:]
                else:
                    raise TypeError, '6 values expected for kpts with \'automatic\' flag'
                
            else:
                raise TypeError, "wrong type for kpts"

        else:
            entries = 0;
            if isinstance(kpts, list):
                kpts_num = kpts[0]
                self.kpts = [];
                for el in range(len(kpts)):
                    if el > 0:
                        if isinstance(kpts[el], list) and (len(kpts[el]) == 4):
                            self.kpts.append(kpts[el])
                            entries += 1

                if entries != kpts_num:
                    raise ValueError, "the specified kpoints are malformed, %d found vs %d expected" % (entries, kpts_num)

                self.kpts_num = kpts_num

    def set_occupations(self, occ1, occ2=None):
        """sets the occupations array"""

        if not self.nbands is None:
            if len(occ1) == self.nbands:
                self.occupations1 = occ1[:]
                if (not occ2 is None) and (self.nspin == 2) and (len(occ2) == self.nbands):
                    self.occupations2 = occ2[:]
        else:
            raise ValueError, "impossibile to set up occupations, nbands has not been specified yet"
                    



# =====================================================================================================
#  POST - PROCESSING DRIVING ROUTINES
# =====================================================================================================

    def PP(self, plot_num, input_args, plot_args):
        """
        """

        # specifies what parameters must be present in input_args, depending on the value of plot_num
        input_deps = { '0' : ['spin_component'],
                       '1' : ['spin_component'],
                       '5' : ['sample_bias'],
                       '7' : ['kpoint', 'kband', 'lsign', 'spin_component'],
                       '10' : ['emin', 'emax', 'spin_component'],
                       '13' : ['spin_component'],
                       '17' : ['spin_component'] }

        # specifies what parameters must be present in plot_args, depending on the value of iflag
        # iflag is an entry of the plot_args dictionary
        plot_deps = { '0' : ['e1', 'x0', 'nx'],
                      '1' : ['e1', 'x0', 'nx'],
                      '2' : ['e1', 'e2', 'x0', 'nx', 'ny'],
                      '3' : ['e1', 'e2', 'e3', 'x0', 'nx,', 'ny', 'nz'] }
        
                      
        
        fh = open('inputpp.in', 'w')

        fh.write(' &INPUTPP\n')
        fh.write('  prefix = %s,\n' % self.prefix)
        fh.write('  outdir = %s,\n' % self.outdir)
        if not input_args_has_key('filename'):
            fh.write('  filplot = filplot,\n')
        else:
            fh.write('  filplot = %s,\n' % input_args['filename'])
        fh.write('  plot_num = %s,\n' % str(plot_num))
        
        for D in input_deps[str(plot_num)]:
            fh.write('  %s = %s,\n' % (D, str(input_args[D])))

        fh.write(' /\n')

        fh.write(' &PLOT\n')
        
        fh.write('  nfile = %s,\n' % str(plot_args['nfiles']))
        
        if plot_args['nfiles'] > 1:
            if not plot_args.has_key('filenames'):
                raise ValueError, "for more than 1 input file you must specify individual filename through plot_args"

            for i in range(plot_args['nfiles']):
                fh.write('  filepp('+str(i+1)+') = %s,\n' % plot_args['filenames'][i])
        else:
            if plot_args.has_key('filenames'):
                fh.write('  filepp(1) = %s,\n' % plot_args['filenames'])
            else:
                fh.write('  filepp(1) = filplot,\n')

        if plot_args.has_key('weights'):
            for i in range(plot_args['nfiles']):
                fh.write('  weight('+str(i+1)+') = %s,\n' % str(plot_args['weights'][i]))


        fh.write('  iflag = %s,\n' % str(plot_args['iflag']))
        fh.write('  output_format = %s,\n' % str(plot_args['formats']))
        fh.write('  fileout = %s,\n' % str(plot_args['output']))
        
        for D in plot_deps[str(plot_args['iflag'])]:

            value = plot_args[D]
            if isinstance(value, list):
                for el in value:
                    fh.write('  %s(%d) = %s,\n' % (D, value.index(el) ,str(value)))
            else:
                fh.write('  %s  = %s,\n' % (D, str(value)))
        
        fh.write(' /\n')
        execbin = "pp.x"
        if 'ESPRESSO_ROOT' in os.environ:
            rootdir = os.environ['ESPRESSO_ROOT']
            bindir = rootdir + '/PP/src/'
        else:
            rootdir= '$HOME/espresso/'
            bindir = rootdir + '/PP/src/'

        if not os.access(bindir+execbin, os.F_OK):
            raise IOError, "binary %s does not exist" % (bindir+execbin)
        if not os.access(bindir+execbin, os.X_OK):
            raise IOError, "you do not have execution permission on the binary %s" % bindir+execbin

        infilename='inputpp.in'
        outfilename='outputpp'
        _commandstring = '%s < %s > %s' % ( bindir+execfile, infilename, outfilename)
        exitcode = os.system(_commandstring)





# conversions among energy units
# all units are given in eV
units_energy = {'ev': 1,                          # ev
                'hr': 27.21138505,                # hartree
                'ry': 13.605692525,               # rydberg
                'j': 6.24150934731e+18}           # joule

# conversions among length units
# all units are given in bohr
units_length = {'bohr': 1,
                'a': 1.89035916824,               # angstrom
                'm': 1.89035916824e10}            # meter

units = {'ev': units_energy,
         'hr': units_energy,
         'ry': units_energy,
         'j' : units_energy,
         'bohr': units_length,
         'a': units_length,
         'm': units_length}

def convert(v1, v2):
    try:
        isinstance(v1, str) & isinstance(v2, str)
    except TypeError:
        print "convert argument must be strings"
    myv1=v1.lower()
    myv2=v2.lower()
    if (v1 in units.keys()) & (v2 in units.keys()):
        try:
            units[v1] == units[v2]
        except ValueError:
            print "units dimensions must be the same for conversion"

        return (units[v1])[v1] / (units[v2])[v2]
    else:
        return -1
