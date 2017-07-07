from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons

# The qe-util package
from qeutil import QuantumEspresso

# Setup crystal and EMT calculator
atoms = bulk('Al', 'fcc', a=4.05)

qe=QuantumEspresso(label='Al',                    # Label for calculations
                       wdir='.',                   # Working directory
                       pseudo_dir='/usr/share/espresso/pseudo',   # Directory with pseudopotentials
                       kpts=[2,2,2],               # K-space sampling for the SCF calculation
                       xc='pz',                    # Exchange functional type in the name of the pseudopotentials
                       pp_type='vbc',              # Variant of the pseudopotential
                       pp_format='UPF',            # Format of the pseudopotential files
                       ecutwfc=15,                 # Energy cut-off (in Rydberg)
                       occupations='smearing', 
                       smearing='methfessel-paxton', 
                       degauss=0.0198529,
                       mixing_beta=0.7, mixing_mode = 'plain', diagonalization = 'cg',
                       restart_mode='from_scratch', tstress=False,
                       use_symmetry=False, procs=8)          # Use symmetry in the calculation ?
                       
print qe.directory

# Assign the calculator to our system
atoms.set_calculator(qe)

ph = Phonons(atoms, qe, supercell=(2,2,2), delta=0.05)
ph.run()
# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)

print 'All done and functional: QE, phonons.'
