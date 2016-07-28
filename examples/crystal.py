# Import the basic libraries

# ASE system
import ase
from ase import Atom, Atoms
from ase import io
from ase.lattice.spacegroup import crystal

# Spacegroup/symmetry library
from pyspglib import spglib

# The qe-util package
from qeutil import QuantumEspresso

# iPython utility function
from IPython.core.display import Image

# Configure qe-util for local execution of the Quantum Espresso on four processors
# QuantumEspresso.pw_cmd='mpiexec -n 4 pw.x < pw.in > pw.out'

a=4.3596                                         # Lattice constant in Angstrom
cryst = crystal(['Si', 'C'],                     # Atoms in the crystal
                [(0, 0, 0), (0.25, 0.25, 0.25)], # Atomic positions (fractional coordinates)
                spacegroup=216,                  # International number of the spacegroup of the crystal
                cellpar=[a, a, a, 90, 90, 90])   # Unit cell (a, b, c, alpha, beta, gamma) in Angstrom, Degrees
                
# Write the image to disk file
ase.io.write('crystal.png',       # The file where the picture get stored
             cryst,               # The object holding the crystal definition
             format='png',        # Format of the file
             show_unit_cell=2,    # Draw the unit cell boundaries
             rotation='115y,15x', # Rotate the scene by 115deg around Y axis and 15deg around X axis
             scale=30)            # Scale of the picture

# Display the image
Image(filename='crystal.png')

print 'Space group:', spglib.get_spacegroup(cryst)

# Create a Quantum Espresso calculator for our work. 
# This object encapsulates all parameters of the calculation, 
# not the system we are investigating.
qe=QuantumEspresso(label='SiC',                    # Label for calculations
                       wdir='.',                   # Working directory
                       pseudo_dir='/usr/share/espresso/pseudo',   # Directory with pseudopotentials
                       kpts=[2,2,2],               # K-space sampling for the SCF calculation
                       xc='pz',                    # Exchange functional type in the name of the pseudopotentials
                       pp_type='vbc',              # Variant of the pseudopotential
                       pp_format='UPF',            # Format of the pseudopotential files
                       ecutwfc=20,                 # Energy cut-off (in Rydberg)
                       use_symmetry=False, procs=8)          # Use symmetry in the calculation ?
                       
print qe.directory

# Assign the calculator to our system
cryst.set_calculator(qe)

# Run the calculation to get stress tensor (Voigt notation, in eV/A^3) and pressure (in kBar)
print "Stress tensor   (Voigt notation eV/A^3):", cryst.get_stress()
# print "External pressure                (kBar):", cryst.get_isotropic_pressure(cryst.get_stress())*1e-3

print 'All done and functional: QE, stress.'
