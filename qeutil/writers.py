
import os
from tempfile import mkdtemp
from pyspglib import spglib
from numpy.linalg import norm
from math import sqrt
from ase import Atoms

# PWscf input file
pw_in='''
&CONTROL
  calculation = '%(calc)s',
    prefix = '%(prefix)s',
    tstress = .true.,
  pseudo_dir = '../pspot',
  outdir = './%(outdir)s/',
/
&SYSTEM
    ecutwfc = %(ecut)f,
    ibrav = %(ibrav)d,
    nat = %(nat)d,
    ntyp = %(ntyp)d,
/
&ELECTRONS
/
ATOMIC_SPECIES
    %(at_species)s
ATOMIC_POSITIONS
    %(at_positions)s
K_POINTS automatic
   %(kx)d %(ky)d %(kz)d   %(shift)d %(shift)d %(shift)d
'''

# PH input file
ph_in='''
%(calc)s
&INPUTPH
    prefix = '%(prefix)s',
    outdir = './%(outdir)s/',
    ldisp = .true.,
    fildrho = 'fildrho',
    fildvscf = 'fildvscf',
    nq1 = %(nq)d , nq2 = %(nq)d , nq3 = %(nq)d
 /
'''

q2r_in='''
&INPUT
    fildyn='matdyn',
    zasr='%(asr)s',
    flfrc='%(prefix)s.fc'
/
'''

matdyn_in='''
&input
    asr='%(asr)s',  
    flfrc='%(prefix)s.fc', 
    flfrq='%(prefix)s.freq',
    q_in_band_form=.true.
/
'''

phdos_in='''
&input
    dos=.true.
    asr='%(asr)s',  
    flfrc='%(prefix)s.fc', 
    fldos='%(prefix)s.dos',
    nk1=%(doskx)d,
    nk2=%(dosky)d,
    nk3=%(doskz)d,
    ndos=%(ndos)d
 /
'''




def write_section_params(fh, k, p):
    """ Automatically writes all the specified parameters for this section.
        Multidimensional parameters are taken into account properly.
    """
    for key, value in p.items():
        # Nothing set skip
        if value is None:
            continue
        # Not this section skip
        if not key in k :
            continue

        # Write the rest
        if isinstance(value, list):
            if inp_writeall[key] == 1:
                for element in value:
                    fh.write('    %s(%d) = %s' % (key, value.index(element), str(value)))
            else:
                for i in range(len(value)):
                    if inp_writeall[key][i] == 1:
                        # i+1 instead of i because enumeration must start from 1
                        fh.write('    %s(%d) = %s,\n' % (key, i+1, str(value[i])))   

        elif isinstance(value, bool):
            fh.write('    %s = .%s.,\n' % (key, str(value).lower()) )
        elif isinstance(value, float):
            fh.write('    %s = %g,\n' % (key, value))
        elif isinstance(value, int):
            fh.write('    %s = %d,\n' % (key, value))
        elif isinstance(value, basestring):  
            fh.write("    %s = '%s',\n" % (key, value))

bravais_lattice_to_ibrav = {"free":0,
                            "cubic p (sc)": 1,
                            "cubic f (fcc)": 2,
                            "cubic i (bcc)": 3,
                            "hexagonal and trigonal p": 4,
                            "trigonal r": 5,
                            "tetragonal p (st)": 6,
                            "tetragonal i (bct)": 7,
                            "orthorhombic p": 8,
                            "orthorhombic base-centered(bco)": 9,
                            "orthorhombic face-centered": 10,
                            "orthorhombic body-centered": 11,
                            "monoclinic p": 12,
                            "monoclinic base-centered": 13,
                            "triclinic p": 14}


def write_cell_params(fh, a, p):
    '''
    Write the specification of the cell using a traditional A,B,C, 
    alpha, beta, gamma scheme. The symmetry and parameters are determined
    from the atoms (a) object. The atoms object is translated into a 
    primitive unit cell but *not* converted. This is just an internal procedure.
    If you wish to work with primitive unit cells in ASE, you need to make 
    a conversion yourself. The cell params are in angstrom.
    
    Input
    -----
        fh  file handle of the opened pw.in file
        a   atoms object
        p   parameters dictionary
        
    Output
    ------
        Primitive cell tuple of arrays: (lattice, atoms, atomic numbers)
    '''
    
    assert(p['use_symmetry'])
    
    # Get a spacegroup name and number for the a
    sg,sgn=spglib.get_spacegroup(a).split()
    # Extract the number
    sgn=int(sgn[1:-1])
    # Extract the lattice type
    ltyp=sg[0]
    
    # Find a primitive unit cell for the system
    # puc is a tuple of (lattice, atoms, atomic numbers)
    puc=spglib.find_primitive(a)
    cell=puc[0]
    apos=puc[1]
    anum=puc[2]

    # Select appropriate ibrav
    if sgn >= 195 :
        # Cubic lattice
        if   ltyp=='P':  p['ibrav']=1  # Primitive
        elif ltyp=='F':  p['ibrav']=2  # FCC
        elif ltyp=='I':  p['ibrav']=3  # FCC
        else :
            print 'Impossible lattice symmetry! Contact the author!'
            raise NotImplementedError
        fh.write('      A = %f,\n' % (sqrt(2)*norm(cell[0]),))
    elif sgn >= 143 :
        # Hexagonal and trigonal 
        if   ltyp=='P' :  p['ibrav']=4  # Primitive
        elif ltyp=='R' :  p['ibrav']=5  # Trigonal rhombohedral
        else :
            print 'Impossible lattice symmetry! Contact the author!'
            raise NotImplementedError
    elif sgn >= 75 :
        raise NotImplementedError
    else :
        raise NotImplementedError
    return puc


def write_pw_in(d,a,p):

    fh=open(os.path.join(d,'pw.in'),'w')
    # ----------------------------------------------------------
    # CONTROL section
    # ----------------------------------------------------------

    fh.write(' &CONTROL\n')
    fh.write("    calculation = '%s',\n" % p['calc'])

    pwin_k=['tstress','nstep','pseudo_dir','outdir','wfcdir']
    write_section_params(fh, pwin_k, p)
    fh.write(' /\n')
    
    # ----------------------------------------------------------
    # SYSTEM section
    # ----------------------------------------------------------

    fh.write(' &SYSTEM\n')
    if p['use_symmetry'] :
        # Need to use symmetry properly
        # create a dummy atoms object for primitive cell
        primcell=write_cell_params(fh,a,p)
        cr=Atoms(cell=primcell[0],scaled_positions=primcell[1],numbers=primcell[2],pbc=True)
    else :
        cr=a
    fh.write("    nat = %d\n" % (cr.get_number_of_atoms()))
    fh.write("    ntyp = %d\n" % (len(set(cr.get_atomic_numbers()))))
    pwin_k=['ecutwfc','ibrav']
    write_section_params(fh, pwin_k, p)
    fh.write(' /\n')

    # ----------------------------------------------------------
    # ELECTRONS section
    # ----------------------------------------------------------

    fh.write(' &ELECTRONS\n')
    fh.write(' /\n')

    
    # ----------------------------------------------------------
    # Card: ATOMIC_SPECIES
    # ----------------------------------------------------------
        
    fh.write('ATOMIC_SPECIES\n')
    
    xc=p['xc']
    pp_type=p['pp_type']
    pp_format=p['pp_format']
    for nm, mass in set(zip(cr.get_chemical_symbols(),cr.get_masses())):
        fh.write("    %s %g %s_%s_%s.%s \n" % (nm, mass, nm, xc, pp_type, pp_format))

    # ----------------------------------------------------------
    # Card: ATOMIC_POSITIONS
    # ----------------------------------------------------------
    
    fh.write('ATOMIC_POSITIONS crystal\n')
    # Now we can write it out
    for n,v in zip(cr.get_chemical_symbols(),cr.get_scaled_positions()):
        fh.write("    %s %g %g %g\n" % (n, v[0], v[1], v[2]))

    # ----------------------------------------------------------
    # Card: CELL_PARAMETERS
    # ----------------------------------------------------------
    
    fh.write('CELL_PARAMETERS angstrom\n')
    
    for v in cr.get_cell():
        fh.write('    %f %f %f\n' % tuple(v))
        
    # ----------------------------------------------------------
    # Card: K_POINTS
    # ----------------------------------------------------------

    fh.write('K_POINTS %s\n' % p['kpt_type'])
    if p['kpt_type'] is 'automatic':
        fh.write('   %d %d %d   %d %d %d\n' % (tuple(p['kpts'])+tuple(p['kpt_shift'])))
    elif not p['kpt_type'] is 'gamma':
        fh.write('  %d\n' % len(p['kpts']))
        for v in p['kpts']:
            fh.write('   %.14f %.14f %.14f %.14f\n' % tuple(v))
            
    fh.close()
    

def make_calc_dir(prefix, bdir='./'):
    '''
    Builds a calculation directory in the given base directory.
    Creates a unique calculation directory and builds all ".in" 
    files for basic implemented functionality.
    Prepares ".in" files for the following calculations:
    
        - pw.x scf calculation
        - ph.x phonon calculation
        - q2r.x real-space force constants calculation
        - matdyn.x phonon dispersion calculation
        - matdyn.x phonon density of states calculation
    
    INPUT
    -----
        params - calculation parameters dictionary
        bdir - base directory
    
    OUTPUT
    ------
        cdir - name of the created calculation directory
    
    '''
    
    bcdir=mkdtemp(prefix=prefix + '.',dir=bdir)
    #cdir=os.path.split(bcdir)[1]
    #print cdir, bcdir
    
#    open(os.path.join(bcdir,'ph.in'),'w').write(ph_in % params)
#    open(os.path.join(bcdir,'q2r.in'),'w').write(q2r_in % params)
#    f=open(os.path.join(bcdir,'matdyn.in'),'w')
#    f.write(matdyn_in % params)
#    f.write('%d\n'%qp.shape[0])
#    for v in qp:
#        f.write("%f %f %f %d\n" % (v[0],v[1],v[2],params['points']))
#    f.close()
#    open(os.path.join(bcdir,'phdos.in'),'w').write(phdos_in % params)
    return bcdir

def make_D3_files(bdir,clc,cdat):
    open(os.path.join(bdir,clc,'phD3G.in'),'w').write(ph3G_in % cdat)
    open(os.path.join(bdir,clc,'D3G.in'),'w').write(d3G_in % cdat)
    open(os.path.join(bdir,clc,('phD3_%(sufix)s.in' % cdat)),'w').write(ph3any_in % cdat)
    open(os.path.join(bdir,clc,('D3_%(sufix)s.in' % cdat)),'w').write(d3any_in % cdat)
    
    
    

