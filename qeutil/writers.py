
from __future__ import division

import os
from tempfile import mkdtemp
from pyspglib import spglib
from numpy.linalg import norm
from numpy import array
from math import sqrt
from ase import Atoms

# PWscf input file
# not used - just a comment
#pw_in='''
#&CONTROL
#  calculation = '%(calc)s',
#    prefix = '%(prefix)s',
#    tstress = .true.,
#    tprnfor = .true.,
#    verbosity = 'high',
#  pseudo_dir = '../pspot',
#  outdir = './%(outdir)s/',
#/
#&SYSTEM
#    ecutwfc = %(ecut)f,
#    ibrav = %(ibrav)d,
#    nat = %(nat)d,
#    ntyp = %(ntyp)d,
#/
#&ELECTRONS
#/
#ATOMIC_SPECIES
#    %(at_species)s
#ATOMIC_POSITIONS
#    %(at_positions)s
#K_POINTS automatic
#   %(kx)d %(ky)d %(kz)d   %(shift)d %(shift)d %(shift)d
#'''





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
    icell=a.get_cell()
    A=norm(icell[0])
    B=norm(icell[1])
    C=norm(icell[2])

    # Select appropriate ibrav
    if sgn >= 195 :
        # Cubic lattice
        if   ltyp=='P':  
            p['ibrav']=1  # Primitive
            qepc=array([[1,0,0],[0,1,0],[0,0,1]])
        elif ltyp=='F':  
            p['ibrav']=2  # FCC
            qepc=array([[-1,0,1],[0,1,1],[-1,1,0]])/2.0
        elif ltyp=='I':  
            p['ibrav']=3  # BCC
            qepc=array([[1,1,1],[-1,1,1],[-1,-1,1]])/2.0
        else :
            print 'Impossible lattice symmetry! Contact the author!'
            raise NotImplementedError
        #a=sqrt(2)*norm(cell[0])
        qepc=A*qepc
        fh.write('      A = %f,\n' % (A,))
    elif sgn >= 143 :
        # Hexagonal and trigonal 
        if   ltyp=='P' :  
            p['ibrav']=4  # Primitive
            qepc=array([[1,0,0],[-1/2,sqrt(3)/2,0],[0,0,C/A]])
        elif ltyp=='R' :  
            p['ibrav']=5  # Trigonal rhombohedral
            raise NotImplementedError
        else :
            print 'Impossible lattice symmetry! Contact the author!'
            raise NotImplementedError
        qepc=A*qepc
        fh.write('      A = %f,\n' % (A,))
        fh.write('      C = %f,\n' % (C,))
    elif sgn >= 75 :
        raise NotImplementedError
    elif sgn ==1 :
        # P1 symmetry - no special primitive cell signal to the caller
        p['ibrav']=0
        return None
    else :
        raise NotImplementedError
    cp=Atoms(cell=puc[0], scaled_positions=puc[1], numbers=puc[2], pbc=True)
    qepc=Atoms(cell=qepc, 
        positions=cp.get_positions(), 
        numbers=cp.get_atomic_numbers(), pbc=True)
    return qepc.get_cell(), qepc.get_scaled_positions(), qepc.get_atomic_numbers()


def write_pw_in(d,a,p):

    if 'iofile' in p.keys() :
        # write special varians of the pw input
        fh=open(os.path.join(d,p['iofile']+'.in'),'w')
    else :
        # write standard pw input
        fh=open(os.path.join(d,'pw.in'),'w')

    # ----------------------------------------------------------
    # CONTROL section
    # ----------------------------------------------------------

    fh.write(' &CONTROL\n')
    fh.write("    calculation = '%s',\n" % p['calc'])

    pwin_k=['tstress', 'tprnfor','nstep','pseudo_dir','outdir',
            'wfcdir', 'prefix','forc_conv_thr', 'etot_conv_thr']
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
        if primcell :
            # primitive cell has been found - let us use it
            cr=Atoms(cell=primcell[0],scaled_positions=primcell[1],numbers=primcell[2],pbc=True)
        else :
            # no primitive cell found - drop the symmetry
            cr=a
    else :
        cr=a
        p['ibrav']=0
    fh.write("    nat = %d,\n" % (cr.get_number_of_atoms()))
    fh.write("    ntyp = %d,\n" % (len(set(cr.get_atomic_numbers()))))
    pwin_k=['ecutwfc','ibrav','nbnd','occupations']
    write_section_params(fh, pwin_k, p)
    fh.write(' /\n')

    # ----------------------------------------------------------
    # ELECTRONS section
    # ----------------------------------------------------------

    fh.write(' &ELECTRONS\n')
    fh.write(' /\n')

    # ----------------------------------------------------------
    # | IONS section
    # ----------------------------------------------------------

    if p['calc'] in ['vc-relax', 'vc-md', 'md', 'relax']:
        fh.write(' &IONS\n')
        pwin_k=['ion_dynamics','ion_positions', 'phase_space', 'pot_extrapolation']
        write_section_params(fh, pwin_k, p)
        
        fh.write(' /\n')


    # ----------------------------------------------------------
    # | CELL section
    # ----------------------------------------------------------

    if p['calc'] in ['vc-relax', 'vc-md']:
        fh.write('&CELL\n')
        pwin_k=['cell_dynamics','press', 'cell_dofree']
        write_section_params(fh, pwin_k, p)
    
        fh.write('/\n')


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
    # Write out only if ibrav==0 - no symmetry used
    if not p['use_symmetry'] or p['ibrav']==0:
        fh.write('CELL_PARAMETERS angstrom\n')
        
        for v in cr.get_cell():
            fh.write('    %f %f %f\n' % tuple(v))
        
    # ----------------------------------------------------------
    # Card: K_POINTS
    # ----------------------------------------------------------

    if p['kpt_type'] in ['automatic','edos'] :
        fh.write('K_POINTS %s\n' % 'automatic')
    else :
        fh.write('K_POINTS %s\n' % p['kpt_type'])
        
    if p['kpt_type'] is 'automatic':
        fh.write('   %d %d %d   %d %d %d\n' % (tuple(p['kpts'])+tuple(p['kpt_shift'])))
    elif p['kpt_type'] is 'edos':
        fh.write('   %d %d %d   %d %d %d\n' % (tuple(p['nkdos'])+tuple([0,0,0])))
    elif not p['kpt_type'] is 'gamma':
        write_q_path(fh,p['qpath'],p['points'])
    fh.close()
    
    
def write_bands_in(d,a,p):
    p=p.copy()
    p.update({'iofile':'bands',
                'calc':'bands',
                'kpt_type':'tpiba_b' })
    write_pw_in(d,a,p)
    p.update({'kpt_type':'automatic'})
    
    
def write_edos_in(d,a,p):

    edos_pp_in='''&DOS
    prefix = '%(prefix)s',
    outdir = '%(outdir)s',
    fildos = '%(prefix)s.edos'
    /
    '''
    
    p=p.copy()
    if not ('occupations' in p.keys()) :
        p.update({'occupations':'tetrahedra'})
    
    p.update({'iofile':'edos',
                'calc':'nscf',
                'kpt_type':'edos'})
    write_pw_in(d,a,p)
    
    fh=open(os.path.join(d,'edos_pp.in'),'w')
    fh.write(edos_pp_in % p)
    fh.close()
    



def write_ph_in(d,a,p):
    # PH input file
    ph_in='''
    %(calc)s
    &INPUTPH
        prefix = '%(prefix)s',
        outdir = './%(outdir)s/',
        ldisp = .true.,
        fildrho = 'fildrho',
        fildvscf = 'fildvscf',
        nq1 = %(nq1)d , nq2 = %(nq2)d , nq3 = %(nq3)d
     /
    '''
    p=p.copy()
    qpts=p['qpts']
    p.update({'nq1':qpts[0],'nq2':qpts[1],'nq3':qpts[2]})
    fh=open(os.path.join(d,'ph.in'),'w')
    fh.write( ph_in % p )
    fh.close()


def write_q2r_in(d,a,p):
    q2r_in='''
    &INPUT
        fildyn='matdyn',
        zasr='%(asr)s',
        flfrc='%(prefix)s.fc'
    /
    '''

    p=p.copy()
    qpts=p['qpts']
    p.update({'nq1':qpts[0],'nq2':qpts[1],'nq3':qpts[2]})
    fh=open(os.path.join(d,'q2r.in'),'w')
    fh.write( q2r_in % p )
    fh.close()


def write_q_path(fh,path,points):
    path=array(path)
    fh.write( ' %d\n' % (len(path)))
    t=[]
    dt=[]
    s=0
    for x in [0]+map(norm,path[1:]-path[:-1]):
        s+=x
        t.append(s)
        dt.append(x)
    s=sum(dt)
    step=s/points
    for dist,q in zip(dt[1:]+[0],path):
        #print tuple(list(q)+[p['points']])
        fh.write(' %f %f %f   %d\n' % tuple(list(q)+[round(dist/step)]))
    

def write_matdyn_in(d,a,p):
    matdyn_in='''
    &input
        asr='%(asr)s',  
        flfrc='%(prefix)s.fc', 
        flfrq='%(prefix)s.freq',
        q_in_band_form=.true.
    /
    '''
    p=p.copy()
    qpts=p['qpts']
    p.update({'nq1':qpts[0],'nq2':qpts[1],'nq3':qpts[2]})
    fh=open(os.path.join(d,'matdyn.in'),'w')
    fh.write( matdyn_in % p )
    write_q_path(fh,p['qpath'],p['points'])
    fh.close()


def write_phdos_in(d,a,p):
    phdos_in='''
    &input
        dos=.true.
        asr='%(asr)s',  
        flfrc='%(prefix)s.fc', 
        fldos='%(prefix)s.dos',
        nk1=%(nk1)d,
        nk2=%(nk2)d,
        nk3=%(nk3)d,
        ndos=%(ndos)d
     /
    '''
    p=p.copy()
    qpts=p['qpts']
    nkdos=p['nkdos']
    p.update({'nq1':qpts[0],'nq2':qpts[1],'nq3':qpts[2]})
    p.update({'nk1':nkdos[0],'nk2':nkdos[1],'nk3':nkdos[2]})
    fh=open(os.path.join(d,'phdos.in'),'w')
    fh.write( phdos_in % p )
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
    
    
    

