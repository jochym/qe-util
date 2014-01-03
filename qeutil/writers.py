from tempfile import mkdtemp
import os

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


def write_pw_in(d,a,p):

    fh=open(os.path.join(d,'pw.in'),'w')
    # ----------------------------------------------------------
    # | CONTROL section
    # ----------------------------------------------------------

    fh.write(' &CONTROL\n')
    fh.write("    calculation = '%s',\n" % p['calc'])

    pwin_k=['tstress','nstep','pseudo_dir','outdir']
    write_section_params(fh, pwin_k, p)
    fh.write(' /\n')
    
    # ----------------------------------------------------------
    # | SYSTEM section
    # ----------------------------------------------------------

    pwin_k=['ecutwfc','ibrav']
    fh.write(' &SYSTEM\n')
    write_section_params(fh, pwin_k, p)
    fh.write("    nat = %d\n" % (a.get_number_of_atoms()))
    fh.write("    ntyp = %d\n" % (len(set(a.get_atomic_numbers()))))
    fh.write(' /\n')

    # ----------------------------------------------------------
    # | ELECTRONS section
    # ----------------------------------------------------------

    fh.write(' &ELECTRONS\n')
    fh.write(' /\n')

    
    # ----------------------------------------------------------
    # | Card: ATOMIC_SPECIES
    # ----------------------------------------------------------
        
    fh.write('ATOMIC_SPECIES\n')
    
    xc=p['xc']
    pp_type=p['pp_type']
    pp_format=p['pp_format']
    for nm, mass in set(zip(a.get_chemical_symbols(),a.get_masses())):
        fh.write("    %s %g %s_%s_%s.%s \n" % (nm, mass, nm, xc, pp_type, pp_format))

    # ----------------------------------------------------------
    # | Card: ATOMIC_POSITIONS
    # ----------------------------------------------------------
    
    fh.write('ATOMIC_POSITIONS {crystal}\n')
    for n,v in zip(a.get_chemical_symbols(),a.get_scaled_positions()):
        fh.write("    %s %g %g %g\n" % (n, v[0], v[1], v[2]))

    # ----------------------------------------------------------
    # | Card: CELL_PARAMETERS
    # ----------------------------------------------------------
    
    fh.write('CELL_PARAMETERS {angstrom}\n')
    
    for v in a.get_cell():
        fh.write('    %f %f %f\n' % tuple(v))
    

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
    
    
    

