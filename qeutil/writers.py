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
    A = %(A)f,
    ecutwfc = %(ecut)f,
    ibrav = 2,
    nat = 3,
    ntyp = 2,
/
&ELECTRONS
/
ATOMIC_SPECIES
 Th 232.038100 Th_%(XC)s_nc.ncpp
 O 15.999400 O_%(XC)s_nc.ncpp
ATOMIC_POSITIONS
 Th 0.000000  0.000000  0.000000
 O 0.250000  0.250000  0.250000
 O -0.250000  -0.250000  -0.250000
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
    nk1=15,
    nk2=15,
    nk3=15,
    ndos=%(ndos)d
 /
'''

def make_calc_dir(bdir,params):
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
        bdir - base directory
        params - calculation parameters dictionary
    
    OUTPUT
    ------
        cdir - name of the created calculation directory
    
    '''
    
    cdir=params['prefix']+'.XXXX'
    cdir=!cd $bdir ; mktemp -d $cdir
    cdir=cdir[0]
    bcdir='calc/'+cdir
    #print cdir, bcdir
    
    open(bdir + cdir +'/pw.in','w').write(pw_in % params)
    open(bdir + cdir +'/ph.in','w').write(ph_in % params)
    open(bdir + cdir +'/q2r.in','w').write(q2r_in % params)
    f=open(bdir + cdir+'/matdyn.in','w')
    f.write(matdyn_in % params)
    f.write('%d\n'%qp.shape[0])
    for v in qp:
        f.write("%f %f %f %d\n" % (v[0],v[1],v[2],params['points']))
    f.close()
    open(bdir+cdir+'/phdos.in','w').write(phdos_in % params)
    return bcdir,cdir

def make_D3_files(bdir,host,clc,cdat):
    cdir = bdir + host + clc
    open(cdir +'/phD3G.in','w').write(ph3G_in % cdat)
    open(cdir +'/D3G.in','w').write(d3G_in % cdat)
    open(cdir +('/phD3_%(sufix)s.in' % cdat),'w').write(ph3any_in % cdat)
    open(cdir +('/D3_%(sufix)s.in' % cdat),'w').write(d3any_in % cdat)
    
    
    

