def read_eigenvect(fn):
    '''Read the normalized eigenvectors/sqrt(mass)=displacement vectors 
    from the ph.x output file. This is different, but similar to the matdyn.modes file.
    Returns a list of lists:
    [ 
        [q, [[omega, ev], [omega, ev], ...]],
        [q, [[omega, ev], [omega, ev], ...]]
        ...
    ]
    where ev is an #atomsx3 complex array of eigenvectors for a given mode with frequency omega.
    The indexes in ev are [atomic,cartesian].
    '''
    
    evl=[]
    start=False
    for ln in open(fn).readlines():
        ln=ln.split()
        if len(ln)==0:
            # skip empty lines
            continue
        if ln[0].find('iagonalizing')==1 :
            # We are in the eigenvectors block
            start=True
        if start and ln[0]=='q' :
            # read the q-vector
            if len(ln)==5 :
                # we are reading the matdyn.modes file
                q=array(ln[2:],dtype=float)
            else :
                # we are reading the fildyn file
                q=array(ln[3:6],dtype=float)
            ql=[]
            evl.append(ql)
            ql.append(q)
            modes=[]
            ql.append(modes)
        if start and ln[0].find('omega')==0 :
            # read the frequency of the mode
            o=float(ln[-5])
            mode=[]
            modes.append(mode)
            mode.append(o)
        if start and ln[0].find('(')==0 :
            # read the atomic displacements line by line converting to complex on the flight
            mode.append(array(ln[1:7:2],dtype=float)+1j*array(ln[2:7:2],dtype=float))
            
    return [[q[0], [[m[0],array(m[1:])] for m in q[1]]] for q in evl]    
    
    
def read_D2(fn='matdyn', wd='.'):
    '''Read the dynamical matrix sampling over the B.Z. from the matdyn files.
    Input
    -----
        fn - base file name of the dynamical matrix, *without* the zero at the end. 
        wd - working directory
    
    Output
    ------
        Dynamical matrix of the crystal sampled over q-points in the B.Z. as a list.
        D2[na, nb, i, j] is a nat x nat x 3 x 3 dynamical matrix at the q point in B.Z.
        nat - number of atoms:
        
            [[q, D2[na, nb, i, j]],
             [q, D2[na, nb, i, j]],
             ...
            ]
    '''
    
    # Read the description file
    
    qvect=loadtxt(os.path.join(wd,fn+'0'),skiprows=2)
    qlatt=array(open(os.path.join(wd,fn+'0')).readline().split(),dtype=int)
    #print "Q-lattice: ", qlatt, ',', len(qvect), " irreducible q-vectors:"
    #print qvect
    
    indm=False
    D2=[]
    for k in range(1,len(qvect)+1):
        # Collect the matrices for each q-vect set
        inhd=True
        lnum=0
        nat=0
        for ln in open(os.path.join(wd,fn+("%d" % (k,)))):
            fld=ln.split()
            lnum+=1
            if inhd and lnum==3 :
                nat=int(fld[1])
            if ln.find('Dynamical  Matrix in cartesian axes') > -1 :
                indm=True
                inhd=False
            if ln.find('Dielectric Tensor:') > -1 :
                indm=False
            if ln.find('Effective Charges') > -1 :
                indm=False
            if ln.find('Diagonalizing the dynamical matrix') > -1 :
                indm=False
            if indm and ln.find('q = (')>-1 :
                q=[]
                q.append(array(fld[-4:-1],dtype=float))
                D2.append(q)
            if indm and len(fld)==2 :
                s_sp=array(fld,dtype=int)-1
            if indm and len(fld)==6 :
                dxyz=array(fld[::2],dtype=float)+1j*array(fld[1::2],dtype=float)
                q.append(dxyz)
                #print dxyz
    # Re-process all D2 matrices into nat x nat x 3 x 3 matrices and put
    return [[ql[0],array(ql[1:]).reshape(nat,nat,3,3)] for ql in D2]
    
def read_D3(fn='d3dyn',wd="."):
    '''Reads d3dyn file from QE d3.x program containing D3 matrix of PRB 87 (2013) 214303.
    Returns q-vectors and re-ordered complex D3 matrix.
    The returned d3 matrix has following indexing (left to right):
    q, index of q-vector in qvect array
    i, cart. index for atom a
    j, cart. index for atom b
    k, cart. index for atom c
    na, index of atom a
    nb, index of atom b
    nc, index of atom c
    and shape (#q, 3, 3, 3, #atoms, #atoms, #atoms)
    '''
    fd3dyn=open(wd+'/'+fn)
    fd3dyn=[ln.split() for ln in fd3dyn]
    nat=int(fd3dyn[2][1])
    qvect=[]
    head=0
    for ln in fd3dyn:
        if not qvect:
            head+=1
        if len(ln)>2 and ln[0]=='q' and ln[1]=='=':
            qvect.append(array([float(x) for x in ln[-4:-1]]))
    modes=[]
    for ln in fd3dyn[head+1:]:
        if len(ln)==3:
            modes.append([float(x) for x in ln])
    qvect=array(qvect)
    modes=array(modes)
    modes=modes.reshape((len(qvect),nat,3,nat,nat,3,3,2))
    # collect real and imaginary parts into complex numbers
    modes=modes[...,0]+1j*modes[...,1]
    # Account for ordering of d3.x output.
    # The d3 matrix has following indexing:
    # d3(
    #    q, index of q-vector in qvect array
    #    i, cart. index for atom a
    #    j, cart. index for atom b
    #    k, cart. index for atom c
    #    na, index of atom a
    #    nb, index of atom b
    #    nc, index of atom c
    # if you just read in the numbers the ordering will be 
    # (slowest to fastest index):
    # q, na, i, nb, nc, j, k, (re/im) 
    d3=modes.transpose(0,2,5,6,1,3,4)
    return qvect, d3
    
    
    
# -------------------------------------------------------------------------------------------------------------------------------
#
#  read infile
#
# -------------------------------------------------------------------------------------------------------------------------------    


# ......................................................................................

def build_atoms(infilename, myprefix="dry_pwscf", verbose = False):
    """
    this function builds an Atoms object in ASE starting
    from a .in file, using the data collected by
    read_in_file()
    blocks must be the dictionary returned by read_in_file()

    it returns an Atoms object
    """

    new_infile = infilename+'.dry'
    
    blocks = read_in_file(infilename)
    # changes the prefix in myprefix
    if blocks.has_key('prefix'):
        _command_string = "sed 's|prefix[ ]*=[ ]*[a-z/0-9]*|prefix = "+myprefix+",|' < "+infilename+" > "+new_infile
    else:
        _command_string = "sed 's|\&control|\&control\\n    prefix = "+myprefix+",|' < "+infilename+" > "+new_infile

    print _command_string
    os.system(_command_string)
        
    
    if blocks['system']['ibrav'] == 0:
        if verbose:
            print "\t[write_atoms] ibrav = 0 found\n"
        
        if blocks['system'].has_key('alat'):
            alat = blocks['system']['alat']
        elif blocks['system'].has_key('celldm(1)'):
            alat = blocks['system']['celldm(1)']
        elif blocks['system'].has_key('a'):
            alat = blocks['system']['a']
        else:
            raise ValueError, "something got wrong: it seems that ibrav=0 but neither 'alat' nor 'celldm(1)' nor 'a' are present"

        # creates the cell as an array
        cell = np.array(blocks['CELL_PARAMETERS']['cell'])
        if verbose:
            print "\t[write_atoms] found cell:\n"
            for R in cell:
                print "\t\t %f %f %f" % tuple(R)
            
        # convert the cell in bohr units
        # (note that alat or celldm(1) are supposed to be already in bohr)
        if blocks['CELL_PARAMETERS']['attrib'] == 'alat' or blocks['CELL_PARAMETERS']['attrib'] == None:
            cell = cell * alat
        elif blocks['CELL_PARAMETERS']['attrib'] == 'angstrom':
            cell = cell * convert('a','bohr')
        if verbose:
            print "\t[write_atoms] cell rescaled is:\n"
            for R in cell:
                print "\t\t %f %f %f" % tuple(R)
            
            
    else:
        # it is needed to start qe in dry mode to obtains the
        # cell parameters, because ibrav >0
        # QE will do this for us

        if verbose:
            print "\t[write_atoms] ibrav > 0 found. Invoking QE\n"
        
        if 'ESPRESSO_ROOT' in os.environ:
                rootdir = os.environ['ESPRESSO_ROOT']
                bindir = rootdir + '/PW/src/'
        else:
            rootdir= '$HOME/espresso/'
            bindir = rootdir + '/PW/src/'

        execbin = "pw.x"
        if not os.access(bindir+execbin, os.F_OK):
            raise IOError, "binary %s does not exist" % (bindir+execbin)
        if not os.access(bindir+execbin, os.X_OK):
            raise IOError, "you do not have execution permission on the binary %s" % bindir+execbin

        # run a dry run: only the xml data file is written, with the information about the cell
        tempfilename = myprefix+".EXIT"
        fh_temp = open(tempfilename, "w")
        fh_temp.close()
        _commandstring = '%s < %s > %s' % ( bindir+execbin, new_infile, new_infile+'.out')
        exitcode = os.system(_commandstring)
        
        # read in cell information
        if verbose:
            print "\t[write_atoms] now read-in the xml output\n"
        root = ET.ElementTree(file=myprefix+'.save/data-file.xml').getroot()
        mydict = {}
        get_cell_xml(root, mydict)

        alat = mydict['alat']
        # the cell is already in bohr
        cell = mydict['lattice_vector']

        if verbose:
            print "\t[write_atoms] found cell:\n"
            for R in cell:
                print "\t\t %f %f %f" % tuple(R)

    # now actually creates the Atoms structure
        
    atoms = Atoms()
    for A in blocks['ATOMIC_POSITIONS']['atoms']:
        print "\t[write_atoms] adding atom: %s\n" % A['symbol']
        atoms.append(Atom(A['symbol'], tuple(A['pos'])))
    print "\t[write_atoms] attaching cell:\n"
    for R in cell:
        print "\t\t %f %f %f" % tuple(R)
    
    atoms.set_cell(cell, scale_atoms=True)

    return atoms


# ......................................................................................

def read_in_file(infilename):
    """
    this function read a quantum espresso  .in file
    and returns a dictionary for all the namelists
    and the cards found.
    it correctly deals lines with more than one
    comma-separated parameters

    it returns a dictionary, whos keys are the name
    of the namelists and cards.
    In turn, each value is a dictionary with all the
    values found.

    as for the namecards:
    (+) atomic_species is a dictionary whos keys are progeressive
        integers and whose values are dictionaries with the following
        keys:
        'symbol', 'attrib', 'mass', 'PP'
        example:
        {'1': {'symbol: 'H', 'attrib': 1, 'mass': 1.0008, 'PP': "pseudo_file_for_H"},
        '2': {'symbol: 'C', 'attrib': None, 'mass': 1, 'PP': "pseudo_file_for_C"},
        '3': {'symbol: 'H', 'attrib': 2, 'mass': 1.0008, 'PP': "different_pseudo_file_for_H"},
        ...}
        the key 'attrib' is set in case there is more than 1 entry for a specie, for instance
        to specify different PP:
        ATOMIC_SPECIES
         H1 1.0008 pseudo_file_for_H
         H2 1.0008 different_pseudo_file_for_H
         C  12 pseudo_file_for_C

    (+) atomic_positions is a dictionary whos keys are progressive
    integers and whose values are dictionaries with the following keys:

    example: 'symbol', 'pos', 'if_pos'
    {'1': {'symbol: 'H', 'pos': [1.0, 1.0, 1.0], 'if_pos': [0.1, 0.1, 0.1]},
     ...}
     the key 'if_pos' may have value None if it wasn't present

     (+) all the other namecards are dictionaries with a unique key, 'elements',
        that is a list of all the lines inside that card
     
    """

    if not os.access(infilename, os.F_OK):
        raise IOError, "the .in file %s is not present\n" % infilename

    fh = open(infilename)
    alltext = fh.read()
    alllines = alltext.split('\n')
    
    card_labels = ['ATOMIC_SPECIES',
                   'ATOMIC_POSITIONS',
                   'CELL_PARAMETERS',
                   'K_POINTS',
                   'CONSTRAINTS'
                   'OCCUPATIONS']

    blocks = {}
    blocks['sparse']=[]
    isnamelist = False
    iscard = False

    for line in alllines:

        # determines whether we are inside a namelist or a card
        if line.strip().startswith('&'):
            isnamelist = True
            blockname = line.strip()[1:]
            blocks[blockname] = {}

        elif line.strip().startswith('/'):
            isnamelist = False

        elif line.strip():
            if line.split()[0].strip() in card_labels:
                iscard = True
                
                blockname = line.strip().split()[0]
                blocks[blockname] = {}

                try:
                    attrib = line[line.find(line.strip().split()[1]):].strip()
                    if attrib.startswith('{'):
                        attrib = attrib[1:]
                    if attrib.endswith('}'):
                        attrib = attrib[:-1]
                except:
                    attrib = None
                blocks[blockname]["attrib"] = attrib

            else:     # if in a namelist, isolate keywords and values

                if isnamelist:
                    tokens = line.split(',')
                    for t in tokens:
                        if t.strip():

                            key, svalue = t.strip().split('=')

                            key = key.strip()
                            svalue = svalue.strip()
                            if svalue.endswith(','):
                                value = svalue[:-1].strip()
                            try:
                                value=int(svalue)
                            except ValueError:
                                try:
                                    value=float(svalue)
                                except ValueError:
                                    if svalue.lower() in ['.true.', 'true' , 't']:
                                        value = True
                                    elif svalue.lower() in ['.false.', 'false' , 'f']:
                                        value = False
                                    else:
                                        value = str(svalue)

                            blocks[blockname][key] = value

                elif iscard:

                    if blockname == "ATOMIC_SPECIES":                        # -- SPECIES
                        tokens = line.split()
                        if(blocks[blockname].has_key('count')):
                            blocks[blockname]['count'] = blocks[blockname]['count']+1
                        else:
                            blocks[blockname]['count'] = 1
                            
                        symbol = tokens[0].strip()
                        
                        if symbol.find('-')> 0 or symbol.find('_')>0:
                            ssymbol = symbol.split('-')
                            if len(ssymbol) == 1:
                                ssymbol = symbol.split('_')
                            symbol=ssymbol[0]
                            symbol_attrib=ssymbol[1]
                        elif not symbol.isalpha():
                            digits="0123456789"
                            found=min([symbol.index(d) for d in digits if d in symbol])
                            symbol_attrib = int(symbol[found:])
                            symbol = symbol[:found]
                        else:
                            symbol_attrib = None
                            
                        blocks[blockname][str(blocks[blockname]['count'])] = {'symbol': symbol,
                                                                              'attrib': symbol_attrib,
                                                                               'mass' : float(tokens[1]),
                                                                               'PP' : tokens[2].strip()}

                    elif blockname == "ATOMIC_POSITIONS":                    # -- POSITIONS
                        tokens = line.split()
                        if(blocks[blockname].has_key('count')):
                            blocks[blockname]['count'] = blocks[blockname]['count']+1
                        else:
                            blocks[blockname]['count'] = 1
                            blocks[blockname]['atoms'] = []

                        idx = blocks[blockname]['count']
                        if len(tokens) > 4:
                            if_pos = [float(eval(s)) for s in tokens[4:]]
                        else:
                            if_pos = None

                        blocks[blockname]['atoms'].append({'symbol': tokens[0].strip(),
                                                          'pos': [float(eval(s)) for s in tokens[1:4]],
                                                          'if_pos': if_pos})

                    elif blockname == "CELL_PARAMETERS":                     # -- CELL PARAMETERS

                        tokens = line.split()

                        if not blocks[blockname].has_key('cell'):
                            blocks[blockname]['cell'] = []

                        blocks[blockname]['cell'].append([float(s) for s in tokens[0:3]])
                        if(len(blocks[blockname]['cell'])) == 3:
                            iscard = False                                                

                    else:
                        if not blocks[blockname].has_key('elements'):
                            blocks[blockname]['elements'] = []
                        blocks[blockname]['elements'].append(line)

                else:
                    blocks['sparse'].append(line)


    fh.close()
    return blocks
