import os, universe, utils, md, protein

def gen_constantpH(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE, cal=False, lambdaInit=0.5):
    # HARDCODED STUFF FOR PROTONATABLE RESIDUES ################################
    GLU_pKa   = 4.25
    GLU_atoms = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'] # atoms part of model
    GLU_qqA   = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    GLU_qqB   = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge
    
    ASP_pKa   = 3.65
    ASP_atoms = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'] # atoms part of model
    ASP_qqA   = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    ASP_qqB   = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge

    # HEAD #####################################################################

    # Skip this entire step if ph_constantpH is false.
    if (not universe.get('ph_constantpH')):
        utils.update("gen_constantpH", "ph_constantpH is False --> skipping...")
        return

    # If we use a charge-leveling scheme, we need the charge states:
    if (universe.get('ph_QQleveling') in [1, 2]):
        # BUF_qqA = [-0.0656, 0.5328, 0.5328] # previously hardcoded for water buffer
        # BUF_qqB = [-0.8476, 0.4238, 0.4238] # previously hardcoded for water buffer
        BUF_qqA = universe.get('ph_bufqqA')
        BUF_qqB = universe.get('ph_bufqqB')

    # If we use "charge-coupling" (1) scheme, extend charge states of GLU, ASP:
    if (universe.get('ph_QQleveling') == 1):
        GLU_qqA += BUF_qqB
        GLU_qqB += BUF_qqA
        ASP_qqA += BUF_qqB
        ASP_qqB += BUF_qqA

    # If we use "charge-constraining" (2) scheme, we need ph_BUF_dvdl:
    if (universe.get('ph_QQleveling') == 2):
        BUF_dvdl = universe.get('ph_BUF_dvdl')

    # Load other dV/dl coefficients:
    GLU_dvdl = universe.get('ph_GLU_dvdl')
    ASP_dvdl = universe.get('ph_ASP_dvdl')

    # Check whether MD.mdp exists.
    if (not os.path.isfile("MD.mdp")):
        utils.update("gen_constantpH", "MD.mdp does not exist, creating...")
        md.gen_mdp('MD', universe.get('d_nsteps'), universe.get('d_nstxout'))

    # Check whether index.ndx exists.
    if (not os.path.isfile("index.ndx")):
        utils.update("gen_constantpH", "index.ndx does not exist, creating...")
        utils.generate_index()

    ############################################################################
    
    file = open('MD.mdp', 'a')

    # Formatting function.
    def addParam(name, value, comment = "NUL"):
        if (comment == "NUL"):
            file.write("{:54s} = {:13s}\n".format(name, str(value)))
        else:            
            file.write("{:54s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    file.write("\n; CONSTANT PH\n")

    # PART 1 - WRITE GENERAL PARAMETERS ########################################
    
    # Update user.
    utils.update("gen_constantpH", "Writing general parameters:")
    utils.update("gen_constantpH", "ph_pH={}, ph_lambdaM={}, ph_nstout={}, ph_barrierE={}...".format(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE))

    addParam('lambda-dynamics', 'yes')
    addParam('lambda-dynamics-simulation-ph', ph_pH)
    addParam('lambda-dynamics-lambda-particle-mass', ph_lambdaM)
    addParam('lambda-dynamics-update-nst', ph_nstout)
    addParam('lambda-dynamics-tau', 2.0) # hardcoded

    if cal:
        addParam('lambda-dynamics-calibration', 'yes')

    if (universe.get('ph_QQleveling') == 2):
        addParam('lambda-dynamics-charge-constraints', 'yes')

    # Compile a list of acidic residues and their ResIDs.
    acidicResidueNameList = []; acidicResidueNumberList = []
    acidicResidueTypeList = []
    
    for residue in universe.get('d_residues'):
        if (residue.d_resname == 'GLU'):
            acidicResidueNameList.append('GLU')
            acidicResidueNumberList.append(residue.d_resid)
        
        if (residue.d_resname == 'ASP'):
            acidicResidueNameList.append('ASP')
            acidicResidueNumberList.append(residue.d_resid)

    if ('GLU' in acidicResidueNameList):
        acidicResidueTypeList.append('GLU')
    
    if ('ASP' in acidicResidueNameList):
        acidicResidueTypeList.append('ASP')

    # If we use "charge-restraining" we also have the BUF residue-type:
    if (universe.get('ph_QQleveling') == 2):  
        acidicResidueTypeList.append('BUF')

    addParam('lambda-dynamics-number-lambda-residues', len(acidicResidueTypeList))
    
    # If we use "charge-restraining" we also have the BUF lambda-group:
    if (universe.get('ph_QQleveling') == 2):
        addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList) + 1)
    else:
        addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList))

    file.write('\n')

    # print(acidicResidueNameList)   # debug
    # print(acidicResidueNumberList) # debug
    # print(acidicResidueTypeList)   # debug

    # PART 2 - WRITE RESIDUE-TYPE SPECIFIC STUFF ###############################

    utils.update("gen_constantpH", "Writing residue-type specific stuff...")

    def writeBlock(number, name, dvdl, pKa, ph_barrierE, qqA, qqB):

        def to_string(Input):
            string = ""
            for element in Input:
                string += "{:.3f}".format(element)
                string += ' '
            return string

        addParam('lambda-dynamics-residue%s-name'              % (number), name)
        addParam('lambda-dynamics-residue%s-dvdl-coefficients' % (number), to_string(dvdl))
        addParam('lambda-dynamics-residue%s-reference-pka'     % (number), pKa)
        addParam('lambda-dynamics-residue%s-barrier'           % (number), ph_barrierE)
        addParam('lambda-dynamics-residue%s-charges-state-A'   % (number), to_string(qqA))
        addParam('lambda-dynamics-residue%s-charges-state-B'   % (number), to_string(qqB))
        
        file.write('\n')

    for idx in range(0, len(acidicResidueTypeList)):
        if (acidicResidueTypeList[idx] == 'GLU'):
            writeBlock(idx + 1, 'GLU', GLU_dvdl, GLU_pKa, ph_barrierE, GLU_qqA, GLU_qqB)

        if (acidicResidueTypeList[idx] == 'ASP'):
            writeBlock(idx + 1, 'ASP', ASP_dvdl, ASP_pKa, ph_barrierE, ASP_qqA, ASP_qqB)

        if (acidicResidueTypeList[idx] == 'BUF'):
            # If number of protonatable residues != number of buffer molecules,
            # we need to increase buffer charge in state A by the ratio:
            nLams   = protein.countRes('ASP') + protein.countRes('GLU')
            nBufs   = universe.get('ph_bufnmol')
            BUF_qqA = [(nLams/float(nBufs)) * i for i in BUF_qqA]

            writeBlock(idx + 1, 'BUF', BUF_dvdl, 0, 0, BUF_qqA, BUF_qqB)

    # PART 3 - WRITE INDIVIDUAL RESIDUE/LAMBDA-GROUP STUF ######################

    utils.update("gen_constantpH", "Writing individual lambda groups...")

    def writeResBlock(number, name, indexLambda, indexName):
        addParam('lambda-dynamics-atom-set%s-name'                  % (number), name)
        addParam('lambda-dynamics-atom-set%s-lambda-residues-index' % (number), indexLambda)
        addParam('lambda-dynamics-atom-set%s-index-group-name'      % (number), indexName)
        addParam('lambda-dynamics-atom-set%s-initial-lambda'        % (number), lambdaInit)
        
        if (universe.get('ph_QQleveling') == 2):
            addParam('lambda-dynamics-atom-set%s-charge-restraint-group-index' % (number), 1)

        if (name == 'BUF'):
            addParam('lambda-dynamics-atom-set%s-buffer-residue' % (number), 'yes')
            addParam('lambda-dynamics-atom-set%s-buffer-residue-multiplier' % (number), universe.get('ph_bufnmol'))

        file.write('\n')

    for idx in range(0, len(acidicResidueNameList)):
        writeResBlock(
                        idx + 1, 
                        acidicResidueNameList[idx],
                        acidicResidueTypeList.index(acidicResidueNameList[idx]) + 1,
                        'LAMBDA%s' % (idx + 1)
                        )

    if (universe.get('ph_QQleveling') == 2):
        writeResBlock(
                        len(acidicResidueNameList) + 1,
                        'BUF',
                        acidicResidueTypeList.index('BUF') + 1,
                        'LAMBDA%s' % (len(acidicResidueNameList) + 1)
                        )

    file.close() # MD.mdp

    # PART 4 - APPEND THE LAMBDA INDEX GROUPS TO INDEX.NDX #####################

    utils.update("gen_constantpH", "Writing lambda index groups to index.ndx...")

    # If we use a charge-leveling scheme, we need a list of atomIndices of the BUFs:
    if (universe.get('ph_QQleveling') in [1, 2]):
        bufferAtomIndexList = []

        atomCount = 1
        for residue in universe.get('d_residues'):
            for atom in residue.d_atoms:
                if (residue.d_resname == 'BUF'):
                    bufferAtomIndexList.append(atomCount)
            
                atomCount += 1

    file = open('index.ndx', 'a') # Append to existing index.ndx

    # Function for adding an atomIndexList to index.ndx
    def writeTheGroup(number, atomIndexList):
        file.write('\n[ LAMBDA{} ]\n'.format(number))
        for index in atomIndexList:
            file.write('{} '.format(index))
        file.write('\n')

    grpNum = 1 # Keeps track of the group (the LAMBDA%s).

    atomCount = 1
    for residue in universe.get('d_residues'):  # loop through all residues

        atomIndexList = []                  # clear atomIndexList

        for atom in residue.d_atoms:        # for each residue, loop through the atoms

            if (residue.d_resname == 'GLU' and atom in GLU_atoms):
                atomIndexList.append(atomCount)

            elif (residue.d_resname == 'ASP' and atom in ASP_atoms):
                atomIndexList.append(atomCount)

            atomCount += 1                  # increment atomCount

        if (len(atomIndexList) > 0):
            # If we use "charge-coupling" (1), assign the atomIndices of one BUF
            # to one protonatable lambda residue (use clever list slicing):
            if (universe.get('ph_QQleveling') == 1):
                start = (grpNum - 1) * len(BUF_qqA)
                stop  = start + len(BUF_qqA)
                atomIndexList += bufferAtomIndexList[start:stop]

            writeTheGroup(grpNum, atomIndexList)
            grpNum += 1

    # If we use "charge-restraining" (2), add everything in bufferAtomIndexList
    # to the last lambda index group:
    if (universe.get('ph_QQleveling') == 2):
        writeTheGroup(grpNum, bufferAtomIndexList)

    file.close() # index.ndx

    # Put relevant pH variables in universe
    universe.add('ph_pH', ph_pH)
    universe.add('ph_lambdaM', ph_lambdaM)
    universe.add('ph_nstout', ph_nstout)
    universe.add('ph_barrierE', ph_barrierE)
