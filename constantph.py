import os, universe, utils, md

def gen_constantpH(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE, cal=False, lambdaInit=0.5):
    # Skip this entire step if ph_constantpH is false.
    if (not universe.get('ph_constantpH')):
        utils.update("gen_constantpH", "ph_constantpH is False --> skipping...")
        return

    # Check whether MD.mdp exists and create if necessary.
    if (not os.path.isfile("MD.mdp")):
        utils.update("gen_constantpH", "MD.mdp does not exist, creating...")
        md.gen_mdp('MD', universe.get('d_nsteps'), universe.get('d_nstxout'))

    # Check whether index.ndx exists, and create if necessary.
    if (not os.path.isfile("index.ndx")):
        utils.update("gen_constantpH", "index.ndx does not exist, creating...")
        utils.generate_index()

    ############################################################################
    file = open('MD.mdp', 'a')

    def addParam(name, value): # Formatting function.
            file.write("{:54s} = {:13s}\n".format(name, str(value)))

    # Update user.
    utils.update("gen_constantpH", "Writing general parameters:")
    utils.update("gen_constantpH", "ph_pH={}, ph_lambdaM={}, ph_nstout={}, ph_barrierE={}...".format(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE))

    # PART 1 - WRITE GENERAL PARAMETERS ########################################
    file.write("\n; CONSTANT PH\n")

    addParam('lambda-dynamics', 'yes')
    addParam('lambda-dynamics-simulation-ph', ph_pH)
    addParam('lambda-dynamics-lambda-particle-mass', ph_lambdaM)
    addParam('lambda-dynamics-update-nst', ph_nstout)
    addParam('lambda-dynamics-tau', 2.0) # hardcoded

    if cal: # If we are in calibration mode:
        addParam('lambda-dynamics-calibration', 'yes')

    # If we use "charge-constraining" (2) scheme:
    if (universe.get('ph_QQleveling') == 2):
        addParam('lambda-dynamics-charge-constraints', 'yes')

    # Gather a list of the different lambda residue-type names that were 
    # specified IN THE UNIVERSE by the user.
    lambdaTypeNamesSpecifed = []
    for obj in universe.get('ph_lambdaTypes'):
        lambdaTypeNamesSpecifed.append(obj.d_resname)

    # Gather a list of the names of the lambda residues in the protein.
    lambdaResidueNameList = []
    for residue in universe.get('d_residues'):
        if (residue.d_resname in lambdaTypeNamesSpecifed):
            lambdaResidueNameList.append(residue.d_resname)

    # Gather a list of the names of the lambda residue-types IN THE PROTEIN,
    # And make sure the order is the same as in ph_lambdaTypes = LambdaTypeNamesSpecified.
    lambdaResidueTypeList = []
    for obj in universe.get('ph_lambdaTypes'):
        if (obj.d_resname in set(lambdaResidueNameList)):
            lambdaResidueTypeList.append(obj.d_resname)

    # If we use the charge leveling scheme "charge-restraining" (2) we also have
    # the BUF residue-type as well as one extra lambda group containing all the BUFs.
    if (universe.get('ph_QQleveling') == 2):
        addParam('lambda-dynamics-number-lambda-residues', len(lambdaResidueTypeList) + 1)
        addParam('lambda-dynamics-number-atom-collections', len(lambdaResidueNameList) + 1)
    else:
        addParam('lambda-dynamics-number-lambda-residues', len(lambdaResidueTypeList))
        addParam('lambda-dynamics-number-atom-collections', len(lambdaResidueNameList))

    file.write('\n')

    # PART 2 - WRITE RESIDUE-TYPE SPECIFIC STUFF ###############################
    def writeBlock(number, name, dvdl, pKa, ph_barrierE, qqA, qqB):
        def to_string(Input):
            string = ""
            for element in Input:
                string += "{:.3f} ".format(element)
            return string

        addParam('lambda-dynamics-residue%s-name'              % (number), name)
        addParam('lambda-dynamics-residue%s-dvdl-coefficients' % (number), to_string(dvdl))
        addParam('lambda-dynamics-residue%s-reference-pka'     % (number), pKa)
        addParam('lambda-dynamics-residue%s-barrier'           % (number), ph_barrierE)
        addParam('lambda-dynamics-residue%s-charges-state-A'   % (number), to_string(qqA))
        addParam('lambda-dynamics-residue%s-charges-state-B'   % (number), to_string(qqB))

        file.write('\n')

    # If we use a charge-leveling scheme, we need the buffer charge states.
    # These are added to universe when you run protein.add_buffer()
    if (universe.get('ph_QQleveling') in [1, 2]):
        BUF_qqA = universe.get('ph_bufqqA') # BUF_qqA = [-0.0656, 0.5328, 0.5328] previously hardcoded for water buffer
        BUF_qqB = universe.get('ph_bufqqB') # BUF_qqB = [-0.8476, 0.4238, 0.4238] previously hardcoded for water buffer

    idx = 1
    for obj in universe.get('ph_lambdaTypes'):
        # This if-statement prevents writing a block when there are no residues of this type.
        if (obj.d_resname in lambdaResidueTypeList):
            # If we use "charge-coupling" (1) scheme, extend the charge states:
            if (universe.get('ph_QQleveling') == 1):
                writeBlock(idx, obj.d_resname, obj.d_dvdl[::-1], obj.d_pKa, ph_barrierE, obj.d_qqA + BUF_qqB, obj.d_qqB + BUF_qqA)
            else:
                writeBlock(idx, obj.d_resname, obj.d_dvdl[::-1], obj.d_pKa, ph_barrierE, obj.d_qqA, obj.d_qqB)
            idx += 1

    if (universe.get('ph_QQleveling') == 2):
        # f we use "charge-constraining" (2) scheme, we additionaly need ph_BUF_dvdl.
        writeBlock(idx, 'BUF', universe.get('ph_BUF_dvdl')[::-1], 0, 0, BUF_qqA, BUF_qqB)

    # PART 3 - WRITE INDIVIDUAL RESIDUE/LAMBDA-GROUP STUF ######################
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

    utils.update("gen_constantpH", "Writing individual lambda groups...")

    idx = 1
    for name in lambdaResidueNameList:
        writeResBlock(idx, name, lambdaResidueTypeList.index(name) + 1, 'LAMBDA{}'.format(idx))
        idx += 1

    if (universe.get('ph_QQleveling') == 2):
        writeResBlock(idx, 'BUF', len(lambdaResidueTypeList) + 1, 'LAMBDA{}'.format(len(lambdaResidueNameList) + 1))

    file.close() # MD.mdp

    # PART 4 - APPEND THE LAMBDA INDEX GROUPS TO INDEX.NDX #####################
    utils.update("gen_constantpH", "Writing lambda index groups to index.ndx...")

    # If we use a charge-leveling scheme, we need a list of atomIndices of the BUFs:
    if (universe.get('ph_QQleveling') in [1, 2]):
        bufferAtomIndexList = []

        count = 1
        for residue in universe.get('d_residues'):
            for atom in residue.d_atoms:
                if (residue.d_resname == 'BUF'):
                    bufferAtomIndexList.append(count)

                count += 1

    file = open('index.ndx', 'a') # Append to existing index.ndx

    def writeTheGroup(number, atomIndexList):
        file.write('\n[ LAMBDA{} ]\n'.format(number))
        for index in atomIndexList:
            file.write('{} '.format(index))
        file.write('\n')

    ph_lambdaTypes = universe.get('ph_lambdaTypes')

    atomCount = 1; groupNumber = 1
    for residue in universe.get('d_residues'):
        if residue.d_resname in lambdaResidueTypeList:

            atomIndexList = []
            obj = ph_lambdaTypes[lambdaTypeNamesSpecifed.index(residue.d_resname)]

            for atom in residue.d_atoms:
                if atom in obj.d_atoms:
                    atomIndexList.append(atomCount)

                atomCount += 1

            # If we use "charge-coupling" (1), assign the atomIndices of one BUF
            # to one protonatable lambda residue (use clever list slicing):            
            if (universe.get('ph_QQleveling') == 1):
                start = (groupNumber - 1) * len(BUF_qqA)
                stop  = start + len(BUF_qqA)
                atomIndexList += bufferAtomIndexList[start:stop]                

            writeTheGroup(groupNumber, atomIndexList)
            groupNumber += 1

        else:
            for atom in residue.d_atoms:
                atomCount += 1

    # If we use "charge-restraining" (2), add everything in bufferAtomIndexList
    # to the last lambda index group:
    if (universe.get('ph_QQleveling') == 2):
        writeTheGroup(groupNumber, bufferAtomIndexList)

    file.close() # index.ndx

    # Put relevant pH variables in universe
    universe.add('ph_pH', ph_pH)
    universe.add('ph_lambdaM', ph_lambdaM)
    universe.add('ph_nstout', ph_nstout)
    universe.add('ph_barrierE', ph_barrierE)
