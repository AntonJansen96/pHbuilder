import os, universe, utils, protein

def add_mol(itpfname, comment, molname=None, molcount=None):
    # Get the contents of current topol.top.
    topList = []
    with open("topol.top") as file:
        for line in file.readlines():
            topList.append(line)
    
    # Add the .itp file (line saying: #include "blabla.itp")
    with open("topol.top", 'w') as file:
        try:
            for line in range(0, len(topList)):
                file.write(topList[line])
                
                if "[ system ]\n" == topList[line + 1]:
                    file.write("; {0}\n".format(comment))
                    file.write("#include \"{0}\"\n\n".format(itpfname))

        except IndexError:
            pass

    # if molcount not present, add it, otherwise do nothing.
        if molname != None and molcount != None and molname not in topList[-1]:
            file.write("{0}\t\t\t{1}\n".format(molname, molcount))

def generate(d_modelFF, d_modelWater, d_terministring=""):
    # Internal helper function.
    def rebuild_topol():
        # If we have only one chain, gromacs will put everything in topol.top.
        # If we have more than one chain, gromacs will do it for us.
        if (len(universe.get('d_chain')) <= 1):
            readingProtein = False
            
            file = open("topol_Protein_chain_A.itp", 'w')

            for line in open("topol.top").readlines():
                if (not readingProtein and line == "[ moleculetype ]\n"):
                    readingProtein = True
            
                if (readingProtein and line == "; Include water topology\n"):
                    readingProtein = False

                if (readingProtein):
                    file.write(line)
    
            file.close()
    
        with open('topol.top', 'w') as file:
            file.write("; Include forcefield parameters\n")
            file.write("#include \"{0}.ff/forcefield.itp\"\n\n".format(universe.get('d_modelFF')))

            file.write("; Include protein topology\n")
            for letter in universe.get('d_chain'):
                file.write("#include \"topol_Protein_chain_{0}.itp\"\n".format(letter))
            file.write('\n')

            file.write('[ system ]\n')
            file.write('{0}\n\n'.format(universe.get('d_pdbName')))

            file.write('[ molecules ]\n')
            file.write('; Compounts \t\t #mols\n')
            for letter in universe.get('d_chain'):
                file.write("Protein_chain_{0}\t\t1\n".format(letter))

    # ADD RELEVANT PARAMETERS TO UNIVERSE ######################################
    universe.add('d_modelFF', d_modelFF)
    universe.add('d_modelWater', d_modelWater)
    universe.add('d_terministring', d_terministring)

    # PRE-TREAT THE .PDB FILE ##################################################
    lambdaTypeNames = []
    for lambdaType in universe.get('ph_lambdaTypes'):
        lambdaTypeNames.append(lambdaType.d_resname)

    lambdaTypeBaseNames = []
    for lambdaTypeName in lambdaTypeNames:
        lambdaTypeBaseNames.append(lambdaTypeName[0:3])

    # print("lambdaTypeNames", lambdaTypeNames)         # debug
    # print("LambdaTypeBaseNames", lambdaTypeBaseNames) # debug

    if (universe.get("ph_constantpH")):
        # If we use our default force field:
        if (d_modelFF == "charmm36-mar2019"):
            utils.update("generate", "using our default ({}) force field...".format(d_modelFF))
        # If we use our custom / modified force field:
        elif (d_modelFF == "charmm36-mar2019-m4"):
            utils.update("generate", "using our custom/modified ({}) force field...".format(d_modelFF))
        # Warn user if we use something different than charmm36-mar2019
        else:
            utils.warning("generate", "using an unknown ({}) force field!".format(d_modelFF))
            utils.warning("generate", "pHbuilder was made for charmm36-mar2019(-m4). All bets are off...".format(d_modelFF))

        residues = universe.get('d_residues')

        for residue in residues:
            if (residue.d_resname)[0:3] in lambdaTypeBaseNames:
                residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname[0:3])]

        # Update d_residues and overwrite curren _PR1.pdb
        universe.add('d_residues', residues)
        protein.write("{}_PR1.pdb".format(universe.get('d_pdbName')))

    # USER UPDATE STUFF ########################################################
    countACID = 0
    for resname in lambdaTypeNames:
        countACID += protein.countRes(resname)

    # If constant-pH is on,
    if universe.get('ph_constantpH'):
        utils.update("generate", 'constant-pH is turned on...')
        
        # and we have at least one protonatable reside,
        if countACID > 0:
            utils.update("generate", "detected {} acidic residue(s):".format(countACID))

            for letter in universe.get('d_chain'):
                count = 0

                for residue in universe.get('d_residues'):
                    if residue.d_chain == letter:
                        count += 1

                        if residue.d_resname in lambdaTypeNames:
                            utils.update("generate", "{:3s}-{:<4d} in chain {}".format(residue.d_resname, count, letter))

            utils.update("generate", "(setting protonation state to True (option 1) for all of these)")

        else:
            utils.warning("generate", "No acidic residues detected. Did you update the lambdaTypeName(s)? Turning off constant-pH...")
            universe.add('ph_constantpH', False)
            universe.add('ph_QQleveling', 0)

    else:
        utils.update("generate", 'constant-pH is turned off...')
        universe.add('ph_QQleveling', 0) # If ph_constantpH is False then this is also False.

    utils.update("generate", "using the {} force field with the {} water model...".format(d_modelFF, d_modelWater))

    if (d_terministring != ""):
        utils.update("generate", "Using options {} for termini...".format(d_terministring))
    else:
        utils.warning("generate", "No termini specified, using gmx default (00 = NH3+ and COO-)...")
    
    utils.update("generate", "running pdb2gmx to create {}_PR2.pdb and topol.top...".format(universe.get('d_pdbName')))

    # RUN ACTUAL PDB2GMX COMMAND ###############################################

    # If constant-pH is turned on AND we have a lambdaType ASPH or GLUH,
    # this means we need to use pdb2gmx interactive to set the protonation
    # states:
    options = ""
    if (universe.get('ph_constantpH')):
        if ("ASPH" in lambdaTypeNames):
            options += "-asp "
        if ("GLUH" in lambdaTypeNames):
            options += "-glu "

    # If the user specified a terministring, add -ter to pdb2gmx.
    if (d_terministring != ""):
        options += "-ter"

    xstr = "<< EOF"
    for chain in universe.get('d_chain'):
        # If constant-pH is turned on AND we have a lambdaType ASPH or GLUH, set
        # protonation state of these to True.
        if (universe.get('ph_constantpH')):
            for residue in universe.get('d_residues'):
                if residue.d_resname in ["ASPH", "GLUH"] and residue.d_chain == chain:
                    xstr += "\n1"

        # If the user specified a terministring, add those to the EOFstring as well:
        if (d_terministring != ""):
            xstr += "\n{}".format(d_terministring[0])
            xstr += "\n{}".format(d_terministring[1])
    
    xstr += "\nEOF" # End EOFstring.
    # print(options) # debug
    # print(xstr)    # debug

    os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -ff {2} -water {3} -ignh {4} >> builder.log 2>&1 {5}".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_modelFF, d_modelWater, options, xstr))

    # WRAPUP ###################################################################

    # Rebuild topology.
    rebuild_topol()

    # To update d_residues.
    protein.load("{0}_PR2.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_PR2.pdb".format(universe.get('d_pdbName')))

def restrain_dihedrals(resName, atomNameList, Type, phi, dphi, fc):
    utils.update("restrain_dihedrals", "will add restraints for {0} (all chains)...".format(resName))

    # Every chain has its own .itp file, so we loop through every file:
    for letter in universe.get('d_chain'):
        # This is to make sure we don't have multiple headers when we add multiple different restraints.
        first = False
        if not "[ dihedral_restraints ]" in open("topol_Protein_chain_{0}.itp".format(letter)).read():
            first = True

        # Append to the end of relevant .itp file:
        with open("topol_Protein_chain_{0}.itp".format(letter), 'a') as file:
            if first:
                file.write("[ dihedral_restraints ]\n")
                file.write("; ai aj ak al type phi dphi fc\n")

            # Atomcount resets for every separate .itp file.
            atomCount = 0
            for residue in universe.get('d_residues'):
                # Dictionary resets for every residue, as we can of course have
                # multiple ASPs or GLUs in one chain.
                dictionary = {}
                for atom in residue.d_atoms:
                    # Only increase atomcount when we read the relevant chain:
                    if residue.d_chain == letter:
                        atomCount += 1

                        if residue.d_resname == resName and atom in atomNameList:
                            dictionary[atom] = atomCount

                if len(dictionary) == 4:
                    utils.update("restrain_dihedrals", "adding restraints for chain {0} {1}-{2}...".format(residue.d_chain, resName, residue.d_resid))
                    for atom in atomNameList:
                        file.write("{:<6d} ".format(dictionary[atom]))
                    file.write(" {}  {}  {}  {}\n".format(Type, phi, dphi, fc))

def restrain_dihedrals_by_idx(indices, Type, phi, dphi, fc):
    utils.update("restrain_dihedrals", "will restrain dihedral {} {} {} {} (Type={}, phi={}, dphi={}, fc={})".format(indices[0], indices[1], indices[2], indices[3], Type, phi, dphi, fc))

    # Every chain has its own .itp file, so we loop through every file:
    for letter in universe.get('d_chain'):
        # This is to make sure we don't have multiple headers when we add multiple different restraints.
        first = False
        if not "[ dihedral_restraints ]" in open("topol_Protein_chain_{0}.itp".format(letter)).read():
            first = True

        # Append to the end of relevant .itp file:
        with open("topol_Protein_chain_{0}.itp".format(letter), 'a') as file:
            if first:
                file.write("[ dihedral_restraints ]\n")
                file.write("; ai aj ak al type phi dphi fc\n")

            for idx in indices:
                file.write("{:<6d} ".format(idx))
                
            file.write(" {}  {}  {}  {}\n".format(Type, phi, dphi, fc))
