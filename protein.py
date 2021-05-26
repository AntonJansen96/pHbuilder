import os
import universe, utils, topol

class Residue: # Stores a single residue's data.
    def __init__(self, atoms, ali, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms      # list      holds atom names
        self.d_ali     = ali        # list      holds alternative loc. ind.
        self.d_resname = resname    # string    holds residue name
        self.d_chain   = chain      # string    holds chain name (A, B)
        self.d_resid   = resid      # int       holds residue number
        self.d_x       = x          # list      holds x-coordiantes
        self.d_y       = y          # list      holds y-coordinates
        self.d_z       = z          # list      holds z-coordinatees

def process(fname, d_model=1, d_ALI='A', d_chain=[], resetResId=False):
    basename = os.path.basename(fname)
    universe.add('d_pdbName', basename[0:len(basename)-4])
    universe.add('d_model', d_model)
    universe.add('d_ALI', d_ALI)

    load(fname, d_model, d_ALI, d_chain)

    # Update d_chain from [] to a list of actual chains used:
    d_chain = []                                    
    for residue in universe.get('d_residues'):
        d_chain.append(residue.d_chain)
    d_chain = list(set(d_chain))
    d_chain.sort()
    universe.add('d_chain', d_chain)

    utils.update("process", 'file={0}, MODEL#={1}, ALI={2}, chain(s)={3}...'.format(fname, d_model, d_ALI, d_chain))

    # Write processed.pdb to file:
    utils.update("process", "writing processed .pdb file to {}_PR1.pdb...".format(universe.get('d_pdbName')))
    write("{0}_PR1.pdb".format(universe.get('d_pdbName')))

# Load a .pdb file into the universe.
def load(fname, d_model=1, d_ALI='A', d_chain=[]):
    d_residues = []
    atomLines  = []
    with open(fname) as file:       # Read .pdb line-by-line.
        read = True                 # True if no specific MODEL specified.

        for line in file.readlines():
            # This is to make sure we only import the specified MODEL.
            if ((line[0:6]) == "MODEL "):
                if ("MODEL {:8d}".format(d_model) in line):
                    read = True
                else:
                    read = False
            
            # Get title.
            if ((line[0:6]) == "TITLE "):
                d_title = line[7:80].rstrip(); universe.add('d_title', d_title)

            # Get periodic box information (if any).
            if ((line[0:6]) == "CRYST1"):
                d_box   = line[7:80].rstrip(); universe.add('d_box', d_box)
                
            # if our line specifies an ATOM,
            if (line[0:6] == "ATOM  "):                 
                # and we are currently reading the correct MODEL,
                if (read == True):
                    # and our line contains the correct specified alternate 
                    # location specifier (if any at all),
                    if (line[16:17] in [d_ALI, " "]):
                        # and we want all chains,
                        if (d_chain == []):
                            # then load the atom.
                            atomLines.append(line)
                        # or if we want a selection of chains,
                        elif (line[21:22] in d_chain):
                            # load that selection.
                            atomLines.append(line)

    # Add one line of padding to prevent IndexError
    atomLines.append("0000000000000000000000000000000000000000000000000000")

    # Loop through lines and create list of Residue objects.
    atoms = []; ali = []; xCoord = []; yCoord = []; zCoord = []
    
    for idx in range(0, len(atomLines) - 1):
        atoms.append(atomLines[idx][12:16])
        ali.append(atomLines[idx][16:17])
        xCoord.append(float(atomLines[idx][30:38]))
        yCoord.append(float(atomLines[idx][38:46]))
        zCoord.append(float(atomLines[idx][46:54]))

        # If the resid of the next line is different, we are at end
        if (atomLines[idx][22:26] != atomLines[idx + 1][22:26]):
            # put the data in a Residue object and append to d_residues:
            d_residues.append(
                Residue
                (
                    atoms, 
                    ali,
                    atomLines[idx][17:20], 
                    atomLines[idx][21:22], 
                    int(atomLines[idx][22:26]), 
                    xCoord, 
                    yCoord,
                    zCoord
                ))
            # Reset.
            atoms = []; ali = []; xCoord = []; yCoord = []; zCoord = []

    universe.add('d_residues', d_residues)

# Write (modified) .pdb file.
def write(name):            
    with open(name, 'w') as file:
        if universe.has('d_title'):
            file.write("TITLE {0}\n".format(universe.get('d_title')))
        
        if universe.has('d_box'):
            file.write("CRYST1{0}\n".format(universe.get('d_box')))
        
        file.write("MODEL {:8d}\n".format(universe.get('d_model')))

        atomNumber = 1
        for residue in universe.get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

    utils.add_to_nameList(name)

# Count number of residues with a specific residue name.
def countRes(resName):
    count = 0
    for residue in universe.get('d_residues'):
        if (residue.d_resname == resName):
            count += 1
    
    return count

def add_box(d_boxMargin, d_boxType='cubic'):
    utils.update("add_box", "adding box using gmx editconf (boxMargin = {0}, boxType = {1})...".format(d_boxMargin, d_boxType))

    os.system("gmx editconf -f {0} -o {1}_BOX.pdb -d {2} -bt {3} >> builder.log 2>&1".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_boxMargin, d_boxType))

    # To set d_boxMargin and d_boxType.
    universe.add('d_boxMargin', d_boxMargin)
    universe.add('d_boxType', d_boxType)

    # To update d_box.
    load("{0}_BOX.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_BOX.pdb".format(universe.get('d_pdbName')))

def add_buffer(ph_bufpdbName="", ph_bufitpName="", ph_bufqqA=[1], ph_bufqqB=[0], ph_bufMargin=2.0, attempts=100000):
    # This function writes a dummy .pdb containing the default buffer (ion).
    def writeDefaultPDB():
        with open("defaultBuffer.pdb", 'w') as file:
            file.write("TITLE     BUFFER PARTICLE\n")
            file.write("MODEL        1\n")
            file.write("ATOM      1  NA  BUF     1     110.896   2.872  68.855  1.00  0.00            \n")
            file.write("TER\n")
            file.write("ENDMDL\n")

    # This function writes the topology for the default buffer (ion).
    # Note: charge should not be 0 because then some interactions are not generated???
    def writeDefaultITP():
        with open("defaultBuffer.itp", 'w') as file:
            file.write("[ moleculetype ]\n")
            file.write("; molname	nrexcl\n")
            file.write("BUF		1\n\n")
            file.write("[ atoms ]\n")
            file.write("; id	at type		res nr	residu name at name  cg nr	charge	 \n")
            file.write("1		SOD			1		   BUF			NA		   1		0	 \n\n")
            file.write("#ifdef POSRES_BUF\n")
            file.write("; Position restraint for each buffer ion\n")
            file.write("[ position_restraints ]\n")
            file.write(";  i funct       fcx        fcy        fcz\n")
            file.write("   1    1       1000       1000       1000\n")
            file.write("#endif\n")

    # Skip this whole step if we don't need it.
    if not (universe.get('ph_constantpH') and universe.get('ph_QQleveling') in [1, 2]):
        utils.update("add_buffer", "either ph_constantpH is False or ph_QQleveling = 0 --> skipping...")
        return

    # Make sure that the sum of the charges in state A and B are correct.
    if (sum(ph_bufqqA) != 1 or sum(ph_bufqqB) != 0):
        utils.update("add_buffer", "buffer charges incorrectly specified! sums must be 1 and 0")
        universe.get('ph_bufqqA')
        universe.get('ph_bufqqB')

    # Determine whether we use the default or a custom buffer.
    useDefault = False
    if (ph_bufpdbName == "" and ph_bufitpName == ""):
        utils.update("add_buffer", "using default (builtin) buffer...")
        useDefault = True
    elif (ph_bufpdbName == "" and ph_bufitpName != ""):
        utils.update("add_buffer", "ph_bufitpName not specified, resorting to default buffer!")
        useDefault = True
    elif (ph_bufpdbName != "" and ph_bufitpName == ""):
        utils.update("add_buffer", "ph_bufpdbName not specified, resorting to default buffer!")
        useDefault = True

    if (useDefault):
        # Check to make sure that the charges for the default buffer are correct.
        if (ph_bufqqA != [1] or ph_bufqqB != [0]):
            utils.update("add_buffer", "buffer charges incorrectly specified for default buffer!")
            universe.get('ph_bufqqA')
            universe.get('ph_bufqqB')

        # Generate the files for the default buffer and update data members.
        writeDefaultPDB(); ph_bufpdbName = "defaultBuffer.pdb"
        writeDefaultITP(); ph_bufitpName = "defaultBuffer.itp"

    # Get the number of buffer molecules we need.
    ph_bufnmol = 0
    for lambdaType in universe.get('ph_lambdaTypes'):
        ph_bufnmol += countRes(lambdaType.d_resname)

    utils.update("add_buffer", "will attempt to add {0} buffer molecule(s)...".format(ph_bufnmol))

    # RUN GROMACS INSERT-MOLECULES COMMAND
    os.system("touch vdwradii.dat") # we need this dummy file for this to work.

    os.system("gmx insert-molecules -f {0} -o {1}_BUF.pdb -ci {2} -nmol {3} -scale 1.0 -radius {4} -try {5} >> builder.log 2>&1".format(
        universe.get('d_nameList')[-1],
        universe.get('d_pdbName'),
        ph_bufpdbName,
        ph_bufnmol,
        0.5 * ph_bufMargin,
        int(attempts / ph_bufnmol)))

    os.remove("vdwradii.dat") # clean dummy file.

    # To update d_residues.
    load("{0}_BUF.pdb".format(universe.get('d_pdbName')))    

    # Give user a warning if there wasn't enough space.
    actual = countRes('BUF')
    if actual < ph_bufnmol:
        utils.update("add_buffer", "warning: only {0}/{1} requested buffer molecules inserted after {2} attempts,".format(actual, ph_bufnmol, attempts))
        utils.update("add_buffer", "warning: try decreasing ph_bufMargin (={0}nm) or increasing d_boxMargin (={1}nm)...".format(ph_bufMargin, universe.get('d_boxMargin')))
    else:
        utils.update("add_buffer", "succesfully added {0} buffer molecule(s)...".format(actual))

    # To add buffer topology to topol.top.
    utils.update("add_buffer", "updating topology...")

    if (useDefault):
        os.remove("defaultBuffer.pdb")              # Remove dummy .pdb file.
    else:
        os.system("cp {} .".format(ph_bufitpName))  # Copy to working dir.

    topol.add_mol(os.path.basename(ph_bufitpName), "Include buffer topology", 'BUF', actual)

    # Set some parameters in the universe.
    universe.add('ph_bufpdbName', ph_bufpdbName)
    universe.add('ph_bufitpName', ph_bufitpName)
    universe.add('ph_bufqqA', ph_bufqqA)
    universe.add('ph_bufqqB', ph_bufqqB)
    universe.add('ph_bufMargin', ph_bufMargin)
    universe.add('ph_bufnmol', actual)

    # To update d_nameList.
    utils.add_to_nameList("{0}_BUF.pdb".format(universe.get('d_pdbName')))

def add_water():
    utils.update("add_water", "running gmx solvate...")
    
    os.system("gmx solvate -cp {0} -o {1}_SOL.pdb >> builder.log 2>&1".format(universe.get('d_nameList')[-1], universe.get('d_pdbName')))

    # To update d_residues.
    load("{0}_SOL.pdb".format(universe.get('d_pdbName')))

    # To update topol.top.
    topol.add_mol("{0}.ff/{1}.itp".format(universe.get('d_modelFF'), universe.get('d_modelWater')), "Include water topology", "SOL", countRes('SOL'))

    # To update d_nameList.
    utils.add_to_nameList("{0}_SOL.pdb".format(universe.get('d_pdbName')))

def add_ions(neutral=True, conc=0, pname='NA', nname='CL'):
    if (not neutral) and (conc == 0):
        utils.update("add_ions", "no ions will be added...")
        return

    if neutral:
        utils.update("add_ions", "will add ions ({}/{}) to neutralize the system...".format(pname, nname))

    if conc > 0:
        utils.update("add_ions", "will add ions ({}/{}) for target concentration = {} mmol/ml...".format(pname, nname, conc))

    # Generate IONS.mdp (just a dummy required).
    os.system('touch IONS.mdp')

    # Add ion topology to topol.top.
    topol.add_mol("{0}.ff/ions.itp".format(universe.get('d_modelFF')), "Include ion topology")

    # Run gmx grompp and genion.
    utils.update("add_ions", "running gmx grompp and genion to add ions...")

    os.system("gmx grompp -f IONS.mdp -c {0} -p topol.top -o IONS.tpr >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))

    if neutral:
        os.system("gmx genion -s IONS.tpr -o {0}_ION.pdb -p topol.top -pname {1} -nname {2} -conc {3} -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF".format(universe.get('d_pdbName'), pname, nname, conc))
    else:
        os.system("gmx genion -s IONS.tpr -o {0}_ION.pdb -p topol.top -pname {1} -nname {2} -conc {3} >> builder.log 2>&1 << EOF\nSOL\nEOF".format(universe.get('d_pdbName'), pname, nname, conc))

    # To update d_residues.
    load("{0}_ION.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_ION.pdb".format(universe.get('d_pdbName')))
