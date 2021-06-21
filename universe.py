import shelve, utils

# Set/update variable to universe.
def add(varName, value):
    with shelve.open('universe') as shelf:
        shelf[varName] = value

# Check whether universe contains a certain varName
def has(varName):
    with shelve.open('universe') as shelf:
        return varName in shelf

# Retrieve variable from universe.
def get(varName):
    if has(varName):
        return shelve.open('universe')[varName]

    data = eval(input("couldn't retrieve var \"{0}\" from universe. Enter manually: ".format(varName)))
    print("add {0} = {1} {2}".format(varName, data, type(data)))
    add(varName, data)
    return data

# Display all variables (name, data, type) stored in the universe.
def inspect():
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in sorted(shelf):
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}] ({3}) {4}".format(item, shelf[item][0], shelf[item][-1], len(shelf[item]), type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1} {2}".format(item, shelf[item], type(shelf[item]), arg=longest))

# Stores the information of a lambda residue-type.
class LambdaType:
    def __init__(self, resname, pKa, atoms, qqA, qqB, dvdl):
        self.d_resname = resname
        self.d_pKa     = pKa
        self.d_atoms   = atoms
        self.d_qqA     = qqA
        self.d_qqB     = qqB
        self.d_dvdl    = dvdl

# Add a lambda residue-type to universe.
def defineLambdaType(resname, pKa, atoms, qqA, qqB, dvdl):
    if (len(resname) != 4):
        resname = input("Name of lambdaType \"{}\" must be four letters (e.g. ASPH, GLUT): ".format(resname))

    NewLambdaType = LambdaType(resname, pKa, atoms, qqA, qqB, dvdl)
    if has('ph_lambdaTypes'):
        temp = get('ph_lambdaTypes')
        
        for entry in temp:
            if entry.d_resname == NewLambdaType.d_resname:
                utils.update("defineLambdaType", "LambdaType {} is already defined in ph_lambdaTypes. Skipping...".format(NewLambdaType.d_resname))
                break
        else:
            temp.append(NewLambdaType)
            add('ph_lambdaTypes', temp)
    else:
        add('ph_lambdaTypes', [NewLambdaType])
