import os
import universe

# Prints an update for the user.
def update(tool, message):
    print("{:18s} : {:s}".format(tool, message))

# Prints a warning for the user.
def warning(tool, message):
    print("{:18s} : WARNING - {:s}".format(tool, message))

# Generates the index.ndx.
def generate_index():
    os.system("gmx make_ndx -f {0} >> builder.log 2>&1 << EOF\nq\nEOF".format(universe.get('d_nameList')[-1]))

# Adds a .pdb file name to the namelist.
def add_to_nameList(name):
    if universe.has('d_nameList'):
        temp = universe.get('d_nameList')
        temp.append(name)
        universe.add('d_nameList', temp)
    else:
        universe.add('d_nameList', [name])
