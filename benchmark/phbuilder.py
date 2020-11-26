#!/bin/python3

import os, sys
from lib import sim

# Set up environment/external files.
gromPath = "/home/anton/GIT/phbuilder/grom"
os.system("cp {0}/buffer.itp {0}/buffer.pdb {0}/IONS.mdp ./".format(gromPath))

################################################################################

bool_pH = bool(int(sys.argv[1]))
opMode  = sys.argv[2]
nsteps  = int(sys.argv[3])

# print(bool_pH)    # debug
# print(opMode)     # debug
# print(nsteps)     # debug

sim = sim()
sim.setconstantpH(value=bool_pH, restrain=False)
sim.processpdb("1cvo.pdb")

sim.protein_add_forcefield("charmm36-mar2019", "tip3p")
sim.protein_add_box(boxSizeMargin=1.0)
sim.protein_add_buffer(minSep=1.5)
sim.protein_add_water()
sim.protein_add_ions()

sim.generate_index()

sim.generate_mdp('EM')
sim.generate_mdp('NVT')
sim.generate_mdp('NPT')
sim.generate_mdp('MD', nsteps=nsteps, nstxout=0)

sim.generate_phdata(pH=4.5, lambdaM=5.0, nstOut=10000, barrierE=5.0)
# sim.generate_phdata_legacy(pH=4.5, lambdaM=5.0, nstOut=10000, barrierE=5.0)

sim.write_run("/usr/local/gromacs_test2", mode=opMode)
sim.write_reset()

sim.energy_minimize()
sim.energy_tcouple()
sim.energy_pcouple(skip=True)