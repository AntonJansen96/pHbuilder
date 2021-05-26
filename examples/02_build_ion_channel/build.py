#!/bin/python3

import phbuilder # Note: preparing this simulation might take a while.

# Set some pH-related data members in universe:
phbuilder.universe.add('ph_constantpH', True)
phbuilder.universe.add('ph_QQleveling', 2)

phbuilder.universe.defineLambdaType(
    resname = 'ASP', 
    pKa     = 3.65,
    atoms   = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'],
    qqA     = [-0.21, 0.75, -0.55, -0.61, 0.44],
    qqB     = [-0.28, 0.62, -0.76, -0.76, 0.00],
    dvdl    = [49.066, -563.785, -290.092, 1313.046, -2439.157, 2010.080, -647.785]
    )

phbuilder.universe.defineLambdaType(
    resname = 'GLU', 
    pKa     = 4.25,
    atoms   = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'],
    qqA     = [-0.21, 0.75, -0.55, -0.61, 0.44],
    qqB     = [-0.28, 0.62, -0.76, -0.76, 0.00],
    dvdl    = [26.339, -535.805, -78.575, -472.196, 1744.287, -1927.134, 655.915]
    )

phbuilder.universe.add('ph_BUF_dvdl', [672.41, -702.45, -63.10, 695.67, -1214.43, 537.14])

################################################################################

phbuilder.protein.process("../../proteins/4hfi.pdb")

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019", "tip3p", d_terministring="11")

phbuilder.protein.add_box(d_boxMargin=1.5, d_boxType='triclinic')
phbuilder.protein.add_buffer("../../proteins/buffer.pdb", "../../proteins/buffer.itp", ph_bufqqA=[-0.0656, 0.5328, 0.5328], ph_bufqqB=[-0.8476, 0.4238, 0.4238], ph_bufMargin=1.5, attempts=200000)
phbuilder.protein.add_water()
phbuilder.protein.add_ions()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()

phbuilder.md.gen_mdp('MD', nsteps=25000000, nstxout=10000)
phbuilder.md.gen_constantpH(ph_pH=4.0, ph_lambdaM=5.0, ph_nstout=500, ph_barrierE=7.5)
phbuilder.write.run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")

phbuilder.universe.inspect()
