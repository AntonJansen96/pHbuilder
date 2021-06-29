#!/bin/python3

import phbuilder

# Set some pH-related data members in universe:
phbuilder.universe.add('ph_constantpH', True)
phbuilder.universe.add('ph_QQleveling', 2)

phbuilder.universe.defineLambdaType(
    resname = 'ASPT', 
    pKa     = 3.65,
    atoms   = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'],
    qqA     = [-0.21, 0.75, -0.55, -0.61, 0.44],
    qqB     = [-0.28, 0.62, -0.76, -0.76, 0.00],
    dvdl    = [44.621, -554.27, -146.03, 278.55, -144.28, -54.078] # m4_nooradV/dl
    )

phbuilder.universe.defineLambdaType(
    resname = 'GLUT', 
    pKa     = 4.25,
    atoms   = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'],
    qqA     = [-0.21, 0.75, -0.55, -0.61, 0.44],
    qqB     = [-0.28, 0.62, -0.76, -0.76, 0.00],
    dvdl    = [27.157, -558.21, -204.91, 514.37, -429.77, 59.636] # m4_nooradV/dl
    )

phbuilder.universe.add('ph_BUF_dvdl', [672.41, -702.45, -63.10, 695.67, -1214.43, 537.14]) # m4_water_cal

################################################################################

phbuilder.protein.process('../../proteins/ASP_tri.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("../../ffields/charmm36-mar2019-m4.ff", "tip3p", d_terministring="34")

phbuilder.protein.add_box(d_boxMargin=2.0)
phbuilder.protein.add_buffer(ph_bufpdbName="../../proteins/buffer.pdb", ph_bufitpName="../../proteins/buffer.itp", ph_bufqqA=[-0.0656, 0.5328, 0.5328], ph_bufqqB=[-0.8476, 0.4238, 0.4238])
phbuilder.protein.add_water()
phbuilder.protein.add_ions()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()

phbuilder.md.gen_mdp('MD', nsteps=100000, nstxout=10000)
phbuilder.md.gen_constantpH(ph_pH=3.65, ph_lambdaM=5.0, ph_nstout=500, ph_barrierE=7.5)
phbuilder.write.run(gmxPath="/usr/local/gromacs_constantph")

phbuilder.universe.inspect()
