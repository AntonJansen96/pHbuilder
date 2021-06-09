#!/bin/python3

import phbuilder

phbuilder.universe.add('ph_constantpH', True)
phbuilder.universe.add('ph_QQleveling', 0)

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

phbuilder.universe.add('ph_BUF_dvdl', [837.234, -888.419, -70.346, -402.684, 1031.590, -547.705]) # m4_ion_cal

################################################################################

phbuilder.protein.process('/home/anton/GIT/phbuilder/proteins/ASP_tri.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("/home/anton/GIT/phbuilder/ffields/charmm36-mar2019-m4.ff", "tip3p", d_terministring="34")

phbuilder.protein.add_box(d_boxMargin=1.8)
phbuilder.protein.add_buffer()
phbuilder.protein.add_water()
phbuilder.protein.add_ions()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()

phbuilder.md.gen_mdp('MD', nsteps=50000, nstxout=10000, sameSeed=True)
phbuilder.write.run(gmxPath="/usr/local/gromacs_test2")
phbuilder.md.gen_constantpH(ph_pH=3.65, ph_lambdaM=5.0, ph_nstout=100, ph_barrierE=0.0)
