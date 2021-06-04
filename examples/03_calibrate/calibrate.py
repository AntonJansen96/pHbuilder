#!/bin/python3

import phbuilder
import os, numpy, load

# Set some pH-related data members in universe:
phbuilder.universe.add('ph_constantpH', True)
phbuilder.universe.add('ph_QQleveling', 0)

phbuilder.universe.defineLambdaType(
    resname = 'ASPT', 
    pKa     = 3.65,
    atoms   = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'],
    qqA     = [-0.21, 0.75, -0.55, -0.61, 0.44],
    qqB     = [-0.28, 0.62, -0.76, -0.76, 0.00],
    dvdl    = [1, 2, 3]
    )

phbuilder.universe.defineLambdaType(
    resname = 'GLUT', 
    pKa     = 4.25,
    atoms   = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'],
    qqA     = [-0.21, 0.75, -0.55, -0.61, 0.44],
    qqB     = [-0.28, 0.62, -0.76, -0.76, 0.00],
    dvdl    = [1, 2, 3]
    )

phbuilder.universe.add('ph_BUF_dvdl', [1, 2, 3])

################################################################################

phbuilder.protein.process('../../proteins/GLU_tri.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019-m4", "tip3p", d_terministring="34")

phbuilder.protein.add_box(d_boxMargin=2.0)
phbuilder.protein.add_buffer("../../proteins/buffer.pdb", "../../proteins/buffer.itp", ph_bufqqA=[-0.0656, 0.5328, 0.5328], ph_bufqqB=[-0.8476, 0.4238, 0.4238])
phbuilder.protein.add_water()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()
phbuilder.write.run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")

# The part where we do the loop to get the mean and standard deviations:

dVdlInitList = [ii / 10.0 for ii in range(-1, 12)]
dVdlMeanList = []
dVdlStdList  = []

for init in dVdlInitList:
    phbuilder.md.energy_tcouple()
    phbuilder.md.energy_pcouple()
    phbuilder.md.gen_mdp('MD', nsteps=50000, nstxout=10000)
    phbuilder.md.gen_constantpH(ph_pH=4.25, ph_lambdaM=0.0, ph_nstout=1, ph_barrierE=0.0, cal=True, lambdaInit=init)

    os.system("./run.sh")

    dVdlList = load.Col('lambda_1.dat', 3)
    dVdlMeanList.append(numpy.mean(dVdlList))
    dVdlStdList.append(numpy.std(dVdlList))

phbuilder.universe.add('ph_dvdl_initList', dVdlInitList)
phbuilder.universe.add('ph_dvdl_meanList', dVdlMeanList)
phbuilder.universe.add('ph_dvdl_stdList', dVdlStdList)

phbuilder.universe.inspect()
phbuilder.analyze.fitCalibration(order=5)
