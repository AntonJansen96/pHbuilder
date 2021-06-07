#!/bin/python3

# phbuilder can also be used to quickly build "normal" simulations, 
# with as little commands as:

import phbuilder 

phbuilder.universe.add('ph_constantpH', False)

phbuilder.protein.process('../../proteins/1cvo.pdb')
phbuilder.topol.generate("../../ffields/charmm36-mar2019-m4.ff", "tip3p")

phbuilder.protein.add_box(d_boxMargin=1.0)
phbuilder.protein.add_water()
phbuilder.protein.add_ions()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()

phbuilder.md.gen_mdp('MD', nsteps=50000, nstxout=500)
phbuilder.write.run()
phbuilder.write.reset()
