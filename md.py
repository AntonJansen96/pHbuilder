import os, universe, utils

from mdp import gen_mdp
from constantph import gen_constantpH

def energy_minimize():
    gen_mdp('EM')

    utils.update("energy_minimize", "running gmx grompp and mdrun for energy minimization...")

    os.system("gmx grompp -f EM.mdp -c {0} -p topol.top -o EM.tpr -r {0} >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx mdrun -deffnm EM -c EM.pdb >> builder.log 2>&1")

    utils.add_to_nameList("EM.pdb")

def energy_tcouple():
    gen_mdp('NVT')

    utils.update("energy_tcouple", "running gmx grompp and mdrun for temperature coupling...")

    os.system("gmx grompp -f NVT.mdp -c {0} -p topol.top -o NVT.tpr -r {0} >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx mdrun -deffnm NVT -c NVT.pdb >> builder.log 2>&1")

    utils.add_to_nameList("NVT.pdb")

def energy_pcouple():
    gen_mdp('NPT')

    utils.update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")

    os.system("gmx grompp -f NPT.mdp -c {0} -p topol.top -o NPT.tpr -r {0} >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx mdrun -deffnm NPT -c NPT.pdb >> builder.log 2>&1")

    utils.add_to_nameList("NPT.pdb")
