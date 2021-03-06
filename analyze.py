import os
import matplotlib.pyplot as plt
import numpy as np
import universe, load

def titrate(lambdaFileName, cutoff=0.80):
    lambda_proto   = 0
    lambda_deproto = 0
    
    for x in load.Col(lambdaFileName, 2):
        if (x > cutoff):
            lambda_deproto += 1
        if (x < 1 - cutoff):
            lambda_proto   += 1

    fraction = float(lambda_deproto) / (lambda_proto + lambda_deproto)

    return fraction

def plotlambda(plotBUF=False):
    resnameList = []    # Get the names and such of all the ASPs and GLUs.
    residList   = []
    for residue in universe.get('d_residues'):
        if residue.d_resname in ["ASP", "ASPH", "ASPT", "GLU", "GLUH", "GLUT"]:
            resnameList.append(residue.d_resname)
            residList.append(residue.d_resid)

    plt.figure()
    for idx in range(1, len(resnameList) + 1):
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        plt.plot(t, x, label="%s-%s" % (resnameList[idx-1], residList[idx - 1]), linewidth=0.5)

    if (plotBUF):
        t = load.Col("lambda_{0}.dat".format(len(resnameList) + 1), 1)
        x = load.Col("lambda_{0}.dat".format(len(resnameList) + 1), 2)

        plt.plot(t, x, label="Buffer", linewidth=0.5)

    plt.xlabel("Time (ps)")
    plt.ylabel(r"$\lambda$-coordinate")
    plt.ylim(-0.1, 1.1)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))

    plt.legend()
    plt.grid()
    plt.show()

def glicphstates():
    # EXPERIMENTAL DATA ON PROTONATION STATES AT VARIOUS PH ####################
    biophys = { # also prevost2012
        'ASPT-13'  : 1,
        'ASPT-31'  : 1,
        'ASPT-32'  : 1,
        'ASPT-49'  : 1,
        'ASPT-55'  : 1,
        'ASPT-86'  : 0,
        'ASPT-88'  : 0,
        'ASPT-91'  : 1,
        'ASPT-97'  : 1,
        'ASPT-115' : 1,
        'ASPT-122' : 1,
        'ASPT-136' : 1,
        'ASPT-145' : 1,
        'ASPT-153' : 1,
        'ASPT-154' : 1,
        'ASPT-161' : 1,
        'ASPT-178' : 1,
        'ASPT-185' : 1,
        'GLUT-14'  : 1,
        'GLUT-26'  : 0,
        'GLUT-35'  : 0,
        'GLUT-67'  : 0,
        'GLUT-69'  : 1,
        'GLUT-75'  : 0,
        'GLUT-82'  : 0,
        'GLUT-104' : 1,
        'GLUT-147' : 1,
        'GLUT-163' : 1,
        'GLUT-177' : 0,
        'GLUT-181' : 1,
        'GLUT-222' : 1,
        'GLUT-243' : 0,
        'GLUT-272' : 1,
        'GLUT-282' : 1
    }

    nury2010 = { # this is also cheng2010, calimet2013
        'ASPT-13'  : 1,
        'ASPT-31'  : 1,
        'ASPT-32'  : 1,
        'ASPT-49'  : 1,
        'ASPT-55'  : 1,
        'ASPT-86'  : 0,
        'ASPT-88'  : 0,
        'ASPT-91'  : 1,
        'ASPT-97'  : 1,
        'ASPT-115' : 1,
        'ASPT-122' : 1,
        'ASPT-136' : 1,
        'ASPT-145' : 1,
        'ASPT-153' : 1,
        'ASPT-154' : 1,
        'ASPT-161' : 1,
        'ASPT-178' : 1,
        'ASPT-185' : 1,
        'GLUT-14'  : 1,
        'GLUT-26'  : 0,
        'GLUT-35'  : 0,
        'GLUT-67'  : 0,
        'GLUT-69'  : 0,
        'GLUT-75'  : 0,
        'GLUT-82'  : 0,
        'GLUT-104' : 1,
        'GLUT-147' : 1,
        'GLUT-163' : 1,
        'GLUT-177' : 0,
        'GLUT-181' : 1,
        'GLUT-222' : 1,
        'GLUT-243' : 0,
        'GLUT-272' : 1,
        'GLUT-282' : 1
    }

    fritsch2011 = {
        'ASPT-13'  : 0,
        'ASPT-31'  : 0,
        'ASPT-32'  : 1,
        'ASPT-49'  : 1,
        'ASPT-55'  : 0,
        'ASPT-86'  : 0,
        'ASPT-88'  : 0,
        'ASPT-91'  : 0,
        'ASPT-97'  : 0,
        'ASPT-115' : 1,
        'ASPT-122' : 1,
        'ASPT-136' : 1,
        'ASPT-145' : 0,
        'ASPT-153' : 0,
        'ASPT-154' : 0,
        'ASPT-161' : 0,
        'ASPT-178' : 0,
        'ASPT-185' : 0,
        'GLUT-14'  : 0,
        'GLUT-26'  : 0,
        'GLUT-35'  : 0,
        'GLUT-67'  : 0,
        'GLUT-69'  : 0,
        'GLUT-75'  : 0,
        'GLUT-82'  : 0,
        'GLUT-104' : 1,
        'GLUT-147' : 0,
        'GLUT-163' : 0,
        'GLUT-177' : 0,
        'GLUT-181' : 0,
        'GLUT-222' : 1,
        'GLUT-243' : 0,
        'GLUT-272' : 0,
        'GLUT-282' : 0
    }

    lev2017 = {
        'ASPT-13'  : 1,
        'ASPT-31'  : 1,
        'ASPT-32'  : 1,
        'ASPT-49'  : 1,
        'ASPT-55'  : 1,
        'ASPT-86'  : 1,
        'ASPT-88'  : 1,
        'ASPT-91'  : 1,
        'ASPT-97'  : 1,
        'ASPT-115' : 1,
        'ASPT-122' : 1,
        'ASPT-136' : 1,
        'ASPT-145' : 1,
        'ASPT-153' : 1,
        'ASPT-154' : 1,
        'ASPT-161' : 1,
        'ASPT-178' : 1,
        'ASPT-185' : 1,
        'GLUT-14'  : 1,
        'GLUT-26'  : 0,
        'GLUT-35'  : 0,
        'GLUT-67'  : 0,
        'GLUT-69'  : 0,
        'GLUT-75'  : 0,
        'GLUT-82'  : 0,
        'GLUT-104' : 1,
        'GLUT-147' : 1,
        'GLUT-163' : 1,
        'GLUT-177' : 0,
        'GLUT-181' : 1,
        'GLUT-222' : 1,
        'GLUT-243' : 0,
        'GLUT-272' : 1,
        'GLUT-282' : 1
    }

    nemecz2017 = { # also Hu2018
        'ASPT-13'  : 1,
        'ASPT-31'  : 1,
        'ASPT-32'  : 1,
        'ASPT-49'  : 1,
        'ASPT-55'  : 1,
        'ASPT-86'  : 0,
        'ASPT-88'  : 0,
        'ASPT-91'  : 1,
        'ASPT-97'  : 1,
        'ASPT-115' : 1,
        'ASPT-122' : 1,
        'ASPT-136' : 1,
        'ASPT-145' : 1,
        'ASPT-153' : 1,
        'ASPT-154' : 1,
        'ASPT-161' : 1,
        'ASPT-178' : 1,
        'ASPT-185' : 1,
        'GLUT-14'  : 1,
        'GLUT-26'  : 0,
        'GLUT-35'  : 0,
        'GLUT-67'  : 1,
        'GLUT-69'  : 1,
        'GLUT-75'  : 1,
        'GLUT-82'  : 1,
        'GLUT-104' : 1,
        'GLUT-147' : 1,
        'GLUT-163' : 1,
        'GLUT-177' : 1,
        'GLUT-181' : 1,
        'GLUT-222' : 0,
        'GLUT-243' : 0,
        'GLUT-272' : 1,
        'GLUT-282' : 1
    }

    ullman = { # unpublished
        'ASPT-13'  : 1,
        'ASPT-31'  : 1,
        'ASPT-32'  : 1,
        'ASPT-49'  : 1,
        'ASPT-55'  : 1,
        'ASPT-86'  : 1,
        'ASPT-88'  : 1,
        'ASPT-91'  : 1,
        'ASPT-97'  : 1,
        'ASPT-115' : 1,
        'ASPT-122' : 1,
        'ASPT-136' : 1,
        'ASPT-145' : 1,
        'ASPT-153' : 1,
        'ASPT-154' : 1,
        'ASPT-161' : 1,
        'ASPT-178' : 1,
        'ASPT-185' : 1,
        'GLUT-14'  : 1,
        'GLUT-26'  : 0,
        'GLUT-35'  : 0,
        'GLUT-67'  : 0,
        'GLUT-69'  : 0,
        'GLUT-75'  : 0,
        'GLUT-82'  : 1,
        'GLUT-104' : 1,
        'GLUT-147' : 0,
        'GLUT-163' : 0,
        'GLUT-177' : 0,
        'GLUT-181' : 1,
        'GLUT-222' : 1,
        'GLUT-243' : 0,
        'GLUT-272' : 0,
        'GLUT-282' : 1
    }

    # DIRECTORY STRUCTURE ######################################################
    dirname = "lambdaplots" 
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        os.system("rm -f {0}/*.png {0}/*.pdf".format(dirname))

    # GET THE RESIDUE NUMBER, NAME, AND CHAIN OF ALL PROTO RESIDUES ############

    resnameList = []
    residList   = []
    chainList   = []
    for residue in universe.get('d_residues'):
        if residue.d_resname in ["ASPT", "GLUT"]:
            resnameList.append(residue.d_resname)
            residList.append(residue.d_resid)
            chainList.append(residue.d_chain)

    # CREATE LAMBDA PLOT FOR EVERY INDIVIDUAL PROTONATABLE RESIDUE #############

    # Loop through all the lambdas:
    for idx in range(1, len(resnameList) + 1):
        
        plt.figure(figsize=(8, 6))

        # Update user
        print("plotting {}/{}".format(idx, len(resnameList)), end='\r')
        
        # Load columns from .dat files
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        # If we had a crash, then the last line may not have been completely written.
        if len(t) > len(x):
            t.pop()
        elif len(t) < len(x):
            x.pop()

        plt.plot(t, x, linewidth=0.5)

        # Title
        plt.title("{0}-{1} in chain {2} in {3}.pdb\npH={4}, nstlambda={5}, deprotonation={6:.2f}".format(
            resnameList[idx - 1],
            residList[idx - 1],
            chainList[idx - 1],
            universe.get('d_pdbName'),
            universe.get('ph_pH'),
            universe.get('ph_nstout'),
            titrate("lambda_{}.dat".format(idx))
        ))

        # Axes and stuff
        plt.ylim(-0.1, 1.1)
        plt.xlabel("Time (ps)")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
        plt.ylabel(r"$\lambda$-coordinate")
        plt.grid()

        # Save.
        fileName = "{}/{}_{}-{:03d}".format(dirname, chainList[idx-1], resnameList[idx-1], residList[idx-1])
        # plt.savefig("{}.pdf".format(fileName)); os.system("pdfcrop {0}.pdf {0}.pdf >> /dev/null 2>&1".format(fileName))
        plt.savefig("{}.png".format(fileName))

        # clf = clear the entire current figure. close = closes a window.
        plt.clf(); plt.close()

    # CREATE HISTOGRAM PLOTS FOR COMBINED PROTO STATE OF ALL FIVE CHAINS #######
    number_of_chains   = len(set(chainList))
    residues_per_chain = int(len(resnameList) / number_of_chains)
    
    for ii in range(1, residues_per_chain + 1):
        data = []        
        for jj in range(0, number_of_chains):
            print(ii + residues_per_chain * jj, end=' ')
            # data += (load.Col('lambda_{}.dat'.format(ii + residues_per_chain * jj), 2, 49713, 124320))
            data += (load.Col('lambda_{}.dat'.format(ii + residues_per_chain * jj), 2))
        print()

        # PLOTTING STUFF #######################################################

        plt.figure(figsize=(8, 6))
        plt.hist(data, density=True, bins=200)
        
        # Title
        plt.title("{0}-{1} (all chains) in {2}.pdb\npH={3}, nstlambda={4}, deprotonation={5:.2f}".format(
            resnameList[ii-1],
            residList[ii-1],
            universe.get('d_pdbName'),
            universe.get('ph_pH'),
            universe.get('ph_nstout'),
            titrate("lambda_{}.dat".format(ii))
            # expVals40["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]
            ))

        # Axes and stuff
        plt.axis([-0.1, 1.1, -0.1, 12])
        plt.xlabel(r"$\lambda$-coordinate")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
        plt.grid()

        # Add green vertical line indicating experimental value
        plt.vlines(x=biophys["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=12, color='r', linewidth=4.0, label="biophysics.se/Prevost2012 = {}".format(biophys["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=nury2010["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=10, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013 = {}".format(nury2010["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=fritsch2011["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=8, color='b', linewidth=4.0, label="Fritsch2011 = {}".format(fritsch2011["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=lev2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=6, color='c', linewidth=4.0, label="Lev2017 = {}".format(lev2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=nemecz2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=4, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018 = {}".format(nemecz2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=ullman["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=2, color='y', linewidth=4.0, label="Ullman (unpublished) = {}".format(ullman["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))

        plt.legend()
        # Save and clear
        fileName = "{}/hist_{}-{:03d}".format(dirname, resnameList[ii-1], residList[ii-1])
        # plt.savefig("{}.pdf".format(fileName)); os.system("pdfcrop {0}.pdf {0}.pdf >> /dev/null 2>&1".format(fileName))
        plt.savefig('{}.png'.format(fileName))
        plt.clf(); plt.close()

def plotpotentials(pKa):
    R   = 8.3145 * 10**-3 # "kJ * mol⁻1 * K^-1"
    T   = 300

    lambda_i = load.Col("lambda_dwp.dat", 1, 942, 2062)
    V_bias   = load.Col("lambda_dwp.dat", 0, 942, 2062)
    
    pH = universe.get('ph_pH')

    V_pH = []
    for i in lambda_i:
        V_pH.append(R * T * np.log(10) * (pKa - pH) * i)

    V_comb = []
    for i in range(0, len(lambda_i)):
        V_comb.append(V_bias[i] + V_pH[i])

    plt.plot(lambda_i, V_bias, color="b", linestyle='--', label="$V_{bias}}$")
    plt.plot(lambda_i, V_pH, color="b", linestyle = ':', label="$V_{pH}$")
    plt.plot(lambda_i, V_comb, color="b", label="$V_{combined}$")

    plt.xlabel(r"$\lambda$-coordinate")
    plt.ylabel(r"$V$ (kJ/mol)")
    plt.grid()
    plt.legend()
    plt.show()

def plotforces(pKa):
    R   = 8.3145 * 10**-3 # "kJ * mol⁻1 * K^-1"
    T   = 300

    lambda_i = load.Col("lambda_dwp.dat", 1, 942, 2062)
    V_bias   = load.Col("lambda_dwp.dat", 0, 942, 2062)
    
    pH = universe.get('ph_pH')

    V_pH = []
    for i in lambda_i:
        V_pH.append(R * T * np.log(10) * (pKa - pH) * i)

    V_bias = np.gradient(V_bias)    # Take derivatives.
    V_pH   = np.gradient(V_pH)
    V_comb = [V_bias[i] + V_pH[i] for i in range(0, len(lambda_i))]

    plt.plot(lambda_i, V_bias, color="b", linestyle='--', label="$F_{bias}}$")
    plt.plot(lambda_i, V_pH, color="b", linestyle = ':', label="$F_{pH}$")
    plt.plot(lambda_i, V_comb, color="b", label="$F_{combined}$")

    plt.ylim(-0.1, 0.1)
    plt.xlabel(r"$\lambda$-coordinate")
    plt.ylabel("Force")
    plt.grid()
    plt.legend()
    plt.show()

def plothistogram(fname, bins=200):
    from scipy.stats import gaussian_kde
    
    data = load.Col(fname, 2)
    
    plt.hist(data, density=True, bins=bins)

    density = gaussian_kde(data)
    xs = np.linspace(-0.1, 1.1, bins)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    plt.plot(xs, density(xs), label="10 ns test")    

    plt.xlabel(r"$\lambda$-coordinate")
    # plt.axis([-0.1, 1.1, 0, 2.5])
    plt.xlim(-0.1, 1.1)
    plt.show()

def fitCalibration(order=5, compare=[]):
    # Get relevant stuff from universe.
    # Note: these data-members are only created when calibrate.py is ran.
    dVdlInitList = universe.get('ph_dvdl_initList')
    dVdlMeanList = universe.get('ph_dvdl_meanList')
    dVdlStdList  = universe.get('ph_dvdl_stdList')

    # Compute dV/dl coefficients.
    coeffs = np.polyfit(dVdlInitList, dVdlMeanList, order)[::-1]

    # Update user with the coefficients.
    print(coeffs)

    # Plot the computed values.
    plt.scatter(dVdlInitList, dVdlMeanList, label="mean dV/dl")
    plt.errorbar(dVdlInitList, dVdlMeanList, xerr=0, yerr=dVdlStdList, fmt='o', capsize=3, color='#1f77b4')

    # Our fit
    fit = []
    for i in dVdlInitList:
        value = 0
        for j in range(0, order + 1):
            value += coeffs[j] * i**j
        fit.append(value)
    plt.plot(dVdlInitList, fit, label="fit")

    # Comparison
    if len(compare) != 0:
        fit = []
        for i in dVdlInitList:
            value = 0
            for j in range(0, len(compare)):
                value += compare[j] * i**j
            fit.append(value)
        plt.plot(dVdlInitList, fit, label="compare")

    plt.title("Calibration for {}.pdb".format(universe.get('d_pdbName')))
    plt.ylabel(r"dV/d$\lambda$")
    plt.xlabel(r"$\lambda$-coordinate")
    plt.legend()
    plt.grid()
    plt.show()

def compareLambdaFiles(namelist):
    # If you accidentally put a string instead of a list, fix it.
    if (type(namelist) == type("")):
        namelist = [namelist]

    # Define (size of) main figure.
    fig = plt.figure(figsize=(24, 10))

    # Define sub-plots.
    plt1 = fig.add_subplot(2, 4, 1)
    plt2 = fig.add_subplot(2, 4, 2)
    plt3 = fig.add_subplot(2, 4, 3)
    plt4 = fig.add_subplot(2, 4, 4)
    plt5 = fig.add_subplot(2, 4, 5)
    plt6 = fig.add_subplot(2, 4, 6)
    plt7 = fig.add_subplot(2, 4, 7)
    plt8 = fig.add_subplot(2, 4, 8)

    # Get the data and plot.
    for name in namelist:
        time        = load.Col(name, 1)
        lambda_x    = load.Col(name, 2)
        lambda_dvdl = load.Col(name, 3)
        lambda_temp = load.Col(name, 4)
        lambda_vel  = load.Col(name, 5)
        F_coulomb   = load.Col(name, 6)
        F_corr      = load.Col(name, 7)
        F_bias      = load.Col(name, 8)
        F_ph        = load.Col(name, 9)

        plt1.plot(time, lambda_x, linewidth=0.5, label="deprotonation = {:.2f}".format(titrate(name)))
        plt2.plot(time, lambda_temp, linewidth=0.5, label="mean = {:.1f} (K)".format(sum(lambda_temp)/len(lambda_temp)))
        plt3.hist(lambda_vel, density=True)
        plt4.scatter(lambda_x, lambda_dvdl, s=5)
        plt5.scatter(lambda_x, F_coulomb, s=5)
        plt6.scatter(lambda_x, F_corr, s=5)
        plt7.scatter(lambda_x, F_bias, s=5)
        plt8.scatter(lambda_x, F_ph, s=5)

    plt1.set_title("$\lambda$-coordinate vs time")
    plt1.set_xlabel("Time (ps)")
    plt1.set_ylabel("$\lambda$-coordinate")
    plt1.set_ylim(-0.1, 1.1)
    plt1.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    plt1.legend()

    plt2.set_title("$\lambda$-temperature vs time")
    plt2.set_xlabel("Time (ps)")
    plt2.set_ylabel("$\lambda$-temperature (K)")
    plt2.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    plt2.legend()

    plt3.set_title("$\lambda$-velocity distribution")
    plt3.set_xlabel("$\lambda$-velocity (km/s)")

    plt4.set_title("Force (dV/dl) on $\lambda$-particle")
    plt4.set_xlabel("$\lambda$-coordinate")
    plt4.set_ylabel("dV/dl")
    plt4.set_xlim(-0.1, 1.1)

    plt5.set_title("Coulomb-force on $\lambda$-particle")
    plt5.set_xlabel("$\lambda$-coordinate")
    plt5.set_ylabel("$F_{Coulomb}$")
    plt5.set_xlim(-0.1, 1.1)

    plt6.set_title("Reference-force on $\lambda$-particle")
    plt6.set_xlabel("$\lambda$-coordinate")
    plt6.set_ylabel("$F_{corr}$")
    plt6.set_xlim(-0.1, 1.1)

    plt7.set_title("Bias-force on $\lambda$-particle")
    plt7.set_xlabel("$\lambda$-coordinate")
    plt7.set_ylabel("$F_{bias}$")
    plt7.axis([-0.1, 1.1, -200, 200])

    plt8.set_title("pH-force on $\lambda$-particle")
    plt8.set_xlabel("$\lambda$-coordinate")
    plt8.set_ylabel("$F_{pH}$")
    plt8.set_xlim(-0.1, 1.1)

    # Stuff we do in all subplots we can do in a loop:
    for plot in [plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8]:
        plot.grid()

    fig.legend(namelist, loc="upper center")
    # plt.tight_layout()
    plt.show()
