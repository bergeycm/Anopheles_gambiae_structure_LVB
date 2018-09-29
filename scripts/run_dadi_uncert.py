#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Script to run demographic inference with dadi
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- Installation of dadi:
# module load python/2.7.8
# export PYTHONPATH=/storage/home/cxb585/local_python
# export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi:/usr/global/python/2.7.8/lib/python2.7/site-packages
# pip install --no-cache-dir --ignore-installed --target=$HOME/local_python/ --upgrade numpy scipy
# cd $HOME/bin
# git clone https://bitbucket.org/gutenkunstlab/dadi.git; cd dadi
# module load gcc
# python setup.py build_ext --inplace
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- Installation of dadi on new cluster
# module load python/2.7.14-anaconda5.0.1
# python -m pip install --target=$HOME/local_python/ --upgrade --force-reinstall --no-cache-dir numpy scipy matplotlib
# export PYTHONPATH=$PYTHONPATH:$HOME/local_python/
# # /opt/aci/sw/python/2.7/data:/opt/aci/sw/python/2.7/lib:/storage/home/cxb585/local_python:/storage/home/cxb585/local_python/lib64/python2.7/site-packages:/storage/home/cxb585/local_python/
# cd $HOME/bin
# git clone https://bitbucket.org/gutenkunstlab/dadi.git; cd dadi
# module load gcc
# python setup.py build_ext --inplace
# export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi
# ----------------------------------------------------------------------------------------

# Before run:
# module load python/2.7.8
# export PYTHONPATH=/storage/home/cxb585/local_python
# export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi:/usr/global/python/2.7.8/lib/python2.7/site-packages

# Usage example:
# python scripts/run_dadi.py \
#    -i data/dadi/dadi.mainlandisland.thin.data \
#    -p ISLAND > \
#    results/dadi/dadi.mainlandisland.thin.ISLAND.out

# matplotlib becomes upset when it cannot connect to an X server for GTK
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import sys, getopt
from numpy import array
import numpy as np
import dadi
import pylab
import os
import re

# Import custom model for this problem (with uncertainty)
sys.path.append(os.path.abspath("scripts"))
sys.path.append(os.path.abspath("/gpfs/scratch/cxb585/ssese-wgs/ssese_wgs"))
#from scripts import *

import Demographics1D_mod

do_boostrapping = 1

def main(argv):
    infile = ''
    pop = ''
    samp = 10
    maxiter = 10
    fold = False
    usage = 'run_dadi.py -i <inputfile> -p <population> -s <sample> -m <maxiter> -f <fold> -e <ending>'
    try:
        opts, args = getopt.getopt(argv,"hi:p:s:m:f:e:",["in=","pop=","sample=","maxiter=","fold=","ending="])
    except getopt.GetoptError:
        sys.stderr.write(usage + '\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.stderr.write(usage + '\n')
            sys.exit()
        elif opt in ("-i", "--in"):
            infile = arg
        elif opt in ("-p", "--pop"):
            pop = arg
        elif opt in ("-s", "--sample"):
            samp = int(arg)
        elif opt in ("-m", "--maxiter"):
            maxiter = int(arg)
        elif opt in ("-f", "--fold"):
            fold = arg == 'True'
        elif opt in ("-e", "--ending"):
            ending = arg
    sys.stderr.write('Input file is ' + infile + '\n')
    sys.stderr.write('Population is '+ pop + '\n')
    sys.stderr.write('Sample is '+ str(samp) + '\n')
    sys.stderr.write('Max. iterations for optimization is ' + str(maxiter) + '\n')
    sys.stderr.write('Do folding is ' + str(fold) + '\n')
    sys.stderr.write('Suffix is ' + ending + '\n')

    if (infile == "" or pop == ""):
        sys.stderr.write(usage + '\n')
        sys.exit(2)

    if (fold == True):
        fold_suffix = ".folded"
    else:
        fold_suffix = ""

    # ------------------------------------------------------------------------------------

    # Parse the data file to generate the data dictionary
    dd = dadi.Misc.make_data_dict(infile)

    # Extract the spectrum for chosen population from data dictionary,
    # projecting down to N samples per population.
    fs = dadi.Spectrum.from_data_dict(dd, [pop], [samp], polarized=(not fold))

    # Mask singletons and doubletons
    fs.mask[1] = True
    fs.mask[2] = True
    fs.mask[-1] = True
    fs.mask[-2] = True

    # While we're here, also compute the SFS without projecting down
    # We use it for the stair
    pop_max = 0
    for i in range(1, 100):
        if (sum(dd[dd.keys()[1]]['calls'][pop]) > pop_max):
            pop_max = sum(dd[dd.keys()[1]]['calls'][pop])

    fs_full = dadi.Spectrum.from_data_dict(dd, [pop], [pop_max], polarized=(not fold))

    fs_full_str = re.sub("-- *", "", re.sub("[\n\[\]]", "", str(fs_full)))

    sfs_file = open(infile.replace(".data", "." + pop + fold_suffix + ending + ".sfs.txt"), 'w')
    sfs_file.write(fs_full_str + "\n")
    sfs_file.close()

    # Plot the fs.
    fs_fig = pylab.figure(1)
    fs_fig.clear()
    dadi.Plotting.plot_1d_fs(fs)
    fs_fig.savefig(infile.replace(".data", "." + pop + fold_suffix + ending + ".1d_fs.pdf"))

    # ====================================================================================

    ns = fs.sample_sizes

    # These are the grid point settings will use for extrapolation.
    # Set to number of chromosomes and then increase by 10
    chr_count = 2 * samp
    pts_l = [chr_count, chr_count + 10, chr_count + 20]

    # ------------------------------------------------------------------------------------
    # --- Compare the four scenarios
    # ------------------------------------------------------------------------------------

    #    snm
    #        [no params]
    #    two_epoch
    #        nu, T
    #    growth
    #        nu, T
    #    bottlegrowth
    #        nuB, nuF, T
    #    three_epoch
    #        nuB, nuF, TB, TF

    models = ['snm_uncert', 'two_epoch_uncert', 'growth_uncert',
              'bottle_uncert', 'three_epoch_uncert']

    params = {}
    params['snm_uncert']         = [                         'p_misid']
    params['two_epoch_uncert']   = ['nu',         'T'       ,'p_misid']
    params['growth_uncert']      = ['nu',         'T'       ,'p_misid']
    params['bottle_uncert']      = ['nuB', 'nuF', 'T'       ,'p_misid']
    params['three_epoch_uncert'] = ['nuB', 'nuF', 'TB', 'TF','p_misid']

    func = {}
    func['snm_uncert']         = Demographics1D_mod.snm_uncert
    func['two_epoch_uncert']   = Demographics1D_mod.two_epoch_uncert
    func['growth_uncert']      = Demographics1D_mod.growth_uncert
    func['bottle_uncert']      = Demographics1D_mod.bottlegrowth_uncert
    func['three_epoch_uncert'] = Demographics1D_mod.three_epoch_uncert

    # Make the extrapolating versions of model functions
    func_ex = {}
    for mod in models:
        func_ex[mod] = dadi.Numerics.make_extrap_log_func(func[mod])

    # --- Optimize

    bound = {}
    bound['upper'] = {}
    bound['lower'] = {}

    bound['upper']['snm_uncert']         = [                                1]
    bound['lower']['snm_uncert']         = [                             1e-8]

    bound['upper']['two_epoch_uncert']   = [10000,          0.1,            1]
    bound['lower']['two_epoch_uncert']   = [ 1e-2,         1e-8,         1e-8]

    bound['upper']['growth_uncert']      = [10000,          0.1,            1]
    bound['lower']['growth_uncert']      = [ 1e-2,         1e-8,         1e-8]

    bound['upper']['bottle_uncert']      = [10000, 10000,   0.1,            1]
    bound['lower']['bottle_uncert']      = [ 1e-2,  1e-2,  1e-8,         1e-8]

    bound['upper']['three_epoch_uncert'] = [10000, 10000,   0.1,   0.1,     1]
    bound['lower']['three_epoch_uncert'] = [ 1e-2,  1e-2,  1e-8,  1e-8,  1e-8]

    # Arbitrary initial guess
    p0 = {}
    p0['snm_uncert']         = [                  0.05]
    p0['two_epoch_uncert']   = [1,    0.01,       0.05]
    p0['growth_uncert']      = [1,    0.01,       0.05]
    p0['bottle_uncert']      = [1, 1, 0.01,       0.05]
    p0['three_epoch_uncert'] = [1, 1, 0.01, 0.01, 0.05]

    # Perturb parameters before optimization
    for mod in models:
        p0[mod] = dadi.Misc.perturb_params(p0[mod], fold=1.5,
                                  upper_bound=bound['upper'][mod],
                                  lower_bound=bound['lower'][mod])

    # Do the optimization.
    # Value for maxiter should be greater than 10
    popt = {}
    for mod in models:
        popt[mod] = dadi.Inference.optimize_log(p0[mod],
                                       fs, func_ex[mod], pts_l,
                                       lower_bound=bound['upper'][mod],
                                       upper_bound=bound['lower'][mod],
                                       verbose=len(p0[mod]), maxiter=maxiter)

    # Fit model
    model = {}
    for mod in models:
        model[mod] = func_ex[mod] (popt[mod], ns, pts_l)

    # Find optimal value of theta given the model
    theta = {}
    for mod in models:
        theta[mod] = dadi.Inference.optimal_sfs_scaling(model[mod], fs)

    # Compute log composite likelihood
    ll = {}
    max_ll = -999999999999
    max_ll_model = ''
    for mod in models:
        ll[mod] = dadi.Inference.ll_multinom(model[mod], fs)
        sys.stdout.write("log L of " + mod + ": " + str(ll[mod]) + "\n")
        if (ll[mod] > max_ll):
            max_ll = ll[mod]
            max_ll_model = mod
    sys.stdout.write("\tBest model by log L: " + max_ll_model + "\n")

    # ------------------------------------------------------------------------------------

    # Calculate AIC **assuming** SNPs are unlinked
    AIC = {}
    for mod in models:
        AIC[mod] = 2 * len(popt[mod]) - 2 * ll[mod]

    # ------------------------------------------------------------------------------------

    # Estimate parameter uncertainties using the Godambe Information Matrix, to
    # account for linkage in the data.

    #    matrix_file = infile.replace(".data", "." + pop + fold_suffix + ending +
    #                       ".GIM_uncert_matrix.txt")

    #    try:
    #        os.remove(matrix_file)
    #    except OSError:
    #        pass

    if (do_boostrapping):
        mode = os.path.splitext(os.path.splitext(infile)[0])[1].replace(".", "")

        total_bs_iter = 100
        if mode == "ag1000g":
            total_bs_iter = 100

        # Load in data dictionaries from bootstrapping
        all_dd = [dadi.Misc.make_data_dict('data/dadi_bootstraps/bs_{0}.{1}.data'.format(ii, mode))
                for ii in range(total_bs_iter)]

        # Extract the spectra for all the bootstrapped data dictionaries
        all_boot = [dadi.Spectrum.from_data_dict(all_dd[ii], [pop], [samp], polarized=(not fold))
                for ii in range(total_bs_iter)]

        # Estimate standard deviations of each parameter (value of theta is last)
        uncerts = {}
        for mod in models:
            uncerts[mod] = dadi.Godambe.GIM_uncert(func_ex[mod], pts_l, all_boot,
                                      popt[mod], fs, multinom=True, return_GIM=False)
            #   # Save uncertainty matrix
            #   f = open(matrix_file, "a")
            #   np.savetxt(f, uncerts[mod][1])
            #   f.close()

        # Also estimate uncertainties with the Fisher Information Matrix,
        # which doesn't account for linkage in the data and thus underestimates uncertainty.
        uncerts_fim = {}
        for mod in models:
            uncerts_fim[mod] = dadi.Godambe.FIM_uncert(func_ex[mod], pts_l, popt[mod], fs,
                                    multinom=True)

    # ------------------------------------------------------------------------------------

    # Print header

    sys.stdout.write("input_file\tgroup\tsample_size\t")

    # Likelihood and AIC
    for mod in models:
        sys.stdout.write("ll_" + mod + "\t")

    sys.stdout.write("max_ll_model\tmax_ll\t")

    for mod in models:
        sys.stdout.write("AIC_" + mod + "\t")

    for mod in models:
        sys.stdout.write("theta_" + mod + "\t")

    # Best fit parameter values
    for mod in models:
        for param in params[mod]:
            sys.stdout.write(mod + "_" + param + "\t")

    # Uncertainties (GIM and FIM output)
    if (do_boostrapping):
        for mod in models:
            for param in params[mod]:
                sys.stdout.write("GIM_" + mod + "_" + param + "\t")
            sys.stdout.write("GIM_" + mod + "_theta\t")
        for mod in models:
            for param in params[mod]:
                sys.stdout.write("FIM_" + mod + "_" + param + "\t")
            sys.stdout.write("FIM_" + mod + "_theta\t")

    sys.stdout.write("max_iter\n")

    # ------------------------------------------------------------------------------------

    # Output basic stats and likelihood info
    print "%s\t%s\t%d\t" % (infile.replace("data/dadi_polarized/", ""),
                                pop, samp),

    for mod in models:
        print "%.6f\t" % (ll[mod]),

    print "%s\t%.6f\t" % (max_ll_model, max_ll),

    for mod in models:
        print "%.6f\t" % (AIC[mod]),

    # Print optimal value of theta
    for mod in models:
        print "%.6f\t" % (theta[mod]),

    # Print optimal parameters
    for mod in models:
        for p in popt[mod]:
            print "%.6f\t" % p,

    # Print GIM uncertainties
    if (do_boostrapping):
        for mod in models:
            for gim in uncerts[mod]:
                print "%.6f\t" % gim,
        for mod in models:
            for fim in uncerts_fim[mod]:
                print "%.6f\t" % fim,

    # Print maximum iterations for this run
    print maxiter

    # ------------------------------------------------------------------------------------

    # # Fix names
    # pop_dic = {'KAZZI':'KAAZI', 'MITYANA':'WAMALA',
    #            'BUGALAI':'BUGALA (I)', 'BUGALAM':'BUGALA (M)'}
    # pop_id_fix = pop_dic.get(pop) if pop in pop_dic else pop

    # Visually compare data and models
    for mod in models:
        cmp_fig = pylab.figure(1)
        cmp_fig.clear()
        dadi.Plotting.plot_1d_comp_multinom(model[mod], fs)
        cmp_fig.savefig(infile.replace(".data", "." + pop + "." + mod + fold_suffix + ending + ".pdf"))

    # Save best model image again
    cmp_fig = pylab.figure(1)
    cmp_fig.clear()
    dadi.Plotting.plot_1d_comp_multinom(model[max_ll_model], fs)
    cmp_fig.savefig(infile.replace(".data", "." + pop + ".bestmodel" + fold_suffix + ending + ".pdf"))

if __name__ == "__main__":
    main(sys.argv[1:])
