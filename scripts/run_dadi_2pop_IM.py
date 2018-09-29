#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Script to run demographic inference with dadi
# ----------------------------------------------------------------------------------------

# Usage example:
# python scripts/run_dadi_2pop_IM.py \
#    -i data/dadi/dadi.chr3.island.data \
#    -p BANDA \
#    -q BUGALAML \
#    -s 10 \
#    -m 100 > \
#    results/dadi/dadi.island.2pop.BANDA.BUGALAML.out

# matplotlib becomes upset when it cannot connect to an X server for GTK
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import sys, getopt
import numpy
from numpy import array
import dadi
import pylab
import os

# Import custom model for this problem
sys.path.append(os.path.abspath("scripts"))
sys.path.append(os.path.abspath("/gpfs/scratch/cxb585/ssese-wgs/ssese_wgs"))
from scripts import *

#import demographic_model_2pop

#from dadi import Demographics2D

import Demographics2D_mod

do_bootsrapping = 1

def main(argv):
    infile = ''
    pop1 = ''
    pop2 = ''
    samp1 = 10
    samp2 = 10
    maxiter = 10
    usage = 'run_dadi.py -i <inputfile> -p <pop1> -q <pop2> -s <sample_1> -z <sample_2> -m <maxiter> -e ending'
    try:
        opts, args = getopt.getopt(argv,"hi:p:q:s:z:m:e:",
            ["in=","pop1=","pop2=","sample_1=","sample_2=","maxiter=","ending="])
    except getopt.GetoptError:
        sys.stderr.write(usage + '\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.stderr.write(usage + '\n')
            sys.exit()
        elif opt in ("-i", "--in"):
            infile = arg
        elif opt in ("-p", "--pop1"):
            pop1 = arg
        elif opt in ("-q", "--pop2"):
            pop2 = arg
        elif opt in ("-s", "--sample_1"):
            samp1 = int(arg)
        elif opt in ("-z", "--sample_2"):
            samp2 = int(arg)
        elif opt in ("-m", "--maxiter"):
            maxiter = int(arg)
        elif opt in ("-e", "--ending"):
            ending = arg
    sys.stderr.write('Input file is ' + infile + '\n')
    sys.stderr.write('Population 1 is ' + pop1 + '\n')
    sys.stderr.write('Population 2 is ' + pop2 + '\n')
    sys.stderr.write('Sample for pop 1 is ' + str(samp1) + '\n')
    sys.stderr.write('Sample for pop 2 is ' + str(samp2) + '\n')
    sys.stderr.write('Max. iterations for optimization is ' + str(maxiter) + '\n')
    sys.stderr.write('Suffix is ' + ending + '\n')

    if (infile == "" or pop1 == "" or pop2 == ""):
        sys.stderr.write(usage + '\n')
        sys.exit(2)

    # For now, let's just use the points of the smallest population
    if (samp1 > samp2):
        samp1 = samp2
    else:
        samp2 = samp1

    # Parse the data file to generate the data dictionary
    dd = dadi.Misc.make_data_dict(infile)

    # Extract the spectrum for chosen population from data dictionary,
    # projecting down to N samples per population.
    fs = dadi.Spectrum.from_data_dict(dd, [pop1, pop2], [samp1,samp2], polarized=True)

    # Mask private singletons and doubletons
    fs.mask[1:3,0] = True
    fs.mask[0,1:3] = True

    # Plot the fs.
    fs_fig = pylab.figure(1)
    fs_fig.clear()
    dadi.Plotting.plot_single_2d_sfs(sfs=fs, vmin=3, pop_ids=[pop1, pop2])
    fs_fig.savefig(infile.replace(".data", "." + pop1 + "-" + pop2 + ending + ".2d_fs.pdf"))

    ns = fs.sample_sizes

    # Grid point settings used for extrapolation.
    # Set to number of chromosomes and then increase by 10
    chr_count = 2 * min(samp1, samp2)
    pts_l = [chr_count, chr_count + 10, chr_count + 20]

    # ------------------------------------------------------------------------------------

    # --- Model without migration

    func = Demographics2D_mod.IM_actually_no_mig

    #####                 s    nu1    nu2     T  p_misid
    ####params = array([0.5,    50,    50,  0.1,    0.05])
    ####upper_bound = [   1, 10000, 10000,    3,       1]
    ####lower_bound = [   0,     0,     0,    0,       0]


    #                 s    nu1    nu2     T  p_misid
    params = array([0.5,     1,     1, 0.01,    0.05])
    upper_bound = [   1, 10000, 10000,  0.1,       1]
    lower_bound = [1e-8,  1e-2,  1e-2, 1e-8,    1e-8]

    #func_ex_nomig = dadi.Numerics.make_extrap_log_func(func)
    func_ex_nomig = dadi.Numerics.make_extrap_func(func) # Switching from log to linear
    model_nomig = func_ex_nomig(params, ns, pts_l)
    ll_model_nomig = dadi.Inference.ll_multinom(model_nomig, fs)
    theta_nomig = dadi.Inference.optimal_sfs_scaling(model_nomig, fs)
    p0_nomig = dadi.Misc.perturb_params(params, fold=1.5, upper_bound=upper_bound)
    popt_nomig = dadi.Inference.optimize_log(p0_nomig, fs, func_ex_nomig, pts_l,
                                       lower_bound=lower_bound,
                                       upper_bound=upper_bound,
                                       verbose=0,
                                       maxiter=maxiter)
    model_nomig = func_ex_nomig(popt_nomig, ns, pts_l)
    theta_opt_nomig = dadi.Inference.optimal_sfs_scaling(model_nomig, fs)
    ll_opt_nomig = dadi.Inference.ll_multinom(model_nomig, fs)
    sys.stderr.write('No-migration model optimized log-likelihood: ' +
        str(ll_opt_nomig) + "\n")

    # Quick test:
    L = 88671774
    mu = 3.5e-9
    g = 1.0 / 11
    test_theta = theta_opt_nomig
    T = popt_nomig[3]
    Na     = test_theta      / (4 * mu * L)
    T_real     = 2 * Na * T     * g

    sys.stderr.write("For no-migration model, initial pop of " + str(round(Na)) +
        " splits " + str(round(T_real)) + " years ago and " +
        "the daughter pops grow by factors of " +
        str(round(popt_nomig[1])) + " and " + str(round(popt_nomig[2])) + ".\n")

    # ------------------------------------------------------------------------------------

    # Calculate AIC **assuming** SNPs are unlinked
    AIC_nomig = 2 * len(popt_nomig) - 2 * ll_opt_nomig
    sys.stderr.write('No-migration model AIC: ' + str(round(AIC_nomig)) + "\n")

    # ----------------------------------------------------------------------------------------

    # --- Model with migration

    func = Demographics2D_mod.IM
    # IM:
    # https://github.com/paulirish/dadi/blob/master/dadi/Demographics2D.py
    # Isolation-with-migration model with exponential pop growth.
    #     s:    *Fraction* of pop 1 after split. (Pop 2 has size 1-s.)
    #     nu1:  Final size of pop 1.
    #     nu2:  Final size of pop 2.
    #     T:    Time in the past of split (in units of 2*Na generations)
    #     m12:  Migration from pop 2 to pop 1 (2*Na*m12)
    #     m21:  Migration from pop 1 to pop 2

    #####                 s    nu1    nu2     T   m12   m21  p_misid
    ####params = array([0.5,    50,    50,  0.1,    3,    3,    0.05])
    ####upper_bound = [   1, 10000, 10000,    3,  300,  300,       1]
    ####lower_bound = [   0,     0,     0,    0,    0,    0,       0]

    #                 s    nu1    nu2     T   m12   m21  p_misid
    params = array([0.5,     1,     1, 0.01,    3,    3,    0.05])
    upper_bound = [   1, 10000, 10000,  0.1,   10,   10,       1]
    lower_bound = [1e-8,  1e-2,  1e-2, 1e-8, 1e-8, 1e-8,    1e-8]

    # Make the extrapolating version of our demographic model function.
    #func_ex = dadi.Numerics.make_extrap_log_func(func)
    func_ex = dadi.Numerics.make_extrap_func(func) # Switching from log to linear
    # Calculate the model AFS.
    model = func_ex(params, ns, pts_l)
    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, fs)
    sys.stderr.write('Model log-likelihood: ' + str(ll_model) + "\n")
    # The optimal value of theta given the model.
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)

    # Perturb our parameter array before optimization. This does so by taking each
    # parameter a up to a factor of two up or down.
    p0 = dadi.Misc.perturb_params(params, fold=1.5, upper_bound=upper_bound)
    # Do the optimization. By default we assume that theta is a free parameter,
    # since it's trivial to find given the other parameters. If you want to fix
    # theta, add a multinom=False to the call.
    # The maxiter argument restricts how long the optimizer will run.
    popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                       lower_bound=lower_bound,
                                       upper_bound=upper_bound,
                                       verbose=0,
                                       maxiter=maxiter)
    sys.stderr.write('Optimized parameters ' + str(repr(popt)) + "\n")

    # Fit model
    model = func_ex(popt, ns, pts_l)

    # Find optimal value of theta given the model
    theta_opt = dadi.Inference.optimal_sfs_scaling(model, fs)

    # Compute log composite likelihood
    ll_opt = dadi.Inference.ll_multinom(model, fs)
    sys.stderr.write('Optimized log-likelihood: ' + str(ll_opt) + "\n")

    # Quick test:
    test_theta = theta_opt
    T = popt[3]
    Na     = test_theta      / (4 * mu * L)
    T_real     = 2 * Na * T     * g
    sys.stderr.write("For migration model, initial pop of " + str(round(Na)) +
        " splits " + str(round(T_real)) + " years ago and " +
        "the daughter pops grow by factors of " +
        str(round(popt[1])) + " and " + str(round(popt[2])) + ".\n")

    # ------------------------------------------------------------------------------------

    # Calculate AIC **assuming** SNPs are unlinked
    AIC = 2 * len(popt) - 2 * ll_opt
    sys.stderr.write('With-migration model AIC: ' + str(AIC) + "\n")

    # ------------------------------------------------------------------------------------

    nomig_better = False

    if (do_bootsrapping):
        mode = os.path.splitext(os.path.splitext(infile)[0])[1].replace(".", "")

        total_bs_iter = 100
        if mode == "ag1000g":
            total_bs_iter = 100

        # Load in data dictionaries from bootstrapping
        all_dd = [dadi.Misc.make_data_dict('data/dadi_bootstraps/bs_{0}.{1}.data'.format(ii, mode))
                for ii in range(total_bs_iter)]

        # Extract the spectra for all the bootstrapped data dictionaries
        all_boot = [dadi.Spectrum.from_data_dict(all_dd[ii], [pop1, pop2], [samp1,samp2],
                polarized=True)
                for ii in range(total_bs_iter)]

        # --- Do LRT to compare models with and without migration.

        # "Since LRT evaluates the complex model using the best-fit parameters from the
        # simple model, we need to create list of parameters for the complex model
        # using the simple (no-mig) best-fit params.  Since evalution is done with more
        # complex model, need to insert zero migration value at corresponding migration
        # parameter index in complex model. And we need to tell the LRT adjust function
        # that the 3rd parameter (counting from 0) is the nested one."

        p_lrt = numpy.copy(popt)
        p_lrt[0:4] = popt_nomig[0:4]
        p_lrt[4:6] = 0
        p_lrt[6] = popt_nomig[4]

        adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, p_lrt, fs,
                                      nested_indices=[4,5], multinom=True)

        D_adj = adj*2*(ll_opt - ll_model_nomig)
        sys.stderr.write('Adjusted D statistic: {0:.4f}'.format(D_adj) + "\n")

        # Because this is test of a parameter on the boundary of parameter space
        # (m cannot be less than zero), our null distribution is an even proportion
        # of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
        # point percent function for a weighted sum of chi^2 dists.
        pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
        sys.stderr.write(
            'p-value for rejecting no-migration model: {0:.4f}'.format(pval) + "\n"
        )

        if (pval >= 0.05):
            nomig_better = True
    else:
        if AIC > AIC_nomig:
            nomig_better = True

    # If no-migration model is better, switch everything to that

    if nomig_better:
        sys.stderr.write("No migration model superior.\n")
        func_ex = func_ex_nomig
        model = model_nomig
        ll_model = ll_model_nomig
        theta = theta_nomig
        p0 = p0_nomig
        popt = popt_nomig
        theta_opt = theta_opt_nomig
        ll_opt = ll_opt_nomig
        param_list = ["s", "nu1", "nu2", "T", "pmisid"]
    else:
        sys.stderr.write("Model with migration superior.\n")
        param_list = ["s", "nu1", "nu2", "T", "m12", "m21", "pmisid"]

    # ------------------------------------------------------------------------------------

    # Plot a comparison of the resulting fs with the data.
    fs_fig2 = pylab.figure(1)
    fs_fig2.clear()

    # Fix names
    pop_dic = {'KAZZI':'KAAZI', 'MITYANA':'WAMALA',
               'BUGALAIS':'BUGALA (I)', 'BUGALAML':'BUGALA (M)'}
    pop_ids = [pop_dic.get(n, n) for n in [pop1, pop2]]

    sys.stderr.write("Fixed pop_ids:\t" + str(tuple(pop_ids)) + "\n")

    dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=3, resid_range=50,
                                        pop_ids = tuple(pop_ids))
    fs_fig2.savefig(infile.replace(".data", "." + pop1 + "-" + pop2 + ending + ".cmp.2d_fs.pdf"))

    sys.stderr.write("\n")

    # ------------------------------------------------------------------------------------

    # Estimate parameter uncertainties using the Godambe Information Matrix, to
    # account for linkage in the data.

    #    matrix_file = infile.replace(".data", "." + pop1 + "-" + pop2 + ending +
    #                       ".GIM_uncert_matrix.txt")

    #    try:
    #        os.remove(matrix_file)
    #    except OSError:
    #        pass

    if (do_bootsrapping):
        #mode = os.path.splitext(os.path.splitext(infile)[0])[1].replace(".", "")

        #total_bs_iter = 300
        #if mode == "ag1000g":
        #    total_bs_iter = 250

        ## Load in data dictionaries from bootstrapping
        #all_dd = [dadi.Misc.make_data_dict('data/dadi_bootstraps/bs_{0}.{1}.data'.format(ii, mode))
        #        for ii in range(total_bs_iter)]

        ## Extract the spectra for all the bootstrapped data dictionaries
        #all_boot = [dadi.Spectrum.from_data_dict(all_dd[ii], [pop1, pop2], [samp1,samp2],
        #        polarized=True)
        #        for ii in range(total_bs_iter)]

        # Estimate standard deviations of each parameter (value of theta is last)
        uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot,
                                popt, fs, multinom=True, return_GIM=False)

        #   # Save uncertainty matrix
        #   f = open(matrix_file, "w")
        #   numpy.savetxt(f, uncerts[1])
        #   f.close()

        # Also estimate uncertainties with the Fisher Information Matrix,
        # which doesn't account for linkage in the data and thus underestimates uncertainty.
        uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, fs,
                                multinom=True)

    # ------------------------------------------------------------------------------------

    # Print header

    sys.stdout.write("input_file\tpop1\tpop2\tsample_size1\tsample_size2\t")

    # Likelihood and AIC
    sys.stdout.write("ll\tAIC\t")

    sys.stdout.write("theta\t")

    # Best fit parameter values
    sys.stdout.write("\t".join(param_list) + "\t")

    # Uncertainties (GIM and FIM output)
    if (do_bootsrapping):
        for param in param_list:
            sys.stdout.write("GIM_" + param + "\t")
        sys.stdout.write("GIM_theta\t")
        for param in param_list:
            sys.stdout.write("FIM_" + param + "\t")
        sys.stdout.write("FIM__theta\t")

    sys.stdout.write("max_iter\n")

    # ------------------------------------------------------------------------------------

    # Output basic stats and likelihood info
    print "%s\t%s\t%s\t%d\t%d\t" % (infile.replace("data/dadi_polarized/", ""),
                                pop1, pop2, samp1, samp2),

    print "%.6f\t%.6f\t" % (ll_opt, AIC),

    # Print optimal value of theta
    print "%.6f\t" % (theta_opt),

    # Print optimal parameters
    for p in popt:
        print "%.6f\t" % p,

    # Print GIM uncertainties
    if (do_bootsrapping):
        for gim in uncerts:
            print "%.6f\t" % gim,
        for fim in uncerts_fim:
            print "%.6f\t" % fim,

    # Print maximum iterations for this run
    print maxiter

if __name__ == "__main__":
    main(sys.argv[1:])
