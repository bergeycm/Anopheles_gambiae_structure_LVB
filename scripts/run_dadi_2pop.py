#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Script to run demographic inference with dadi
# ----------------------------------------------------------------------------------------

# Usage example:
# python scripts/run_dadi_2pop.py \
#    -i data/dadi/dadi.chr3.mainlandisland.data \
#    -p MAINLAND \
#    -q ISLAND \
#    -s 50 \
#    -m 20 > \
#    results/dadi/dadi.island.2pop.MAINLAND.ISLAND.out

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

import demographic_model_2pop

from dadi import Demographics2D

do_bootsrapping = 1

def main(argv):
    infile = ''
    pop1 = ''
    pop2 = ''
    samp = 10
    maxiter = 10
    usage = 'run_dadi.py -i <inputfile> -p <pop1> -q <pop2> -s <sample> -m <maxiter>'
    try:
        opts, args = getopt.getopt(argv,"hi:p:q:s:m:",["in=","pop1=","pop2=","sample=","maxiter="])
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
        elif opt in ("-s", "--sample"):
            samp = int(arg)
        elif opt in ("-m", "--maxiter"):
            maxiter = int(arg)
    sys.stderr.write('Input file is ' + infile + '\n')
    sys.stderr.write('Population 1 is ' + pop1 + '\n')
    sys.stderr.write('Population 2 is ' + pop2 + '\n')
    sys.stderr.write('Sample is ' + str(samp) + '\n')
    sys.stderr.write('Max. iterations for optimization is ' + str(maxiter) + '\n')

    if (infile == "" or pop1 == "" or pop2 == ""):
        sys.stderr.write(usage + '\n')
        sys.exit(2)

    # Parse the data file to generate the data dictionary
    dd = dadi.Misc.make_data_dict(infile)

    # Extract the spectrum for chosen population from data dictionary,
    # projecting down to N samples per population.
    fs = dadi.Spectrum.from_data_dict(dd, [pop1, pop2], [samp,samp], polarized=True)

    # Plot the fs.
    fs_fig = pylab.figure(1)
    fs_fig.clear()
    dadi.Plotting.plot_single_2d_sfs(sfs=fs, pop_ids=[pop1, pop2])
    fs_fig.savefig(infile.replace(".data", "." + pop1 + "-" + pop2 + ".2d_fs.pdf"))

    ns = fs.sample_sizes

    # Grid point settings used for extrapolation.
    # Set to number of chromosomes and then increase by 10
    chr_count = 2 * samp
    pts_l = [chr_count, chr_count + 10, chr_count + 20]

    func = demographic_model_2pop.prior_onegrow_mig
    # prior_onegrow_mig:
    #     https://github.com/paulirish/dadi/blob/master/dadi/Demographics2D.py
    # Model with growth, split, bottleneck in pop2 , exp recovery, migration
    #     nu1F: The ancestral population size after growth. (Its initial size is
    #           defined to be 1.)
    #     nu2B: The bottleneck size for pop2
    #     nu2F: The final size for pop2
    #     m: The scaled migration rate
    #     Tp: The scaled time between ancestral population growth and the split.
    #     T: The time between the split and present

    #             nu1F  nu2B  nu2F     m   Tp    T
    params = array([ 5,    5,    5,  0.1, 0.1, 0.1])
    upper_bound = [ 75,   75,   75,   20,   5,   5]
    lower_bound = [0.1,  0.1,  0.1,    0,   0,   0]

    # func = Demographics2D.split_mig
    # # split_mig:
    # # https://github.com/paulirish/dadi/blob/master/dadi/Demographics2D.py
    # # Split into two populations of specifed size, with migration.
    # #     nu1: Size of population 1 after split.
    # #     nu2: Size of population 2 after split.
    # #     T: Time in the past of split (in units of 2*Na generations)
    # #     m: Migration rate between populations (2*Na*m)

    # #              nu1   nu2     T,    m
    # params = array([ 5,    5,  0.1,  0.1])
    # upper_bound = [100,  100,    3,   30]
    # lower_bound = [  0,    0,    0,    0]

    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    # Calculate the model AFS.
    model = func_ex(params, ns, pts_l)
    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, fs)
    sys.stderr.write('Model log-likelihood: ' + str(ll_model))
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
    sys.stderr.write('Optimized parameters ' + str(repr(popt)))

    # Fit model
    model = func_ex(popt, ns, pts_l)

    # Find optimal value of theta given the model
    theta_opt = dadi.Inference.optimal_sfs_scaling(model, fs)

    # Compute log composite likelihood
    ll_opt = dadi.Inference.ll_multinom(model, fs)
    sys.stderr.write('Optimized log-likelihood: ' + str(ll_opt))

    # Plot a comparison of the resulting fs with the data.
    fs_fig2 = pylab.figure(1)
    fs_fig2.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=3,
                                        pop_ids =(pop1,pop2))
    fs_fig2.savefig(infile.replace(".data", "." + pop1 + "-" + pop2 + ".cmp.2d_fs.pdf"))

    sys.stderr.write("\n")

    # ------------------------------------------------------------------------------------

    # Calculate AIC **assuming** SNPs are unlinked
    AIC = 2 * len(popt) - 2 * ll_opt

    # ------------------------------------------------------------------------------------

    # Estimate parameter uncertainties using the Godambe Information Matrix, to
    # account for linkage in the data.

    if (do_bootsrapping):
        mode = os.path.splitext(os.path.splitext(infile)[0])[1].replace(".", "")

        total_bs_iter = 500
        if mode == "ag1000g":
            total_bs_iter = 250

        # Load in data dictionaries from bootstrapping
        all_dd = [dadi.Misc.make_data_dict('data/dadi_bootstraps/bs_{0}.{1}.data'.format(ii, mode))
                for ii in range(total_bs_iter)]

        # Extract the spectra for all the bootstrapped data dictionaries
        all_boot = [dadi.Spectrum.from_data_dict(all_dd[ii], [pop1, pop2], [samp,samp], polarized=True)
                for ii in range(total_bs_iter)]

        # Estimate standard deviations of each parameter (value of theta is last)
        uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot,
                                popt, fs, multinom=True)


        # Also estimate uncertainties with the Fisher Information Matrix,
        # which doesn't account for linkage in the data and thus underestimates uncertainty.
        uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, fs,
                                multinom=True)

    # ------------------------------------------------------------------------------------

    # Print header

    sys.stdout.write("input_file\tpop1\tpop2\tsample_size\t")

    # Likelihood and AIC
    sys.stdout.write("ll\tAIC\t")

    sys.stdout.write("theta\t")

    # Best fit parameter values
    sys.stdout.write("nu1F\tnu2B\tnu2F\tm\tTp\tT\t")
    #sys.stdout.write("nu1\tnu2\tT\tm\t")

    # Uncertainties (only GIM output)
    if (do_bootsrapping):
        for param in ("nu1F", "nu2B", "nu2F", "m", "Tp", "T"):
            sys.stdout.write("GIM_" + param + "\t")
        sys.stdout.write("GIM_theta\t")

    sys.stdout.write("max_iter\n")

    # ------------------------------------------------------------------------------------

    # Output basic stats and likelihood info
    print "%s\t%s\t%s\t%d\t" % (infile.replace("data/dadi_polarized/", ""),
                                pop1, pop2, samp),

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

    # Print maximum iterations for this run
    print maxiter

if __name__ == "__main__":
    main(sys.argv[1:])
