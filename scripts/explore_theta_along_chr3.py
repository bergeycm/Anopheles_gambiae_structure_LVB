#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Script to see how dadi-inferred Na changes across chromosome 3
# ----------------------------------------------------------------------------------------

# Before run:
# module load python/2.7.14-anaconda5.0.1
# export PYTHONPATH=/storage/home/cxb585/local_python
# export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/local_python//lib64/python2.7/site-packages
# export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi;

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

import Demographics1D_mod

def main(argv):
    infile = ''
    pop = ''
    samp = 10
    maxiter = 10
    fold = False
    usage = 'explore_theta_along_chr3.py -i <inputfile> -p <population> -s <sample> -m <maxiter> -f <fold> -e <ending>'
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

    # ====================================================================================

    ns = fs.sample_sizes

    # These are the grid point settings will use for extrapolation.
    # Set to number of chromosomes and then increase by 10
    chr_count = 2 * samp
    pts_l = [chr_count, chr_count + 10, chr_count + 20]

    mode = os.path.splitext(os.path.splitext(infile)[0])[1].replace(".", "")

    total_bs_iter = 100

    # Load in data dictionaries from bootstrapping
    all_dd = [dadi.Misc.make_data_dict('data/dadi_bootstraps/bs_{0}.{1}.data'.format(ii, mode))
            for ii in range(total_bs_iter)]

    # Extract the spectra for all the bootstrapped data dictionaries
    all_boot = [dadi.Spectrum.from_data_dict(all_dd[ii], [pop], [samp], polarized=(not fold))
            for ii in range(total_bs_iter)]

    # ====================================================================================

    model = 'three_epoch_uncert'
    params = ['nuB', 'nuF', 'TB', 'TF','p_misid']
    func = Demographics1D_mod.three_epoch_uncert
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    bound = {}
    bound['upper'] = [10000, 10000,     3,     3,     1]
    bound['lower'] = [    0,     0,     0,     0,     0]

    out = infile.replace("data/dadi_polarized", "results")
    out = out = out.replace(".data", "." + pop + ".theta_along_chr3.txt")

    theta_out_file = open(out, 'w')

    all_thetas = []

    for this_fs_idx in range(total_bs_iter):

        this_fs = all_boot[this_fs_idx]

        p0 = [50, 50, 0.1, 0.1, 0.05]

        p0 = dadi.Misc.perturb_params(p0, fold=1.5,
                                      upper_bound=bound['upper'],
                                      lower_bound=bound['lower'])

        popt = dadi.Inference.optimize_log(p0,
                                           this_fs, func_ex, pts_l,
                                           lower_bound=bound['upper'],
                                           upper_bound=bound['lower'],
                                           verbose=len(p0), maxiter=maxiter)

        model = func_ex (popt, ns, pts_l)

        theta = dadi.Inference.optimal_sfs_scaling(model, this_fs)

        this_boot = all_dd[this_fs_idx]
        boot_keys = this_boot.keys()

        all_thetas.append(theta)

        theta_out_file.write("%s\t%s\t%s\n" % (boot_keys[1], boot_keys[-1], theta))

    theta_out_file.close()

if __name__ == "__main__":
    main(sys.argv[1:])
