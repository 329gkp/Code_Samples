'''Code sample: software developed by the IceCube Neutrino Observatory is kept in private repos, so here I share a limited example of my work.
Part of my research included developing and utilizing SnowStorm, an efficient method for treating systematic uncertainties IceCube analyses.
For the full physics motivation and procedure behind SnowStorm, find our paper here: https://iopscience.iop.org/article/10.1088/1475-7516/2019/10/048
IceCube data ("events") are particles detected with an energy and direction ("zenith").  
This file is a modified sample of the program that extracts nuisance parameter gradients in energy and zenith from a Monte Carlo set of simulated events.'''

import tables
import numpy as np
import scipy as sp
from scipy.stats import norm as nm
import itertools
import copy
import os
import sys
import argparse

##### ESTABLISH ARGUMENTS TO FOR PROGRAM USE ON DISTRIBUTED COMPUTING SYSTEM
##### NUISANCE PARAMETERS ARE FOURIER SERIES AMPLITUDE AND PHASE "modes" FOR HISTORICAL REASONS
##### GRADIENTS AQUIRED FROM "SPLITTING" NUISANCE SPACE INTO ARBITRARILY-SELECTED "POSITIVE" AND "NEGATIVE" REGIONS
parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('-i',   '--infiles',    dest='infiles',         nargs = '+',                                        help='LIST OF MULTISIM MONTE CARLO SETS TO PROCESS')
parser.add_argument('-o',   '--outpath',    dest='outpath',                                 default = os.getcwd(),      help='NAME OF THE OUTPUT FILE')
parser.add_argument('-n',   '--max_events', dest='nevents',                                 default=-1,                 help='MAXIMUM NUMBER OF MC EVENTS TO USE')
parser.add_argument('-m',   '--modes',      dest='splitmodes',      nargs = '+',            default=[],                 help='LIST OF MODES TO CALCULATE')
parser.add_argument('-s',   '--splitpoint', dest='splitpoint',                              default=0,                  help='POINT IN NUISANCE SPACE TO SPLIT')
parser.add_argument('-p',   '--phases',     dest='SplitPhases',     action='store_true',                                help='SPLIT ALONG PHASES')
parser.add_argument('-a',   '--amplitudes', dest='SplitAmplitudes', action='store_true',                                help='SPLIT ALONG AMPLITUDES')
args = parser.parse_args()

SplitPhases         = args.SplitPhases
SplitAmplitudes     = args.SplitAmplitudes
splitmodes          = [int(i) for i in args.splitmodes]
nevents             = int(args.nevents)
splitpoint          = float(args.splitpoint)
infiles             = args.infiles
outpath             = args.outpath

##### VERIFY ACCEPTABLE INPUT
if (not (args.infiles)):
    raise RuntimeError("YOU MUST SPECIFY AT LEAST ONE INPUT FILE \n")
if SplitPhases and SplitAmplitudes:
    raise RuntimeError("YOU CANNOT SPLIT PHASES AND AMPLITUDES SIMULTANEOUSLY \n")
if not (SplitPhases or SplitAmplitudes):
    raise RuntimeError("YOU MUST SPECIFY PARAMETER TO SPLIT \n")
if len(splitmodes) == 0: 
    print("NO MODES TO SPLIT SPECIFIED, SPLITTING ALL AVAILABLE MODES \n")

for splitmode in splitmodes:            

    if os.path.isfile(outfile_e):
        print("ERROR: THIS GRADIENT FILE EXISTS, ABORTING.")
        sys.exit()

    ##### READ DATA FROM PRELOADED MONTE CARLO SIMULATION. 
    ##### EXAMPLE LIMITED TO FOURIER PHASES, POSITIVELY-PERTURBED MODELS (pos), AND DISTRIBUTIONS IN ENERGY (e)
    print("SPLITTING MODE " + str(splitmode))
    splitindex    =    splitmode        

    if SplitPhases:                    
        mc_pos_e        =    my_MC.MuExEnergy[ my_MC.SnowStormPhases[:,splitindex] > splitpoint][:nevents] 
        weight_pos_e    =    my_MC.weights[    my_MC.SnowStormPhases[:,splitindex] > splitpoint][:nevents]

    ##### CONFIGURE THE HISTOGRAM
    pos_hist_e, bins_e = np.histogram(mc_pos_e, bins=e_bins, weights=weight_pos_e)  
    bin_centers_e = (e_bins[:-1] + e_bins[1:]) / 2
    bin_widths_e  = abs((e_bins[1:] - e_bins[:-1]) / 2)

    ##### CALCULATE WEIGHTED UNCERTAINTIES
    pos_errors_e = np.zeros(len(bin_centers_e)) 
    for i in range(len(mc_pos_e)):
        for j in range(len(e_bins)-1):
            if mc_pos_e[i] > e_bins[j] and mc_pos_e[i] < e_bins[j+1]:
                pos_errors_e[j] = pos_errors_e[j] + (weight_pos_e[i])**2
    
    pos_errors_e = np.sqrt(pos_errors_e)

    ##### GET MONTE CARLO DISTRIBUTIONS
    mc_pos_dist_e, bins_e = np.histogram( mc_pos_e, bins=e_bins, weights=weight_pos_e, normed=False)   

    ##### COLLECT GRADIENT INFORMATION AND SAVE TO FILE
    fname_e = "SplitCounts_Phs_Energy_" + str(splitmode) + ".csv"
    outfile_e = os.path.join(outpath, fname_e) 
    gradient_e=np.array([bin_centers_e,bin_widths_e, mc_pos_dist_e,pos_errors_e,mc_pos_dist_e,pos_errors_e])
    print("SAVING ENERGY GRADIENT TO " +str(outfile_e) + "\n")
    np.savetxt(outfile_e,gradient_e)

print("################## GRADIENT EXTRACTION COMPLETE ################## ")
