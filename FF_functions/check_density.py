#!/bin/env python                                                                                                                                                             

import sys,argparse,os,time,ast,random,re,fnmatch,matplotlib
matplotlib.use('Agg') # Needed for cluster image generation 
from pylab import *
from numpy import *
from numpy.linalg import norm
from scipy.spatial.distance import cdist
from shutil import move,copyfile
from copy import deepcopy

def main(argv):
    
    parser = argparse.ArgumentParser(description='Extracts pair-wise configurations from a LAMMPS trajectory and generates an input for an orca binding energy calculation.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('datafile', help = 'A lammps datafile whose density is to be calculated.')


    #optional arguments
    parser.add_argument('-d', dest="d_mode", default="number", help = "The density to be calculated. 'number' and 'mass' are valid options. (default: number)")

    # Assign stdin arguments to variables
    args=parser.parse_args(argv)
    args.d_mode = args.d_mode.lower()
    if args.d_mode not in ['number','mass']: 
        print("ERROR: {} is not a valid option for -d. Valid options are 'mass' and 'density'. Exiting...".format(args.d_mode))
        quit()
    if args.d_mode == 'mass':
        print("OOPS, the mass -d option isn't implemented yet. Please implement it!")
        quit()
    if os.path.isfile(args.datafile) is False:
        print("ERROR: could not find datafile {}. Exiting...".format(args.datafile))
        quit()

    box = zeros([3,2])
    with open(args.datafile,'r') as f:
        for lines in f:
            fields = lines.split()
            if len(fields) == 2 and fields[1] == "atoms":
                N_atoms = float(fields[0])
                
            if len(fields) == 4:
                if fields[2] == 'xlo': box[0,:] = array([float(fields[0]),float(fields[1])])
                if fields[2] == 'ylo': box[1,:] = array([float(fields[0]),float(fields[1])])
                if fields[2] == 'zlo': box[2,:] = array([float(fields[0]),float(fields[1])]); break
                    
    volume = (box[0,1]-box[0,0]) * (box[1,1]-box[1,0]) * (box[2,1]-box[2,0])
    print(N_atoms/volume)
    return (N_atoms/volume)


if __name__ == "__main__":
   main(sys.argv[1:])
