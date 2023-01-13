#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *

# NOTE: the order of the imports matters until I resolve the namespace issues (XXX DON'T USE import *)
import random
from numpy import *
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in CHARMM output pdb file and transform it to Gromacs compatible pdb which has time and volume ')

    #required (positional) arguments                                                                                                  

    parser.add_argument('-filename', dest='pdbfile', default='traj.pdb',
                        help = 'input pdb filename, default: traj.pdb')


    parser.add_argument('-print_freq', dest='print_freq', default=1000,
                        help = 'print out frequency of charmm, default: 1000 steps (1ps)')

    parser.add_argument('-o', dest='outfile', default='traj_time.pdb',
                        help = 'output pdb filename, default: traj_time.pdb')

    # Make parse inputs
    args=parser.parse_args(argv)

    if os.path.isfile(args.pdbfile) is False:
         print("ERROR: {} file not found".format(args.pdbfile))
         quit()
    if os.path.isfile('volume.txt') is False:
         print("ERROR: volume.txt not found, this is needed")
         quit()

    print_freq = int(args.print_freq)
    volume = []
    with open('volume.txt','r') as f:
         for lines in f:
            fields = lines.split()
            volume += [float(fields[1])]

    with open(args.pdbfile,'r') as f:
         with open(args.outfile,'w') as g:
            count = 0
            for lines in f:
               fields = lines.split()
               if fields[0] == 'REMARK' and fields[1] == 'SETUP':
                  g.write(lines)
                  g.write("TITLE     Ethanol in water t= {:10.5f} step= {:}\n".format(count*print_freq/1000,count*print_freq)) 
                  g.write("REMARK    THIS IS A SIMULATION BOX\n")
                  g.write("CRYST1   {:<6.3f}   {:<6.3f}   {:<6.3f}  90.00  90.00  90.00 P 1           1 \n".format(volume[count],volume[count],volume[count]))
                  g.write("MODEL        1\n")
                  continue
               elif fields[0] == 'ATOM':
                  g.write(lines)
                  continue
               elif fields[0] == 'TER':
                  g.write('TER \n')
                  g.write('ENDMDL\n')
                  count += 1
    print("TOTAL {} frames converted!".format(len(volume)))

    return


if __name__ == "__main__":
   main(sys.argv[1:])
