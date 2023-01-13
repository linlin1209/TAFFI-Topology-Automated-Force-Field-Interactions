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

    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file (usually outputted from an intramolecular mode parse) '+\
                                                 'and writes inputs for a LAMMPS job to generate configurations for the parsing VDW parameters.')


    #required (positional) arguments                                                                                                  
    parser.add_argument('pdbfile', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')

    parser.add_argument('-o', dest='outfile', default='',
                        help = 'output pdb filename, default: name of xyz file')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)
   

    if args.outfile == '':
      args.outfile = args.pdbfile.split('.')[-2].split('/')[-1] + '_pbc.pdb'

    # start writing
    f = open(args.outfile,'w') 
    with open(args.pdbfile,'r') as g:
      for lc,lines in enumerate(g):
         fields = lines.split()

         if len(fields)>2 and fields[0] == 'CRYST1':
            box = [ float(i) for i in fields[1:4]]
      
         if len(fields)>2 and fields[0] == 'ATOM':
           aindex = int(fields[1])
           molindex = int(fields[4])
           aname = fields[2]
           molname = fields[3]
           Geo = [float(i) for i in fields[5:8]]
           for count_i,i in enumerate(Geo):
               if i < 0: Geo[count_i] += box[count_i]
           f.write("{:<4s}  {:>5d}  {:<4s} {:>3s} {:>4d}  {:>8.3f} {:>8.3f} {:>8.3f} {:>6.2f} {:>6.2f}\n".format('ATOM',aindex,aname,molname,molindex,Geo[0],Geo[1],Geo[2],1,0))
         else:
            f.write(lines)
         
         
         
    print("Successfully convert {} to {}".format(args.pdbfile,args.outfile))



    return


if __name__ == "__main__":
   main(sys.argv[1:])
