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
import numpy as np

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file (usually outputted from an intramolecular mode parse) '+\
                                                 'and writes inputs for a LAMMPS job to generate configurations for the parsing VDW parameters.')


    #required (positional) arguments                                                                                                  
    parser.add_argument('xyzfile', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')

    parser.add_argument('-o', dest='outfile', default='',
                        help = 'output pdb filename, default: name of xyz file')

    parser.add_argument('-resname', dest='resname', default='LIN',
                        help = 'resname of the molecule,(no more than 3 characters) default: LIN')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)
   
    Elements,Geometry = xyz_parse(args.xyzfile)
    Adj_mat = Table_generator(Elements,Geometry)
    Hybridizations = Hybridization_finder(Elements,Adj_mat)
    gens =2
    Atom_types = id_types(Elements,Adj_mat,gens)

    resname = ''.join(split_string(args.resname)[:3])

    if args.outfile == '':
      args.outfile = args.xyzfile.split('.')[-2].split('/')[-1] + '.pdb'
    with open(args.outfile,'w') as f:
         f.write("TITLE     Molecule from xyz_to_gro\n")
         f.write("REMARK    THIS IS A SIMULATION BOX\n")
         box = np.amax(Geometry,axis=0) - np.amin(Geometry,axis=0)
         f.write("CRYST1   {:<6.3f}   {:<6.3f}   {:<6.3f}  90.00  90.00  90.00 P 1           1 \n".format(box[0],box[1],box[2]))
         f.write("MODEL        1\n")
         for count_i,i in enumerate(Atom_types):
            f.write("{:<4s}  {:>5d}  {:<4s} {:>3s} {:>4d}  {:>8.3f} {:>8.3f} {:>8.3f} {:>6.2f} {:>6.2f}\n".format('ATOM',count_i+1,Elements[count_i]+str(count_i+1),resname,1,Geometry[count_i][0],Geometry[count_i][1],Geometry[count_i][2],1,0))
         f.write("TER \n")
         f.write("ENDMDL\n")
         
    print("Successfully convert {} to {}".format(args.xyzfile,args.outfile))
    return

def split_string(word): 
    return [char for char in word]  

if __name__ == "__main__":
   main(sys.argv[1:])
