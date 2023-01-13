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

    parser = argparse.ArgumentParser(description='Reads in a xyz file then transform to gro file (fixed format)') 


    #required (positional) arguments                                                                                                  
    parser.add_argument('xyzfile', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')

    parser.add_argument('-o', dest='outfile', default='',
                        help = 'output pdb filename, default: name of xyz file')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)
   
    Elements,Geometry = xyz_parse(args.xyzfile)
    Adj_mat = Table_generator(Elements,Geometry)
    Hybridizations = Hybridization_finder(Elements,Adj_mat)
    gens =2
    Atom_types = id_types(Elements,Adj_mat,gens)
    if args.outfile == '':
      args.outfile = args.xyzfile.split('.')[-2].split('/')[-1] + '.gro'
    with open(args.outfile,'w') as f:
         f.write("{:20s}\n".format("Molecule from xyz_to_gro"))
         if args.polar_flag:
            drude_num = len([i for i in Elements if i != 'H'])
            f.write("{}\n".format(len(Elements)+drude_num))
         else:
            f.write("{}\n".format(len(Elements)))

         for count_i,i in enumerate(Atom_types):
            f.write("{:5d}{:<5s}{:>5s}{:>5d}{:> 8.3f}{:> 8.3f}{:> 8.3f}\n".format(1,'LIN',Elements[count_i]+str(count_i+1),count_i+1,Geometry[count_i][0]/10,Geometry[count_i][1]/10,Geometry[count_i][2]/10))
         if args.polar_flag:
            for count_i,i in enumerate(Atom_types):
               if Elements[count_i] != 'H':
                  f.write("{:5d}{:<5s}{:>5s}{:>5d}{:> 8.3f}{:> 8.3f}{:> 8.3f}\n".format(1,'LIN','D'+Elements[count_i]+str(count_i+1),count_i+1,Geometry[count_i][0]/10,Geometry[count_i][1]/10,Geometry[count_i][2]/10))
         
         box = amax(Geometry,axis=0) - amin(Geometry,axis=0)
         f.write("{} {} {}".format(box[0]/10,box[1]/10,box[2]/10))
    subprocess.call("module load gromacs",shell=True)
    command = 'gmx_mpi editconf -f {}.gro -o {}.gro'.format(args.outfile,args.outfile)
    output = subprocess.Popen(command.split(),stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()[-1]
    print("Successfully convert {} to {}".format(args.xyzfile,args.outfile))
    return


    

if __name__ == "__main__":
   main(sys.argv[1:])
