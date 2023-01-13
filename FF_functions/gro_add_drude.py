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
    parser.add_argument('grofile', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)
   
    title = args.grofile.split(".")[0]

    new_lines = []
    tmp_store = []
    with open(args.grofile,'r') as f:
      atom_count = 0
      pre_resi = 1
      for lc,lines in enumerate(f):
         fields = lines.split()
         if lc == 0:
            new_lines.append(lines)
            #g.write(lines)
            continue
         if lc == 1:
            new_lines.append(lines)
            natoms = int(fields[0])
            continue
         elif len(fields) < 4:
            box = float(fields[0])
            continue
         else:
            resi = int(fields[0].split("L")[0])
            if pre_resi != resi:
               for line in tmp_store:
                  tmp_fields = line.split() 
                  atom_type = tmp_fields[1]
                  tmp_line = "{:5d}{:5s}{:5s}{:5d} {:< 8.3f}{:< 8.3f}{:< 8.3f}\n".format(pre_resi,'LIN','D'+atom_type,atom_count+1,float(tmp_fields[3]),float(tmp_fields[4]),float(tmp_fields[5]))
                  new_lines.append(tmp_line)
                  atom_count += 1
               tmp_store = []
            pre_resi = resi
            atom_type = fields[1]
            element = split_word(fields[1])[0]
            tmp_line = "{:5d}{:5s}{:5s}{:5d} {:< 8.3f}{:< 8.3f}{:< 8.3f}\n".format(resi,'LIN',atom_type,atom_count+1,float(fields[3]),float(fields[4]),float(fields[5]))
            new_lines.append(tmp_line)
            atom_count += 1
            # we want to add Drude after "all" real atoms in a single residue, so we only write Drude after we reach 
            # the end of that residue, before that, we temporarily store the lines of heavt atoms
            if element != 'H':
               tmp_store.append(lines)
      # add Drude for the last residue 
      for line in tmp_store:
         tmp_fields = line.split() 
         atom_type = tmp_fields[1]
         tmp_line = "{:5d}{:5s}{:5s}{:5d} {:< 8.3f}{:< 8.3f}{:< 8.3f}\n".format(pre_resi,'LIN','D'+atom_type,atom_count+1,float(tmp_fields[3]),float(tmp_fields[4]),float(tmp_fields[5]))
         new_lines.append(tmp_line)
         atom_count += 1
          
               
      tmp_line = "{} {} {}".format(box,box,box)
      new_lines.append(tmp_line)
      new_lines[1] = "{}\n".format(atom_count)

    with open(title+'_drude.gro','w') as g:
      for line in new_lines:
         g.write(line)

    
    quit()
               
            
   
    title = args.xyzfile.split('.')[0] + '.gro'
    with open(title,'w') as f:
         f.write("{:20s}\n".format("Molecule from xyz_to_gro"))
         if args.polar_flag:
            drude_num = len([i for i in Elements if i != 'H'])
            f.write("{}\n".format(len(Elements)+drude_num))
         else:
            f.write("{}\n".format(len(Elements)))
         for count_i,i in enumerate(Atom_types):
            f.write("{:5d}{:5s}{:5s}{:5d} {:< 8.3f}{:< 8.3f}{:< 8.3f}\n".format(1,'LIN',Elements[count_i]+str(count_i+1),count_i+1,Geometry[count_i][0]/10,Geometry[count_i][1]/10,Geometry[count_i][2]/10))
         if args.polar_flag:
            for count_i,i in enumerate(Atom_types):
               if Elements[count_i] != 'H':
                  f.write("{:5d}{:5s}{:5s}{:5d} {:< 8.3f}{:< 8.3f}{:< 8.3f}\n".format(1,'LIN','D'+Elements[count_i]+str(count_i+1),count_i+1,Geometry[count_i][0]/10,Geometry[count_i][1]/10,Geometry[count_i][2]/10))
         
         box = amax(Geometry,axis=0) - amin(Geometry,axis=0)
         f.write("{} {} {}".format(box[0]/10,box[1]/10,box[2]/10))
    print("Successfully convert {} to {}".format(args.xyzfile,title))
    return

def split_word(word): 
    return [char for char in word]  


if __name__ == "__main__":
   main(sys.argv[1:])
