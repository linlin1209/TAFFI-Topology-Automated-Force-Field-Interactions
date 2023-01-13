#!/bin/env python                                                                                                                                                             

import os,sys,argparse,fnmatch,time,math
import math
from numpy import *
import numpy as np
from copy import copy
from math import sqrt,sin,cos,tan,factorial,acos
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile
from pylab import *
import random
from copy import deepcopy

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *

def main(argv):

    parser = argparse.ArgumentParser(description='take in xyz file and index(4) will return CHARMM format lone pair positions')


    #required (positional) arguments                                                                                                  
    parser.add_argument('xyzfile', help = 'xyz file name')
    parser.add_argument('-index', dest='index', default='0,1,2,3',
                        help = 'index for atom selection')

    # Make parse inputs
    args=parser.parse_args(argv)
    
    args.index = [ int(i) for i in args.index.split(',')]

    Elements,Geometry = xyz_parse(args.xyzfile)
    #geos[2] = (geos[2]+geos[3])/2
    #distance,angle,dihedral = lonepair(Geometry,[6,1,0,2])
    distance,angle,dihedral = lonepair(Geometry,args.index)
    print("distance: {:>4.2f} angle: {:>4.2f} dihedral: {:>4.2f}".format(distance,angle,dihedral))

        
    return

def lonepair(Geometry,index):
    geos = np.zeros([4,3])
    for count_i,i in enumerate(index):
      geos[count_i] = Geometry[i]
    #geos[2] = (geos[2]+geos[3])/2
    distance = norm(geos[1]-geos[0])
    angle = angle_between_normal(geos[0]-geos[1],geos[2]-geos[1]) 
    dihedral = dihedral_calc(geos)/math.pi*180

    return distance,angle,dihedral
   
    

# Description: Calculates the dihedral angle (in radians) for
#              a quadruplet of atoms
#
# Inputs:      xyz: a 4x3 numpy array, where each row holds the 
#                   cartesian position of each atom and row indices
#                   determine placement in the dihedral ( i.e, which
#                   atoms correspond to 1-2-3-4 in the dihedral )
# 
# Returns:    angle: dihedral angle in radians
#
def dihedral_calc(xyz):
    
    # Calculate the 2-1 vector           
    v1 = (xyz[1]-xyz[0]) 
                                                             
    # Calculate the 3-2 vector           
    v2 = (xyz[2]-xyz[1]) 
                                                             
    # Calculate the 4-3 bond vector      
    v3 = (xyz[3]-xyz[2]) 

    # Calculate dihedral (funny arctan2 formula avoids the use of arccos and retains the sign of the angle)
    angle = arctan2( dot(v1,cross(v2,v3))*(dot(v2,v2))**(0.5) , dot(cross(v1,v2),cross(v2,v3)) )
    
    return angle

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between_normal(v1,v2):
# Return acute and obstuse angle
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))/math.pi*180.0
    #angle = np.arccos(np.dot(v1_u, v2_u))/math.pi*180.0
    
    return angle

def bisector(v1,v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    bi = v1_u + v2_u

    return bi

# Create logger to save stdout to logfile
class Logger(object):

    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/paramgen.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    
    def flush(self):
        pass

if __name__ == "__main__":
   main(sys.argv[1:])
