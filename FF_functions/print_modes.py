#!/bin/env python                                                                                                                                                             

import os,sys,argparse,fnmatch,time,math
from numpy import *
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

    parser = argparse.ArgumentParser(description='Reads in a .xyz file and generates orca inputs files to determine the intramolecular modes')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_file', help = 'The input file. (currently must be either (i) an xyz with the atom types in the fourth column or (ii) a folder holding optimized fragments and mode '+\
                                             'information generated by frag_gen.py)')

    # Set random seed
    random.seed(1234)

    # Make relevant inputs lowercase
    args=parser.parse_args(argv)

    # Extract Element list and Coord list from the file
    Elements,Geometry = xyz_parse(args.coord_file)

    # Generate adjacency table
    Adj_mat = Table_generator(Elements,Geometry)

    # Check the number of molecules
    Num_mol = mol_count(Adj_mat)
    if Num_mol > 1:
        print "ERROR: {} molecules were discovered. Check the geometry of the input file. Exiting...".format(Num_mol)
        quit()

    # Find linear, branched, and cyclic segments
    Structure = Structure_finder(Adj_mat)

    # Find backbone
    Backbone = Dijkstra(Adj_mat)

    # Find Hybridizations
    Hybridizations = Hybridization_finder(Elements,Adj_mat)

    # Find atom_types
    #Atom_types = id_types(Elements,Adj_mat,2,Hybridizations,Geometry)
    Atom_types = id_types(Elements,Adj_mat,2) #id_types v.062520

    # Print diagnostic
    print "\n{}".format("*"*144)
    print "* {:^140s} *".format("Summary of Atom Types and Topology")
    print "{}".format("*"*144)

    for count_i,i in enumerate(sorted(set(Atom_types))):
        print "{:80s} : {}".format(i,Atom_types.count(i))    

    # Cyclic Properties
    Num_3 = int(ceil(len([ i for i in Structure if i == 3 ])/3.0))
    Num_4 = int(ceil(len([ i for i in Structure if i == 4 ])/4.0))
    Num_5 = int(ceil(len([ i for i in Structure if i == 5 ])/5.0))
    Num_6 = int(ceil(len([ i for i in Structure if i == 6 ])/6.0))
    Num_7 = int(ceil(len([ i for i in Structure if i == 7 ])/7.0))
    Num_8 = int(ceil(len([ i for i in Structure if i == 8 ])/8.0))

    Number_of_rings = Num_3 + Num_4 + Num_5 + Num_6 + Num_7 + Num_8

    print '\nRings: {}'.format(Number_of_rings)

    if Number_of_rings == 0:
        print " "
    if Num_3 > 0:
        print "3-membered rings: {}".format(Num_3)
    if Num_4 > 0:
        print "4-membered rings: {}".format(Num_4)
    if Num_5 > 0:
        print "5-membered rings: {}".format(Num_5)
    if Num_6 > 0:
        print "6-membered rings: {}".format(Num_6)
    if Num_7 > 0:
        print "7-membered rings: {}".format(Num_7)
    if Num_8 > 0:
        print "8-membered rings: {}".format(Num_8)

    print "\n{}".format("*"*144)
    print "* {:^140s} *".format("Searching for Intramolecular Modes")
    print "{}".format("*"*144)
    
    # Look up the FF parameters: Angles Bonds and Dihedrals are a list of lists, holding 
    #the indices of the atoms participating in each object.
    Angles,Bonds,Dihedrals = Find_parameters(Adj_mat,Atom_types)

    print "\n{}".format("*"*144)
    print "* {:^140s} *".format("Topological Summary")
    print "{}".format("*"*144)

    print "\n\t{:<40s} {}".format("Number of Bonds:",len(Bonds))
    print "\t{:<40s} {}".format("Number of Angles:",len(Angles))
    print "\t{:<40s} {}".format("Number of Dihedrals:",len(Dihedrals))
    
    print "\n\t{:<40s} {}".format("Bonds:",", ".join([ "{:20}".format(i) for i in Bonds[:5] ]))
    if len(Bonds) > 5:
        for i in range(int(ceil(float(len(Bonds)-5)/5.0))):
            print "\t{:<40s} {}".format(" ",", ".join([ "{:20}".format(i) for i in Bonds[5*(i+1):5*(i+2)] ]))

    print "\n\t{:<40s} {}".format("Angles:",", ".join([ "{:20}".format(i) for i in Bonds[:5] ]))
    if len(Angles) > 5:
        for i in range(int(ceil(float(len(Angles)-5)/5.0))):
            print "\t{:<40s} {}".format(" ",", ".join([ "{:20}".format(i) for i in Angles[5*(i+1):5*(i+2)] ]))

    print "\n\t{:<40s} {}".format("Dihedrals:",", ".join([ "{:20}".format(i) for i in Dihedrals[:5] ]))
    if len(Dihedrals) > 5:
        for i in range(int(ceil(float(len(Dihedrals)-5)/5.0))):
            print "\t{:<40s} {}".format(" ",", ".join([ "{:20}".format(i) for i in Dihedrals[5*(i+1):5*(i+2)] ]))

    return

def Structure_finder(Adj_mat):

    # The Structure variable holds the highest structure factor for each atom
    # (atoms might be part of several rings, only the largest is documented)
    # Values correspond to the following structural features:
    # 0: terminal (e.g. hydrogens)
    # 1: chain (e.g. methylene)
    # 2: branch point (e.g. Carbon attached to 3 or more other carbons; -2 indicates a possible chiral center based on coordination)
    # 3: 3-membered ring
    # 4: 4-membered ring
    # 5: 5-membered ring
    # 6: 6-membered ring
    # 7: 7-membered ring
    # 8: 8-membered ring
    Structure = array([-1]*len(Adj_mat))
    # Remove terminal sites (sites with only a single length 2
    # self walk). Continue until all of the terminal structure 
    # has been removed from the topology. NOTE: THIS VERSION DIFFERS FROM POLYGEN BY ITS INDIFFERENCE TO HEAD AND TAIL
    Adj_trimmed = copy(Adj_mat)
    ind_trim = [ count for count,i in enumerate(diag(dot(Adj_trimmed,Adj_trimmed))) if i == 1 ]
    Structure[ind_trim]=0

    # Remove terminal sites
    Adj_trimmed[:,ind_trim] = 0
    Adj_trimmed[ind_trim,:] = 0
    ind_trim = []

    # Find branch points (at this point all hydrogens have been removed, all remaining atoms with
    # over two connected neighbors are at least branches)
    branch_ind = [ count for count,i in enumerate(diag(dot(Adj_trimmed,Adj_trimmed))) if i > 2 ]

    while( len(ind_trim) > 0 ):

        # Remove remaining terminal sites to reveal cyclic structures (This time, the head and tail can be removed
        ind_trim = ind_trim + [ count for count,i in enumerate(diag(dot(Adj_trimmed,Adj_trimmed))) if i == 1 ]
        Structure[ind_trim] = 1

        # Remove terminal sites
        Adj_trimmed[:,ind_trim] = 0
        Adj_trimmed[ind_trim,:] = 0
        ind_trim = []

    # Label branches. This has to be done here otherwise it would get overwritten during the while loop
    Structure[branch_ind] = 2

    # Label possible chiral centers (narrow down to 4-centered branch sites (remove sp2 carbon type branches))
    Chiral_ind = [ i for i in branch_ind if sum(Adj_mat[i]) == 4 ]
    Structure[Chiral_ind] = -2

    # Find non-repeating looping walks of various lengths to identify rings
    # Algorithm: Non-repeating walks are conducted over the trimmed adjacency matrix to 
    tmp=zeros([len(Adj_trimmed),len(Adj_trimmed)])
    for i in range(len(Adj_trimmed)):

        # If the structure of the current atom is already known, then it is skipped
#        if Structure[i] in [1,2]:
#            continue

        # Instantiate generation lists. These hold tupels where the first entry is the connected
        # vertex and the second site is the previous vertex) 
        Gen_1 = []
        Gen_2 = []
        Gen_3 = []
        Gen_4 = []
        Gen_5 = []
        Gen_6 = []
        Gen_7 = []
        Gen_8 = []

        # Find 1st generation connections to current atom. (connection site, previous site)
        Gen_1 = [ (count_z,i) for count_z,z in enumerate(Adj_trimmed[i,:]) if z == 1 ]

        # Loop over the 1st generation connections and find the connected atoms. Avoid back hops using the previous site information 
        for j in Gen_1:
            Gen_2 = Gen_2 +  [ (count_z,j[0],j[1]) for count_z,z in enumerate(Adj_trimmed[j[0],:]) if (z == 1 and count_z not in j[0:-1]) ]
            
        # Loop over the 2nd generation connections and find the connected atoms. Avoid back hops using the previous site information
        for k in Gen_2:
            Gen_3 = Gen_3 + [ (count_z,k[0],k[1],k[2]) for count_z,z in enumerate(Adj_trimmed[k[0],:]) if (z == 1 and count_z not in k[0:-1]) ]
        # Find complete loops, store structure factor, and remove looping sequences from Gen_3 (avoids certain fallacious loops)
        del_ind = [ count_k for count_k,k in enumerate(Gen_3) if k[0] == i ]
        if len(del_ind) > 0:
            if Structure[i] in [1,2]:
                Structure[i] = int(str(Structure[i])+str(3))
            else:
                Structure[i]=3
            Gen_3 = [ z for count_z,z in enumerate(Gen_3) if count_z not in del_ind ]

        # Loop over the 3rd generation connections and find the connected atoms. Avoid back hops using the previous site information        
        for l in Gen_3:
            Gen_4 = Gen_4 + [ (count_z,l[0],l[1],l[2],l[3]) for count_z,z in enumerate(Adj_trimmed[l[0],:]) if (z == 1 and count_z not in l[0:-1]) ]
        # Find complete loops, store structure factor, and remove looping sequences from Gen_3 (avoids certain fallacious loops)
        del_ind = [ count_l for count_l,l in enumerate(Gen_4) if l[0] == i ]
        if len(del_ind) > 0:
            if Structure[i] in [1,2]:
                Structure[i] = int(str(Structure[i])+str(4))
            else:
                Structure[i]=4
            Gen_4 = [ z for count_z,z in enumerate(Gen_4) if count_z not in del_ind ]

        # Loop over the 4th generation connections and find the connected atoms. Avoid back hops using the previous site information
        for m in Gen_4:
            Gen_5 = Gen_5 + [ (count_z,m[0],m[1],m[2],m[3],m[4]) for count_z,z in enumerate(Adj_trimmed[m[0],:]) if (z == 1 and count_z not in m[0:-1]) ]
        # Find complete loops, store structure factor, and remove looping sequences from Gen_3 (avoids certain fallacious loops)
        del_ind = [ count_m for count_m,m in enumerate(Gen_5) if m[0] == i ]
        if len(del_ind) > 0:
            if Structure[i] in [1,2]:
                Structure[i] = int(str(Structure[i])+str(5))
            else:
                Structure[i]=5
            Gen_5 = [ z for count_z,z in enumerate(Gen_5) if count_z not in del_ind ]

        # Loop over the 5th generation connections and find the connected atoms. Avoid back hops using the previous site information
        for n in Gen_5:
            Gen_6 = Gen_6 + [ (count_z,n[0],n[1],n[2],n[3],n[4],n[5]) for count_z,z in enumerate(Adj_trimmed[n[0],:]) if (z == 1 and count_z not in  n[0:-1]) ]
        # Find complete loops, store structure factor, and remove looping sequences from Gen_3 (avoids certain fallacious loops)
        del_ind = [ count_n for count_n,n in enumerate(Gen_6) if n[0] == i ]
        if len(del_ind) > 0:
            if Structure[i] in [1,2]:
                Structure[i] = int(str(Structure[i])+str(6))
            else:
                Structure[i]=6
            Gen_6 = [ z for count_z,z in enumerate(Gen_6) if count_z not in del_ind ]

        # Loop over the 6th generation connections and find the connected atoms. Avoid back hops using the previous site information
        for o in Gen_6:
            Gen_7 = Gen_7 + [ (count_z,o[0],o[1],o[2],o[3],o[4],o[5],o[6]) for count_z,z in enumerate(Adj_trimmed[o[0],:]) if (z == 1 and count_z not in o[0:-1]) ]
        # Find complete loops, store structure factor, and remove looping sequences from Gen_3 (avoids certain fallacious loops)
        del_ind = [ count_o for count_o,o in enumerate(Gen_7) if o[0] == i ]
        if len(del_ind) > 0:
            if Structure[i] in [1,2]:
                Structure[i] = int(str(Structure[i])+str(7))
            else:
                Structure[i]=7
            Gen_7 = [ z for count_z,z in enumerate(Gen_7) if count_z not in del_ind ]

        # Loop over the 7th generation connections and find the connected atoms. Avoid back hops using the previous site information
        for p in Gen_7:
            Gen_8 = Gen_8 + [ (count_z,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]) for count_z,z in enumerate(Adj_trimmed[p[0],:]) if (z == 1 and count_z not in p[0:-1]) ]
        # Find complete loops, store structure factor, and remove looping sequences from Gen_3 (avoids certain fallacious loops)
        del_ind = [ count_p for count_p,p in enumerate(Gen_7) if p[0] == i ]
        if len(del_ind) > 0:
            if Structure[i] in [1,2]:
                Structure[i] = int(str(Structure[i])+str(8))
            else:
                Structure[i]=8
            Gen_8 = [ z for count_z,z in enumerate(Gen_8) if count_z not in del_ind ]

    # Any remaining atoms with unassigned structure must be chain atoms connecting cyclic units. 
    Structure[Structure==-1] = 1 

    return Structure

def Find_parameters(Adj_mat, Atom_types, verbose=1):

    # Initialize lists of each instance and type of FF object.
    # instances are stored as tuples of the atoms involved 
    # (e.g., bonds between atoms 1 and 13 and 17 and 5 would be stored as [(1,13),(17,5)] 
    # Similarly, types are stored as tuples of atom types.
    Bonds = []
    Bond_types = []
    Angles = []
    Angle_types = []
    Dihedrals = []
    Dihedral_types = []
    VDW_types = []

    # Find bonds #
    if verbose == 1:
        print "Parsing bonds..."
    for count_i,i in enumerate(Adj_mat):        
        Tmp_Bonds = [ (count_i,count_j) for count_j,j in enumerate(i) if j == 1 and count_j > count_i ]

        # Store bond tuple so that lowest atom *type* between the first and the second atom is placed first
        # and avoid redundant placements
        for j in Tmp_Bonds:
            if Atom_types[j[1]] < Atom_types[j[0]] and (j[1],j[0]) not in Bonds and (j[0],j[1]) not in Bonds:
                Bonds = Bonds + [ (j[1],j[0]) ]
                Bond_types = Bond_types + [ (Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1]) not in Bonds and (j[1],j[0]) not in Bonds:
                Bonds = Bonds + [ (j[0],j[1]) ]
                Bond_types = Bond_types + [ (Atom_types[j[0]],Atom_types[j[1]]) ]


    # Find angles #
    if verbose == 1:
        print "Parsing angles..."
    for i in Bonds:        

        # Find angles based on connections to first index of Bonds
        Tmp_Angles = [ (count_j,i[0],i[1]) for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j != i[1] ]

        # Store angle tuple so that lowest atom *type* between the first and the third is placed first
        # and avoid redundant placements
        for j in Tmp_Angles:
            if Atom_types[j[2]] < Atom_types[j[0]] and (j[2],j[1],j[0]) not in Angles and (j[0],j[1],j[2]) not in Angles:
                Angles = Angles + [(j[2],j[1],j[0])]
                Angle_types = Angle_types + [ (Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2]) not in Angles and (j[2],j[1],j[0]) not in Angles:
                Angles = Angles + [(j[0],j[1],j[2])]
                Angle_types = Angle_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]]) ]

        # Find angles based on connections to second index of Bonds
        Tmp_Angles = [ (i[0],i[1],count_j) for count_j,j in enumerate(Adj_mat[i[1]]) if j == 1 and count_j != i[0] ]

        # Store angle tuple so that lowest atom *type* between the first and the third is placed first
        # and avoid redundant placements
        for j in Tmp_Angles:
            if Atom_types[j[2]] < Atom_types[j[0]] and (j[2],j[1],j[0]) not in Angles and (j[0],j[1],j[2]) not in Angles:
                Angles = Angles + [(j[2],j[1],j[0])]
                Angle_types = Angle_types + [ (Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2]) not in Angles and (j[2],j[1],j[0]) not in Angles:
                Angles = Angles + [(j[0],j[1],j[2])]
                Angle_types = Angle_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]]) ]

        
    # Find dihedrals #
    if verbose == 1:
        print "Parsing dihedrals..."
    for i in Angles:
        
        # Find atoms attached to first atom of each angle
        Tmp_Dihedrals = [ (count_j,i[0],i[1],i[2]) for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in [i[1],i[2]] ]
        
        # Store dihedral tuple so that the lowest atom *type* between the first and fourth is placed first
        # and avoid redundant placements        
        for j in Tmp_Dihedrals:

            # If the first and fourth atoms are equal, then sorting is based on the second and third
            if Atom_types[j[3]] == Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                if Atom_types[j[2]] < Atom_types[j[1]]:
                    Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
                else:
                    Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

            elif Atom_types[j[3]] < Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2],j[3]) not in Dihedrals and (j[3],j[2],j[1],j[0]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

        # Find atoms attached to the third atom of each angle
        Tmp_Dihedrals = [ (i[0],i[1],i[2],count_j) for count_j,j in enumerate(Adj_mat[i[2]]) if j == 1 and count_j not in [i[0],i[1]] ]
        
        # Store dihedral tuple so that the lowest atom *type* between the first and fourth is placed first
        # and avoid redundant placements        
        for j in Tmp_Dihedrals:

            # If the first and fourth atoms are equal, then sorting is based on the second and third
            if Atom_types[j[3]] == Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                if Atom_types[j[2]] < Atom_types[j[1]]:
                    Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
                else:
                    Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

            elif Atom_types[j[3]] < Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                Dihedral_types = Dihedral_types+ [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2],j[3]) not in Dihedrals and (j[3],j[2],j[1],j[0]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]
            
    # Find non-bonded interactions #
#    print "Parsing VDW interactions..."

    # Note that LAMMPS easily calculates cross VDW parameters on the fly
    # using the pair_modify keyword, so usually cross terms do not need to 
    # be included explicitly.
#     if VDW_ij == 1:
#         for i in set(Atom_types):
#             for j in set(Atom_types):
#                 if i<=j and (i,j) not in VDW_types:
#                     VDW_types = VDW_types + [(i,j)]

#                 elif j<i and (j,i) not in VDW_types:
#                     VDW_types = VDW_types + [(j,i)]
#     else:
#         for i in set(Atom_types):
#             VDW_types = VDW_types + [(i,i)]

    # Reduce the type lists down to unique entries
    Bond_types = [ i for i in set(Bond_types) ]
    Angle_types = [ i for i in set(Angle_types) ]
    Dihedral_types = [ i for i in set(Dihedral_types) ]
        
    # Print System characteristics
    if verbose == 1:
        print "\nSystem characteristics (look over for anything suspicious):"
        print "\nAtom_types ({}):\n".format(len(set(Atom_types)))
        for i in sorted(set(Atom_types)):
            print "\t{}".format(i)
        print "\nBond types ({}):\n".format(len(set(Bond_types)))
        for i in sorted(set(Bond_types)):
            print "\t{}".format(i)
        print "\nAngle types ({}):\n".format(len(set(Angle_types)))
        for i in sorted(set(Angle_types)):
            print "\t{}".format(i)
        print "\nDihedral types ({}):\n".format(len(set(Dihedral_types)))
        for i in sorted(set(Dihedral_types)):
            print "\t{}".format(i)
    ##############################################################

    return Angles,Bonds,Dihedrals

if __name__ == "__main__":
   main(sys.argv[1:])
