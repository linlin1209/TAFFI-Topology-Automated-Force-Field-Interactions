#!/bin/env python                                                                                                                                                             
# Composed: 7-14-14 
# Description: Reads in a .xyz file and generates a connectivity table

import sys,argparse
import os
from numpy import *
import time
import math
from copy import copy
from math import sqrt,sin,cos,tan,factorial
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile
from pylab import *
import random
from copy import deepcopy

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a .xyz file and prints formatted versions fo the adjacency matrix and atom types to screen.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_file', help = 'The input file. (currently must be an xyz with the atom types in the fourth column.')

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'Controls the bond search depth for identifying unique groups (default: 2)')


    # Make relevant inputs lowercase
    args=parser.parse_args(argv)
    args.gens = int(args.gens)
    # Check that the input is an .xyz file.
    if args.coord_file.split('.')[-1] != 'xyz':
        print "ERROR: Check to ensure that the input file is in .xyz format."
        return

    # Extract Element list and Coord list from the file
    Elements,Geometry,Atom_types = xyz_parse(args.coord_file)

    # Generate adjacency table
    Adj_mat = Table_generator(Elements,Geometry)

    print "\nElements (list-friendly):\n"
    print Elements

    print "\nAdjacency matrix (numpy-friendly):\n"
    label='array(['
    label += ','.join([ '[{}]'.format(','.join([str(int(j)) for j in i])) for i in Adj_mat ])
    label+='])'
    print label

    print "\nAdjacency matrix with elemental basis (formatted):\n"
    print "   "+" ".join(["{:<2s}".format(i) for i in Elements])
    for count_i,i in enumerate(Adj_mat):
        print "{:<2s} ".format(Elements[count_i])+"  ".join([str(int(j)) for j in i ])

    print "\nAdjacency matrix with elemental basis (unformatted):\n"
    print "  "+" ".join([i for i in Elements])
    for count_i,i in enumerate(Adj_mat):
        print Elements[count_i]+" "+" ".join([str(int(j)) for j in i ])

    print "\nAdjacency matrix Latex:"
    print "\\begin{array}{c} "+" \\\ ".join(Elements)+" \end{array}"
    N = len(Elements) + 1 # the syntax is N exclusive, hence + 1
    latex_string =  r"\left[ \begin{array}{@{}*{"+"{}".format(N)+"}{c}@{}} "
    for count_i,i in enumerate(Adj_mat):
        latex_string += " & ".join([str(int(j)) for j in i ]) + " \\\ "
    latex_string = latex_string[:-4]
    latex_string += r"\end{array} \right]"
    print latex_string


    ###### START OF ADJACENCY LIST COMMANDS ########

    # Assemble Adj List
    adj_lst = []
    for count_i,i in enumerate(Adj_mat):
        lst2 = []
        for count_j,j in enumerate(i):
            if j == 1:
                lst2 += [count_j]
        adj_lst += [lst2]

    print "\nAdjacency list (list-friendly):\n"
    label="["
    for count_i,i in enumerate(Adj_mat):
        if count_i == 0:
            label += "["
        else:
            label += ",["
        first = 0
        for count_j,j in enumerate(i):
            if j == 1:
                if first == 0:
                    label += "{:d}".format(count_j)
                    first = 1                    
                else:
                    label += ",{:d}".format(count_j)
        label += "]"
    label += "]"
    print label

    print "\nAdjacency list (formatted):\n"
    for count_i,i in enumerate(adj_lst):
        print "{} {}".format(Elements[count_i],' '.join([ "{}".format(j) for j in i ]))

    print "\nAdjacency list (Latex):"
    print "\\begin{array}{c} "+" \\\ ".join(Elements)+" \end{array}"
    N = len(Elements) + 1 # the syntax is N exclusive, hence + 1
    latex_string =  r"\left[ \begin{array}{@{}*{"+"{}".format(N)+"}{c}@{}} "
    for count_i,i in enumerate(adj_lst):
#        print "\\begin{array}{c} "+" \\\ ".join([ "{}".format(j) for j in i ])+" \end{array}"
        latex_string += " & ".join([str(int(j)) for j in i ]) + " \\\ "
#    latex_string = latex_string[:-4]
    latex_string += r"\end{array} \right]"
    print latex_string


    # Find Hybridizations
    Hybridizations = Hybridization_finder(Elements,Adj_mat)

    # Find atom_types
    Atom_types = id_types(Elements,Adj_mat,args.gens,Hybridizations,Geometry)

    print "\nAtom types (indexed to adjacency matrix):"
    for i in Atom_types:
        print i

    

# # Generates the adjacency matrix 
# # OLD METHOD BASED ON HEURISTIC RADII
# def Table_generator(Elements,Geometry):
    
#     Radii = {'H' :1.4, 'He':1.5, 
#              'Li':2.0, 'Be':1.6,                                                                                                     'B' :1.7, 'C' :1.7, 'N' :1.7, 'O' :1.6, 'F' :1.7, 'Ne':1.6,
#              'Na':2.2, 'Mg':1.9,                                                                                                     'Al':2.0, 'Si':2.1, 'P' :2.1, 'S' :2.0, 'Cl':1.9, 'Ar':1.7,
#              'K' :2.6, 'Ca':2.2, 'Sc':2.4, 'Ti':2.4, 'V' :2.6, 'Cr':2.5, 'Mn':2.6, 'Fe':2.5, 'Co':2.4, 'Ni':2.4, 'Cu':2.5, 'Zn':2.5, 'Ga':2.1, 'Ge':2.1, 'As':2.2, 'Se':2.1, 'Br':2.1, 'Kr':1.8,
#              'Rb':2.9, 'Sr':2.5, 'Y' :2.7, 'Zr':2.6, 'Nb':2.7, 'Mo':2.6, 'Tc':2.5, 'Ru':2.6, 'Rh':2.5, 'Pd':2.5, 'Ag':2.6, 'Cd':2.7, 'In':2.4, 'Sn':2.4, 'Sb':2.4, 'Te':2.4, 'I' :2.3, 'Xe':1.9,
#              'Cs':3.2, 'Ba':2.7, 'La':2.6, 'Hf':2.6, 'Ta':2.7, 'W' :2.6, 'Re':2.6, 'Os':2.6, 'Ir':2.6, 'Pt':2.6, 'Au':2.5, 'Hg':2.6, 'Tl':2.4, 'Pb':2.5, 'Bi':2.5, 'Po':2.5, 'At':2.5, 'Rn':2.0,
#              'default' : 1.8 }

#     # Print warning for uncoded elements.
#     for i in Elements:
#         if i not in Radii.keys():
#             print "ERROR: The geometry contains an element ({}) that the Table_generator function doesn't have bonding information for. This needs to be directly added to the Radii".format(i)+\
#                   " dictionary before proceeding. Exiting..."
#             quit()


#     # Generate distance matrix holding atom-atom separations (only save upper right)
#     Dist_Mat = triu(cdist(Geometry,Geometry))
    
#     # Find plausible connections
#     x_ind,y_ind = where( (Dist_Mat > 0.0) & (Dist_Mat < max([ Radii[i] for i in Radii.keys() ])) )

#     # Initialize Adjacency Matrix
#     Adj_mat = zeros([len(Geometry),len(Geometry)])

#     # Iterate over plausible connections and determine actual connections
#     for count,i in enumerate(x_ind):
        
#         # Determine connection based on the longer of the two radii
#         if Dist_Mat[i,y_ind[count]] < max(Radii[Elements[i]],Radii[Elements[y_ind[count]]]):
#             Adj_mat[i,y_ind[count]]=1

#     # Hermitize Adj_mat
#     Adj_mat=Adj_mat + Adj_mat.transpose()

#     return Adj_mat

# Generates the adjacency matrix 
def Table_generator(Elements,Geometry):

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    Radii = {  'H':2.57, 'He':2.10,\
              'Li':2.18, 'Be':2.45,                                                                                                                'B':3.64,  'C':3.43,  'N':3.26,  'O':3.12,  'F':3.00, 'Ne':2.89,\
              'Na':2.66, 'Mg':2.69,                                                                                                               'Al':4.01, 'Si':3.83,  'P':3.69,  'S':3.59, 'Cl':3.52, 'Ar':3.45,\
               'K':3.40, 'Ca':3.03, 'Sc':2.94, 'Ti':2.83,  'V':2.80, 'Cr':2.69, 'Mn':2.64, 'Fe':2.59, 'Co':2.56, 'Ni':2.52, 'Cu':3.11, 'Zn':2.46, 'Ga':3.90, 'Ge':3.81, 'As':3.77, 'Se':3.75, 'Br':3.73, 'Kr':3.69,\
              'Rb':3.67, 'Sr':3.24,  'Y':2.98, 'Zr':2.78, 'Nb':2.82, 'Mo':2.72, 'Tc':2.67, 'Ru':2.64, 'Rh':2.61, 'Pd':2.58, 'Ag':2.80, 'Cd':2.54, 'In':3.98, 'Sn':3.91, 'Sb':3.94, 'Te':3.98,  'I':4.01, 'Xe':3.92,\
              'Cs':4.02, 'Ba':3.30, 'La':2.80, 'Hf':2.80, 'Ta':2.82,  'W':2.73, 'Re':2.63, 'Os':2.78, 'Ir':2.53, 'Pt':2.45, 'Au':2.93, 'Hg':2.41, 'Tl':3.87, 'Pb':3.83, 'Bi':3.89, 'Po':4.20, 'At':4.23, 'Rn':4.25,\
              'default' : 1.8 }
    
    # Scale factor is used for determining the bonding threshold. The heuristic 0.6 is used
    scale_factor = 0.55

    # Print warning for uncoded elements.
    for i in Elements:
        if i not in Radii.keys():
            print "ERROR: The geometry contains an element ({}) that the Table_generator function doesn't have bonding information for. This needs to be directly added to the Radii".format(i)+\
                  " dictionary before proceeding. Exiting..."
            quit()


    # Generate distance matrix holding atom-atom separations (only save upper right)
    Dist_Mat = triu(cdist(Geometry,Geometry))
    
    # Find plausible connections
    x_ind,y_ind = where( (Dist_Mat > 0.0) & (Dist_Mat < max([ Radii[i]*scale_factor for i in Radii.keys() ])) )

    # Initialize Adjacency Matrix
    Adj_mat = zeros([len(Geometry),len(Geometry)])

    # Iterate over plausible connections and determine actual connections
    for count,i in enumerate(x_ind):
        
        # Assign connection if the ij separation is less than the UFF-sigma value times the scaling factor
        if Dist_Mat[i,y_ind[count]] < (Radii[Elements[i]]+Radii[Elements[y_ind[count]]])/2.0*scale_factor:
            Adj_mat[i,y_ind[count]]=1

    # Hermitize Adj_mat
    Adj_mat=Adj_mat + Adj_mat.transpose()

    return Adj_mat

def xyz_parse(input):

    # Open file and read contents into variable
    with open(input,'r') as f:
            content=f.readlines()

    # Find number of atoms and initialize various matrices
    Atom_Number = int(content[0].split()[0])
    Elements = ['']*Atom_Number
    Geometry = zeros([Atom_Number,3])
    Atom_types = [0]*Atom_Number

    # Iterate over the remainder of contents and read the
    # geometry and elements into variable. Note that the
    # first two lines are considered a header
    count=0
    for lines in content[2:]:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        # Write geometry containing lines to variable
        if len(fields) > 3:
            Elements[count]=fields[0]
            Geometry[count,:]=array([float(fields[1]),float(fields[2]),float(fields[3])])

            # If atom type information is stored in the fifth column it is saved to variable
            if len(fields) >= 5:
                Atom_types[count] = fields[4]

            # If the nuclear charge is stored in the sixth column it is saved to variable
            if len(fields) >= 6:
                Charges[count]=float(fields[5])
            count = count + 1

    return Elements,Geometry,Atom_types

def Hybridization_finder(Elements,Adj_mat):

    Hybridizations = ['X']*len(Elements)
    fail_flag = 0
    
    for i in range(len(Hybridizations)):

        # Find connected atoms and coordination number
        Coord_ind = [ count_j for count_j,j in enumerate(Adj_mat[i]) if j == 1 ]
        Coord_number = len(Coord_ind)

        # Based upon element and coordination number determine hybridization
        if Elements[i] in ['B','Al','C','Si','Ge']:
            if Coord_number == 2:
                Hybridizations[i] = 'sp'
            elif Coord_number == 3:
                Hybridizations[i] = 'sp2'
            elif Coord_number == 4:
                Hybridizations[i] = 'sp3'

        elif Elements[i] in ['N','P','As']:
            if Coord_number == 1:
                Hybridizations[i] = 'sp'
            elif Coord_number == 2:
                Hybridizations[i] = 'sp2'
            elif Coord_number == 3:
                Hybridizations[i] = 'sp3'

        elif Elements[i] in ['O','S','Se']:
            if Coord_number == 1:
                Hybridizations[i] = 'sp2'
            elif Coord_number == 2:
                Hybridizations[i] = 'sp3'

        elif Elements[i] in ['H','Li','Na','K','Rb','Cs','Fr','Be','Mg','Ca','Sr','Ba','Ra']:
            if Coord_number == 1:
                Hybridizations[i] = 's' 

        elif Elements[i] in ['F','Cl','Br','I']:
            if Coord_number == 1:
                Hybridizations[i] = 'sp3'

        else: 
            print "ERROR: hybridization of atom {} could not be determined! Exiting...".format(i)
            fail_flag = 1

    if fail_flag == 1:
        quit()

    return Hybridizations




# Description: This function takes in an element list and adjacency matrix and returns 
#              an index matches list of topological atom types. Here "topological" means
#              the bonding configurations around each atom, going out to "gens" bonds deep.
#              The function uses a canonicalization algorithm to ensure that the same 
#              topological name is returned regardless of the element ordering. 
# 
# Inputs:      Elements:   A list of atomic symbols corresponding to the atom types in the geometry
#                          This list is indexed matched to the adjacency matrix.
#              Adj_mat:    A matrix where each column corresponds to an atom in the geometry, 1's 
#                          correspond to bonds between atoms and 0's to no bonds between atoms.
#
# Returns      Atom_types: A list of topologically named atom types                           
def id_types(Elements,Adj_mat,gens=3,hybridizations=[],geo=[]):

    # Initialize periodic table
    periodic = { "h": 1,  "he": 2,\
                 "li":3,  "be":4,                                                                                                      "b":5,    "c":6,    "n":7,    "o":8,    "f":9,    "ne":10,\
                 "na":11, "mg":12,                                                                                                     "al":13,  "si":14,  "p":15,   "s":16,   "cl":17,  "ar":18,\
                 "k":19,  "ca":20,  "sc":21,  "ti":22,  "v":23,  "cr":24,  "mn":25,  "fe":26,  "co":27,  "ni":28,  "cu":29,  "zn":30,  "ga":31,  "ge":32,  "as":33,  "se":34,  "br":35,  "kr":36,\
                 "rb":37, "sr":38,  "y":39,   "zr":40,  "nb":41, "mo":42,  "tc":43,  "ru":44,  "rh":45,  "pd":46,  "ag":47,  "cd":48,  "in":49,  "sn":50,  "sb":51,  "te":52,  "i":53,   "xe":54,\
                 "cs":55, "ba":56,            "hf":72,  "ta":73, "w":74,   "re":75,  "os":76,  "ir":77,  "pt":78,  "au":79,  "hg":80,  "tl":81,  "pb":82,  "bi":83,  "po":84,  "at":85,  "rn":86}

    # Formatted as a list of lists. Each list corresponds to a generation in the tree
    # and each element in the list is a three element tuple whose first index points to its
    # its parent in the previous generation, whose second index is its index in the
    # Adj_mat/Elements list, and whose third element is its atomic number. The zeroth generation is the atom itself
    trees = [ [ [] ] for i in Elements ]

    for count_i,i in enumerate(Adj_mat):

        # Initialize a list of added atoms so that none get added twice (looping topologies)
        added_atoms = [count_i]

        # Fill in first slot with the parent atom's index
        trees[count_i][0] = [(count_i,count_i,periodic[Elements[count_i].lower()])]

        # Create a list for the first generation connections
        trees[count_i].append([])

        # Iterate over connections and add them to the first generation list
        for count_a,a in enumerate(i):
            
            # If a connection exists add the index to the list
            if a == 1 and count_a not in added_atoms:
                trees[count_i][1] += [(count_i,count_a,periodic[Elements[count_a].lower()])]
                added_atoms += [count_a]

        # Recursively add branch information based on the number of generations requested by user
        for g in range(1,gens):
        
            # Create a list for the second generation connections
            trees[count_i].append([])

            # Iterate over the adj_mat rows for the first generation connections
            for count_j,j in enumerate(trees[count_i][g]):
                for count_a,a in enumerate(Adj_mat[j[1]]):

                    # If a connection exists that is not the connection from the parent generation add the index to the list
                    if a == 1 and count_a != j[0] and count_a not in added_atoms:
                        trees[count_i][g+1] += [(j[1],count_a,periodic[Elements[count_a].lower()])]
                        added_atoms += [count_a]

    # Sort trees: each generation is sorted based on atomic number,
    # where equal, sort is based on the number of attached atoms, 
    # where equal, sort is based on the atomic number of attached atoms
    # where equal, sort continues based on these criteria through all generations
    # until a point of difference arises. If none is found, ordering is left as is.
    for count_i,i in enumerate(trees):
        trees[count_i] = sort_tree(i)

    # Find E/Z stereocenters (if they exist and the necessary 
    # hybridizations and geometry have been supplied) 
    if hybridizations != [] and geo != []:
        EZ_labels = find_EZ(trees,hybridizations,geo)

    # Nest the trees according to topology 
    nested_trees = gen_nested_trees(trees)

    # Created nested atom lists according to topology
    # (same structure as nested tree, except tuples are 
    #  replaced with the atomic number only)
    nested_atoms = list(nested_trees)
    for count_i,i in enumerate(nested_atoms):
        nested_atoms[count_i] = strip_trees(i)

    # Create atom_type labels from the nested atom lists
    atom_types = ['X']*len(nested_atoms)        
    for count_i,i in enumerate(nested_atoms):
        atom_types[count_i] = '['.join(str(i).split(', ['))
        atom_types[count_i] = EZ_labels[count_i] + atom_types[count_i]

    return atom_types

# Sort trees: canonicalization algorithm follows cahn-ingold-prelog rules
# each generation is sorted based on atomic number,
# where equal, each atom is replaced by a list of attached atoms and comparisons are made between each sublist
# where equal, each sublist is replaced by a list of attached atoms in the next generation
# this process continues until a difference is found. 
# Note 1: the common multiple bond rule is omitted as it is unneccesary for canonicalization)
# Note 2: there is subtlety when applying the CIP rules to highly branched structures that is associated 
#         with the order in which subshells are compared. Neither Mcarrey, Wade, or Sorrell's books
#         are explicit on the matter. Here I always compare lists based upon highest MW preceeding
#         atom.
def sort_tree(tree):

    # Sort all generations by atomic weight
    for count_i,i in enumerate(tree):
        tree[count_i].sort(key=lambda x:x[2],reverse=True)

    # Iterate over each generation ranking each node's priority
    for count_i,i in enumerate(tree):        

        # If there are empty generations, avoid sorting
        if len(i) == 0: continue

        # Find nodes with identical atomic weights and the number of connections of each node in the next generation
        # Initialize rank list
        rank = [-1]*len(i)
        rank_counter = 0
        
        for j in range(len(i)):

            # If unassigned, assign a rank
            if rank[j] == -1:
                rank[j] = rank_counter
                
            # Check for matches based on atomic number and assign equal ranks to those
            for k in range(len(i)):
                if i[j][2] == i[k][2] and j != k and rank[k] == -1:
                    rank[k]=rank_counter

            # Increment rank_counter based on the number of unresolved pairs that had the same priority
            rank_counter += len([ k for k in rank if k == rank_counter ])

        # Recursively search through the depth of the tree to resolve priority differences
        depth = count_i
        for j in set(rank):

            # Find connections that need a deeper comparison to resolve priority
            idx = [ count_k for count_k,k in enumerate(rank) if k == j ]

            # Skip already resolved ranks
            if len(idx) == 1:
                continue

            # Initialize tmp_ranks (holds the sublist of ranks for the atoms being resolved)
            tmp_ranks = [0]*len(idx)
            
            # Initialize the sorting list (a list of lists holding sublists of connections for each connection)
            # i.e., the elements in the first list correspond to lists holding the attached atoms
            #       the elements in each sublist hold sublists to compare branches between generations. 
            sorting_list = [ [find_connections(tree,[depth,k])] for k in idx ]

            # Recursively search the tree for differences in priority. 
            # The loop is exited once all of the unresolved branches have been assigned a rank
            while len(idx) > 0:
                                
                # rank the sorting_list (list of lists of lists; each element in the first list
                disp = 0
                for k in sorted(set(tmp_ranks)):

                    # This set of loops guarantees that comparisons are only performed between equally ranked subsets of 
                    # the unresolved set (e.g,, there might be 4 branches that are still unresolved, but 2 are ranked 0,0
                    # and 2 are ranked 2,2. It is important in this case that 0,0 and 2,2 get compared only amongst themselves.
                    loc_list = []
                    sub_sorting_list = []
                    for count_m,m in enumerate(tmp_ranks):
                        if m == k:
                            sub_sorting_list += [sorting_list[count_m]]
                            loc_list += [count_m]

                    # rank_nested_lists: returns the priorities/ranks of the sub_sorting_list
                    sub_ranks = rank_nested_lists(sub_sorting_list)
                    for count_m,m in enumerate(loc_list):
                        tmp_ranks[m] = sub_ranks[count_m]+disp 
                        
                    # the minimum ranking of the next atoms being resolved starts at disp + #of_just_resolved
                    disp += len(sub_ranks)
                
                # Update rank for resolved atoms
                for count_k,k in enumerate(tmp_ranks):
                    if tmp_ranks.count(k) == 1:
                        rank[idx[count_k]] += tmp_ranks[count_k]

                # Update the rank of unresolved atoms:
                # the overall rank of each unresolved atom is incremented for each just
                # resolved atom that is ranked in front of it.
                for count_k,k in enumerate(tmp_ranks):
                    if tmp_ranks.count(k) == 1:
                        for count_m,m in enumerate(tmp_ranks):
                            if tmp_ranks.count(m) != 1 and m > k:
                                rank[idx[count_m]] += 1

                # ID the remaining unresolved branches and remove the
                # rest from sorting_list and idx
                tmp_sorting_list = []
                tmp_idx = []                
                tmp_tmp_ranks = []
                for count_k,k in enumerate(tmp_ranks):
                    if tmp_ranks.count(k) != 1:
                        tmp_sorting_list += [sorting_list[count_k]]
                        tmp_idx += [idx[count_k]]
                        tmp_tmp_ranks += [k]
                sorting_list = tmp_sorting_list
                idx = tmp_idx
                tmp_ranks = tmp_tmp_ranks

                # Repopulate sorting_list from the next generation
                tmp_sorting_list = [ [] for k in range(len(idx)) ]
                for count_k,k in enumerate(sorting_list):
                    for count_m,m in enumerate(k):
                        for count_n,n in enumerate(m):

                            # if n is an empty node then each subranch mutually deadended here and it can be skipped
                            if n == []:
                                continue
                            tmp_sorting_list[count_k] += [ find_connections(tree,n[1]) ]
                sorting_list = tmp_sorting_list

        # Sort the tree based on the final rankings
        tmp = [[]]*len(tree[count_i])
        for count_j,j in enumerate(rank):
            tmp[j] = tree[count_i][count_j]
        tree[count_i] = tmp

    return tree

def rank_nested_lists(lol):    

    # matrix based sort algorithm. This is knowingly suboptimal, but it 
    # gives much more flexibility during the sort than standard algorithms. 
    order=zeros([len(lol),len(lol)])
    rank = [0]*len(lol)
    terminal_list = []
    for count_i,i in enumerate(lol):
        
        # If there were no connections in the current generation then comparisons have 
        # been exhausted and ranking should be set in this round.
        if len(flatten(i)) == 0:
            terminal_list += [count_i]
            order[count_i,:] = 0
            continue

        for count_j,j in enumerate(lol):

            if count_j > count_i:

                flag=0

                # iterate over sublists in i while performing
                # comparisons with corresponding sublists in j
                for count_k,k in enumerate(i):

                    if flag == 1:
                        break

                    # If there are no sublists left in j and difference has been found, i has priority
                    elif count_k > len(j)-1:
                        order[count_i,count_j]=1
                        order[count_j,count_i]=0
                        break

                    for count_m,m in enumerate(k):

                        # If the sublist in j has run out of values, i has priority
                        if count_m > len(j[count_k])-1:
                            order[count_i,count_j]=1
                            order[count_j,count_i]=0
                            flag = 1
                            break                       

                        # If the sublist element is an empty node (with branching empty nodes can be mixed with full nodes)
                        if i[count_k][count_m] == [] and j[count_k][count_m] != []:
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break

                        # If the sublist element is an empty node (with branching empty nodes can be mixed with full nodes) (inverse)
                        if i[count_k][count_m] != [] and j[count_k][count_m] == []:
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break

                        # If both sublist elements are empty nodes (with branching empty nodes can be mixed with full nodes)
                        if i[count_k][count_m] == [] and j[count_k][count_m] == []:
                            continue

                        # If there is a difference in the sublist elements assign priority
                        elif i[count_k][count_m][0] > j[count_k][count_m][0]:
                            order[count_i,count_j]=1
                            order[count_j,count_i]=0
                            flag = 1
                            break

                        # If there is a difference in the sublist elements assign priority (inverse)
                        elif i[count_k][count_m][0] < j[count_k][count_m][0]:
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break

                        # If the sublist in i has run out of values and j still has attachments, j has priority
                        elif count_m == len(i[count_k])-1 and len(j[count_k]) > len(i[count_k]):
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break
                            
                    # i has run out of sublists without finding a difference and j still has lists, j has priority
                    if count_k == len(i)-1 and len(j) > len(i) and flag == 0:
                        order[count_i,count_j]=0
                        order[count_j,count_i]=1                        
                        break

                    # if both i and j have run out of sublists without finding a difference then they are identical up to this point
                    elif count_k == len(i)-1 and len(i) == len(j) and flag == 0:
                        order[count_i,count_j]=1
                        order[count_j,count_i]=1
                        break
                    
    # Assign ranks based on entries in the order matrix
    for count_i,i in enumerate(order):
        rank[count_i]=int(len(order)-sum(i)-1)

    # Assign unique ranks for branches that are at their terminus
    # Reverse order shouldn't matter, but it guarrantees that the
    # first discovered gets higher priority
    for count_i,i in enumerate(terminal_list[::-1]):
        rank[i]=rank[i]-count_i

    return rank

def find_connections(tree,loc):
    if len(tree) > loc[0]+1:
        
        # Find connections
        connections = [ [i[2],[loc[0]+1,count_i]] for count_i,i in enumerate(tree[loc[0]+1]) if i[0] == tree[loc[0]][loc[1]][1] ]

        # If no connections are found then a 0-entry is returned
        if connections == []:
            connections = [[]]
        # If connections are found then proceed with sort
        else:
            connections.sort(key=lambda x:x[0],reverse=True)

    # If a next generation doesn't exist then a 0-entry is returned
    else:
        connections = [[]]

    return connections

def flatten(lol):
    flat = []
    for i in lol:
        if type(i) == type([]):
            flat += flatten(i)
        else:
            flat += [i]
    return flat
    
# lol: nested list of lists, where each element is a 3-element tuple
# This function returns a list of lists with the same structure as lol, except 
# that each tuple is replaced by its third element (atomic number in this case)
def strip_trees(lol):
    tmp = []
    for count_i,i in enumerate(lol):
        if type(i) == type([]):
            tmp += [strip_trees(i)]
        else:
            tmp += [i[2]]
    return tmp

# Description: iterate over the nodes in each generation and search the previous generation for matches (using find_loc function)
# insert all matches into the label as a list in the index following the previous generation node.
# Operates on the total tree structure (so, trees is a list of trees)
def gen_nested_trees(trees):
    nested_trees = [ [] for i in range(len(trees)) ]
    for count_i,i in enumerate(trees):

        # Count the number of occupied generations
        num_gens = len([ count_j for count_j,j in enumerate(i) if len(j)>0])

        # Intiialize the nested list with the first generation atom
        nest = [i[0][0]]

        # Loop over generations and unwind the topology
        for j in range(1,num_gens):

            # Reversed list keeps the ordering so that highest priority ends up at the front
            for k in i[j][::-1]:         
                loc = find_loc(k[0],nest)
                if len(loc) == 1:
                    nest.insert(loc[0]+1,[k])
                elif len(loc) == 2:
                    nest[loc[0]].insert(loc[1]+1,[k])
                elif len(loc) == 3:
                    nest[loc[0]][loc[1]].insert(loc[2]+1,[k])
                elif len(loc) == 4:
                    nest[loc[0]][loc[1]][loc[2]].insert(loc[3]+1,[k])
                elif len(loc) == 5:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]].insert(loc[4]+1,[k])
                elif len(loc) == 6:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]].insert(loc[5]+1,[k])
                elif len(loc) == 7:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]].insert(loc[6]+1,[k])
                elif len(loc) == 8:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]][loc[6]].insert(loc[7]+1,[k])
                elif len(loc) == 9:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]][loc[6]][loc[7]].insert(loc[8]+1,[k])
                elif len(loc) == 10:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]][loc[6]][loc[7]][loc[8]].insert(loc[9]+1,[k])
                else:
                    print "Nesting is too deep, ensure that the requested generations do not exceed 10"
                    quit()
        nested_trees[count_i] = nest

    return nested_trees

# Description: Recursive search to id location of the matching node in the list of lists.
#              Returns a list of indices corresponding to the position within the list of
#              lists of the matching element. (e.g., loc=[2,4,3] means lol[2][4][3] is the match)
#
# Inputs:      idx: scalar, sought match for the second element in the tuples held within the list of lists.
#              lol: list of lists (arbitrary nesting) whose lowest level elements are 3-tuples. 
#
# Returns:     loc: a list of indices specifying the location in lol of the matching tuple
def find_loc(idx,lol):
    loc = []                                   # Initialize list to hold the location indices
    for count_i,i in enumerate(lol):           # Iterate over the elements in the list
        if type(i) == type([]):                # If a list is encountered a recursion is initiated
            tmp = find_loc(idx,i)              # recursive search for a match (yields [] is no match is found)
            if tmp != []:                      # If the recursion yields a match then we have a winner
                loc = loc + [count_i] + tmp    # append the list index (count_i) and the loc list returned by the recursion
        elif idx == i[1]:                      # If a match is found in the element-wise search append its location
            loc += [count_i]      
    return loc

def find_EZ(trees,hybridizations,geo):
    
    EZ_labels = ['']*len(trees)
    for count_i,i in enumerate(trees):

        # Check if the anchor atom is a potential stereocenter
        # NOTE: in the current implementation only C,N,Si,and P are considered potential EZ centers)
        if hybridizations[count_i] != 'sp2' or i[0][0][2] not in [6,7,14,15]:
            continue

        # Find indices of potential EZ bonds
        idx = [ j[1] for j in i[1] if hybridizations[j[1]] == 'sp2' and j[2] in [6,7,14,15] ]

        # Iterate over all of the EZ bonds and determine stereochemistry
        for j in idx:

            # Find anchor atom's highest priority connection (excluding the sp2 bonded group)      
            ind_1 = [ k[1] for k in i[1] if k[1] != j ][0]
            
            # Find bonded atom's highest priority connection (excluding the sp2 bonded group)
            ind_4 = [ k[1] for k in trees[j][1] if k[1] != count_i ][0] 

            # Intialize the coordinate array of the relevant dihedral
            xyz = zeros([4,3])
            xyz[0] = geo[ind_1,:]
            xyz[1] = geo[count_i,:]
            xyz[2] = geo[j,:]
            xyz[3] = geo[ind_4,:]

            # Calculate dihedral using call to dihedral_calc
            angle = 180.0*dihedral_calc(xyz)/pi         # Calculate dihedral
            angle = int(180*round(abs(angle)/180.0))    # Round to 0 or 180

            # ID E/Z based on angle
            if angle == 0:
                EZ_labels[count_i] += "Z"
            elif angle == 180:
                EZ_labels[count_i] += "E"

    return EZ_labels

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

# Description: finds the number of disconnected subnetworks in the 
#              adjacency matrix, which corresponds to the number of 
#              separate molecules.
#
# Inputs:      adj_mat: numpy array holding a 1 in the indices of bonded
#                        atom types. 
#
# Returns:     mol_count: scalar, the number of molecules in the adj_mat
def mol_count(adj_mat):
    
    # Initialize list of atoms assigned to molecules and a counter for molecules
    placed_idx = []    
    mol_count = 0

    # Continue the search until all the atoms have been assigned to molecules
    while len(placed_idx)<len(adj_mat):

        # Use sequential elements of the adj_mat as seeds for the spanning network search
        for count_i,i in enumerate(adj_mat):

            # Only proceed with search if the current atom hasn't been placed in a molecule
            if count_i not in placed_idx:

                # Increment mol_count for every new seed and add the seed to the list of placed atoms
                mol_count += 1               
                placed_idx += [count_i]
                
                # Find connections
                idx = [ count_j for count_j,j in enumerate(i) if j==1 and count_j not in placed_idx ]
                
                # Continue until no new atoms are found
                while len(idx) > 0:
                    current = idx.pop(0)
                    if current not in placed_idx:
                        placed_idx += [current]
                        idx += [ count_k for count_k,k in enumerate(adj_mat[current]) if k == 1 and count_k not in placed_idx ]
    return mol_count



if __name__ == "__main__":
   main(sys.argv[1:])
