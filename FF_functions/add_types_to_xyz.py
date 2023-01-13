#####adddddd
####two
#!/bin/env python                                                                                                                                                             

import sys,argparse,os,time,math
from numpy import *
from math import sqrt,sin,cos,tan,factorial
from scipy import *
from numpy.linalg import *
from shutil import move,copyfile
from fnmatch import fnmatch
from copy import deepcopy
from itertools import permutations

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from kekule import *
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from __builtin__ import any  # any has a namespace conflict with numpy. XXX WE REALLY NEED TO CLEAN UP THESE NAMESPACE ISSUES

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a .xyz file and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    parser.add_argument('-o', dest='outputname', default='frags',
                        help = 'Controls the output folder name for the fragments. (default: frags)')

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'Controls the bond search depth for identifying unique atom types (default: 2)')

    parser.add_argument('-R', dest='recursive_opt', default=False, action='store_const',const=True,
                        help = 'When this flag is present, filenames and wildcards will be searched for recursively.')

    # Make relevant inputs lowercase
    args=parser.parse_args(argv)    

    # Parse additional inputs 
    args.gens = int(args.gens)

    # Find files
    args.coord_files = args.coord_files.split()
    wild_card_files  = [ i for i in args.coord_files if "*" in i ]
    if args.recursive_opt == True:
        args.coord_files = [ i for i in args.coord_files if "*" not in i ]
        args.coord_files += [ dp+"/"+files for i in wild_card_files for dp,dn,fn in os.walk('.') for files in fn if fnmatch(files,i) ] 
    else:
        args.coord_files = [ i for i in args.coord_files if "*" not in i ]
        for i in wild_card_files:
            path = '/'.join(i.split('/')[:-1])
            if len(path) == 0:
                path = "."
            args.coord_files += [ path+"/"+files for files in os.listdir(path) if fnmatch(files,i) ] 

    # Handle "./" condition
    for count_i,i in enumerate(args.coord_files):
        if i[:2] == "./": args.coord_files[count_i] = args.coord_files[count_i][2:]    
    args.coord_files = list(set(args.coord_files))

    # Perform some consistency checks
    if False in [ i.split('.')[-1] == 'xyz' for i in args.coord_files ]:
        print "ERROR: Check to ensure that the input file(s) are in .xyz format."
        return
    if False in [ os.path.isfile(i) for i in args.coord_files ]:
        print "ERROR: Could not find file {}. Exiting...".format(next([ i for i in args.coord_files if os.path.isfile(i) == False ]))

    # Generate base filename for naming output
    Filename = args.outputname

    # Initialize dictionary and lists for holding modes and geometry information
    geo_dict    = {}
    bonds       = []
    bond_types  = []
    bond_files  = []
    angles      = []
    angle_types = []
    angle_files = []
    dihedrals   = []
    dihedral_types = []
    dihedral_files = []

    # Loop over the input structures, generate the all trans-minimized conformer,  and parse the modes
    print "\n"+"*"*120
    print "* {:^116s} *".format("Adding atomtypes to the following xyz files")
    print "*"*120 + "\n"    
    for file in args.coord_files:        

        # Print the name
        print "  {}".format(file)

        # Initialize subdictionary
        geo_dict[file] = {}

        # Extract Element list and Coord list from the file
        geo_dict[file]["elements"],geo_dict[file]["geo"] = xyz_parse(file)

        # Generate adjacency table
        geo_dict[file]["adj_mat"] = Table_generator(geo_dict[file]["elements"],geo_dict[file]["geo"])

        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        hybridizations = Hybridization_finder(geo_dict[file]["elements"],geo_dict[file]["adj_mat"])

        # Find atom types
        #geo_dict[file]["atom_types"] = id_types(geo_dict[file]["elements"],geo_dict[file]["adj_mat"],args.gens,hybridizations,geo_dict[file]["geo"])
        geo_dict[file]["atom_types"] = id_types(geo_dict[file]["elements"],geo_dict[file]["adj_mat"],args.gens) # id_types v.062520
 
        # Write xyz with types
        with open(Filename+'_wtypes.xyz','w') as f:
            f.write("{}\n\n".format(len(geo_dict[file]["elements"])))
            for count_i,i in enumerate(geo_dict[file]["elements"]):
                f.write("{:<20s} {:< 20.8f} {:< 20.8f} {:< 20.8f} {}\n".format(i,geo_dict[file]["geo"][count_i][0],geo_dict[file]["geo"][count_i][1],geo_dict[file]["geo"][count_i][2],geo_dict[file]["atom_types"][count_i]))

# Wrapper function for write commands for *modes files
def write_modelist(name,modes):
    with open(name,'w') as f:
        for count_i,i in enumerate(modes.keys()):
            if len(i) == 2:
                modetype='bond'
            elif len(i) == 3:
                modetype='angle'
            elif len(i) == 4 and "R" in i[1] and "R" in i[2]:
                modetype='harmonic_dihedral'                
            elif len(i) == 4:
                modetype='dihedral'                
            f.write("\n{} start\n".format(modetype))
            f.write('{}\n'.format(len(modes[i]["modes"])))
            f.write("{}\n".format(" ".join([ "{:<60s}".format("_".join([ str(k) for k in j])) for j in modes[i]["modes"] ])))
            for count_j,j in enumerate(modes[i]["atomtypes"][0]):
                f.write('{}\n'.format(" ".join([ "{:<60s}".format(modes[i]["atomtypes"][k][count_j]) for k in range(len(modes[i]["atomtypes"])) ])))
            f.write("{} end\n".format(modetype))
    return 
        

# Returns the geometry and relevant property lists corresponding to the
# smallest fragment that is consistent with the mode being parametrized.
# If the 1 or 4 atom is a ring, then the entire ring is included. 
def mode_frag(M,Geometry,Adj_mat,Elements,Atom_types,gens,keep_terminal=False):

    # Initialize mass_dict (used for identifying the dihedral among a coincident set that will be explicitly scanned)
    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Create copies of input lists/arrays
    N_Elements = copy(Elements)
    N_Atom_types = deepcopy(Atom_types)
    N_Geometry = copy(Geometry)
    N_Adj_mat = copy(Adj_mat)

    if len(M) < 4:
        current = list(M)             # Holds the atoms found in each generation (seeded with the mode atoms)
        atoms = list(M)               # Holds the total set of included atoms (seeded with the mode atoms)
    else:
        current = [M[1],M[2]]         # Holds the atoms found in each generation (seeded with the 2-3 atoms)
        atoms = [M[1],M[2]]           # Holds the total set of included atoms (seeded with the 2-3 atoms)

    # Perform a depth search based on the gens argument
    for i in range(gens):

        # Initialize an empty list to hold the indices of atoms discovered in this generation
        new = []

        # Iterate over the atoms discovered in the previous generation (seeded initially with the bond atoms)
        for j in current:

            # Add bonded atoms to the new generation
            new += [ count_k for count_k,k in enumerate(N_Adj_mat[j]) if k == 1 and count_k not in atoms ]

        # Update the current and atoms lists with the atoms discovered in this generation
        # NOTE: the current list is overwritten while the atoms list is appended. The set call avoids double additions (sometimes occurs with cyclic topologies)
        current = list(set(new))
        atoms += list(set(new))

    # Finish any rings that were included in the initial gens-bond deep search
    rings = set([ i for i in new if "R" in Atom_types[i].split(']')[0].split('[')[0] ])
    while len(rings) > 0:
        atoms += [ i for i in rings if i not in atoms ]
        for i in rings:
            tmp = list(set([ count_j for count_j,j in enumerate(N_Adj_mat[i]) if j == 1 and "R" not in Atom_types[count_j].split(']')[0].split('[')[0] and count_j not in atoms ])) # Add non-ring atoms connected one deep from the ring 
            new += tmp
            atoms += tmp
        rings = set([ count_j for i in rings for count_j,j in enumerate(N_Adj_mat[i]) if j == 1 and "R" in Atom_types[count_j].split(']')[0].split('[')[0] and count_j not in atoms ]) # Add unique ring atoms

    # If keep_terminal is True, then for the last generation, hydrogenate non-terminal bonds of heavy atoms, terminal bonded atoms are kept as is
    if keep_terminal is True:
        for i in new:
            ind = [ count_k for count_k,k in enumerate(N_Adj_mat[i]) if k == 1 and count_k not in atoms ]
            term = [ k for k in ind if sum(N_Adj_mat[k]) == 1 ]
            atoms += ind
            for j in ind:

                # Replace non-terminal link atoms with hydrogens
                if j not in term: 
                    N_Elements[j] = "H"
                    N_Geometry[j] = N_Geometry[i] + ( N_Geometry[j] - N_Geometry[i] ) / norm( N_Geometry[j] - N_Geometry[i] ) * 1.1

                # Assign link atom type. In the FF, 0 is reserved for dummy atoms. Atom types are prepended with a "-link" so that the 
                # the script for extracting parameters can use the proper charges/vdw for these types.
                N_Atom_types[j] = 'link-'+N_Atom_types[j]

        # Canonicalize by sorting the elements based on hashing
        Masses = [ mass_dict[N_Elements[i]] for i in atoms ]
        A = N_Adj_mat[atoms][:,atoms]
        hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_i,A,Masses),i) for count_i,i in enumerate(atoms) ],reverse=True)) ]

        # Update lists/arrays based on atoms
        N_M          = tuple([ atoms.index(i) for i in M ])
        N_Geometry   = N_Geometry[atoms]
        N_Adj_mat    = N_Adj_mat[atoms]
        N_Adj_mat    = N_Adj_mat[:,atoms]
        N_Elements   = [ N_Elements[i] for i in atoms ]
        N_Atom_types = [ N_Atom_types[i] for i in atoms ]

    # Perform full hydrogenation of the fragment.
    else:        

        # Update lists/arrays based on atoms
        N_M          = tuple([ atoms.index(i) for i in M ])
        N_Geometry   = N_Geometry[atoms]
        N_Adj_mat    = N_Adj_mat[atoms]
        N_Adj_mat    = N_Adj_mat[:,atoms]
        N_Elements   = [ N_Elements[i] for i in atoms ]
        N_Atom_types = [ N_Atom_types[i] for i in atoms ]

        # Perform Hydrogenation
        N_Geometry,tmp_Atom_types,N_Elements,N_Adj_mat,added_idx = add_hydrogens(N_Geometry,N_Adj_mat,deepcopy(N_Atom_types),N_Elements)

        # Update the link types
        N_Atom_types += [ 'link-'+tmp_Atom_types[i] for i in added_idx ]
        
        # Canonicalize by sorting the elements based on hashing (NOTE: range(len(N_atom_types)) is used here rather than "atoms" as in the keep_terminal is True option. 
        Masses = [ mass_dict[N_Elements[i]] for i in range(len(N_Atom_types)) ]
        hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_i,N_Adj_mat,Masses),i) for count_i,i in enumerate(range(len(N_Atom_types))) ],reverse=True)) ]
        
        # Update lists/arrays based on the sorted atoms
        N_M          = tuple([ atoms.index(i) for i in N_M ])
        N_Geometry   = N_Geometry[atoms]
        N_Adj_mat    = N_Adj_mat[atoms]
        N_Adj_mat    = N_Adj_mat[:,atoms]
        N_Elements   = [ N_Elements[i] for i in atoms ]
        N_Atom_types = [ N_Atom_types[i] for i in atoms ]

    return N_M,N_Geometry,N_Adj_mat,N_Elements,N_Atom_types,hash_list

# Description: This is a simple implementation of the Dijkstra algorithm for 
#              finding the backbone of a polymer 
def Dijkstra(Adj_mat):

    # Remove terminal sites (sites with only a single length 2
    # self walk). Continue until all of the terminal structure 
    # has been removed from the topology.
    Adj_trimmed = copy(Adj_mat)

    # Initialize Distances, Previous, and Visited lists    
    Distances = array([0] + [100000]*(len(Adj_mat)-1)) # Holds shortest distance to origin from each site
    Previous = array([0]*len(Adj_mat)) # Holds the previous site on the short distance to origin
    Visited = [0]*len(Adj_mat) # Holds which sites have been visited

    # Initialize current site and neighbors list
    i = 0 # current site
    neighbors = []

    # Iterate through sites. At each step find the shortest distance between all of hte 
    # current sites neighbors and the START. Update the shortest distances of all sites
    # and choose the next site to iterate on based on which has the shortest distance of
    # among the UNVISITED set of sites. Once the terminal site is identified the shortest
    # path has been found
    while( 0 in Visited):

        # If the current site is the terminal site, then the algorithm is finished
        if i == len(Adj_trimmed)-1:
            break

        # Add new neighbors to the list
        neighbors = [ count_j for count_j,j in enumerate(Adj_trimmed[i]) if j == 1 ]

        # Remove the current site from the list of unvisited
        Visited[i] = 1

        # Iterate over neighbors and update shortest paths
        for j in neighbors:

            # Update distances for current generation of connections
            if Distances[i] + Adj_trimmed[j,i] < Distances[j]:
                Distances[j] = Distances[i] + Adj_trimmed[j,i]
                Previous[j] = i
        
        # Find new site 
        tmp = min([ j for count_j,j in enumerate(Distances) if Visited[count_j] == 0])
        i = [ count_j for count_j,j in enumerate(Distances) if j == tmp and Visited[count_j] == 0 ][0]
        
    # Find shortest path by iterating backwards
    # starting with the end site.
    Shortest_path = [len(Adj_trimmed)-1]    
    i=len(Adj_trimmed)-1

    while( i != 0):
        Shortest_path = Shortest_path + [Previous[i]]    
        i = Previous[i]

    # Reverse order of the list to go from start to finish
    Shortest_path = Shortest_path[::-1]
    return Shortest_path

def write_xyz(Output,Elements,Geometry,Atom_types):
    
    # Open file for writing and write header
    fid = open(Output+'.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {}\n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],Atom_types[count_i]))

    fid.close()

def write_pdb(Output,Elements,Geometry,Atom_types):

    # Open file for writing and write header
    f = open(Output+'.pdb','w')
    f.write('COMPND    SYSTEM\nAUTHOR    Generated by polygen.py\n')
    count=0
    for count_i,i in enumerate(Elements):
        count=count+1
        f.write('HETATM {:>4} {:<3}  {:>3} 1       {:> 8.3f}{:> 8.3f}{:> 8.3f}  1.00 {: >6}           {:<3}\n'\
                .format(count,i," ",float(Geometry[count_i,0]),float(Geometry[count_i,1]),float(Geometry[count_i,2]),count_i,Atom_types[count_i]))     
    f.close()

# Adds parameters from the FF file(s) to the FF_dict.
def parse_FF_params(FF_files,FF_dict={"masses":{},"charges":{},"bonds":{},"angles":{},"dihedrals":{},"dihedrals_harmonic":{},"vdw":{}}):
                   
    modes_from_FF = []
    for i in FF_files:
        with open(i,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) == 0: continue
                if fields[0].lower() == "atom":   FF_dict["masses"][fields[1]] = float(fields[3])
                if fields[0].lower() == "charge": FF_dict["charges"][fields[1]] = float(fields[2])
                if fields[0].lower() == "bond":   
                    modes_from_FF += [(fields[1],fields[2])]
                    modes_from_FF += [(fields[2],fields[1])]
                    FF_dict["bonds"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                if fields[0].lower() == "angle":
                    modes_from_FF += [(fields[1],fields[2],fields[3])]
                    modes_from_FF += [(fields[3],fields[2],fields[1])]
                    FF_dict["angles"][(fields[1],fields[2],fields[3])] = [fields[4],float(fields[5]),float(fields[6])]
                if fields[0].lower() in ["dihedral","torsion"]: 
                    modes_from_FF += [(fields[1],fields[2],fields[3],fields[4])]
                    modes_from_FF += [(fields[4],fields[3],fields[2],fields[1])]
                    if fields[5] == "opls":       
                        FF_dict["dihedrals"][(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
                    elif fields[5] == "harmonic":
                        FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ] 
                    elif fields[5] == "quadratic":
                        FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),float(fields[7]) ]
                if fields[0].lower() == "vdw":    
                    FF_dict["vdw"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                    FF_dict["vdw"][(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]

    return FF_dict,modes_from_FF

# # This function handles remapping the atomtypes of one fragment on to another such that the modes (bonds
# # angles, dihedrals, etc) in the remapped fragment are the same after remapping. 
def map_frag(fixed_atomtypes,fixed_adj_mat,fixed_hash,remapped_atomtypes,remapped_adj_mat,remapped_hash):

    # The equivalency of the two topologies is based upon the mode types and number each topology possesses
    bond_types_0,angle_types_0,dihedral_types_0,one_five_types_0 = Find_modes(remapped_adj_mat,remapped_atomtypes,return_all=1)[4:]
    comp_obj = ( sorted(bond_types_0),sorted(angle_types_0),sorted(dihedral_types_0) )

    # Seed the mapping by placing the first atom in the mappable set onto the first instance 
    # of its hash type in the fixed set.
    #first = len(remapped_hash)-1
    first = 0
    mapping   = {first:fixed_hash.index(remapped_hash[first])}    
    R_mapping = {fixed_hash.index(remapped_hash[first]):first}    
    remapped_place_me = set(range(len(remapped_atomtypes)))
    fixed_place_me = set(range(len(remapped_atomtypes)))
    remapped_place_me.remove(first)
    fixed_place_me.remove(mapping[first])

    # Seed the connections list (idx in fixed, connection restraint in remapped) and iterate over the connections until no new connections are found
    cons=[ (count_i,first) for count_i,i in enumerate(fixed_adj_mat[mapping[first]]) if i == 1 and count_i in fixed_place_me ]
    for i in cons:        
        if i[0] not in fixed_place_me: continue
        remapped_match = next( j for j in remapped_place_me if remapped_hash[j] == fixed_hash[i[0]] and remapped_adj_mat[i[1],j] == 1 )
        mapping[remapped_match] = i[0]
        R_mapping[i[0]] = remapped_match
        fixed_place_me.remove(i[0])
        remapped_place_me.remove(remapped_match)
        cons += [ (count_j,remapped_match) for count_j,j in enumerate(fixed_adj_mat[i[0]]) if j == 1 and count_j in fixed_place_me ] 

    # Remap the atomtypes
    new_atomtypes = [ "X" ]*len(remapped_atomtypes)
    for i in mapping: new_atomtypes[mapping[i]] = remapped_atomtypes[i]

    # Check if the set of types and number of each mode match the original fragment
    # NOTE: the fixed_adj_mat is used for determining the connectivity
    bond_types_1,angle_types_1,dihedral_types_1,one_five_types_1 = Find_modes(fixed_adj_mat,new_atomtypes,return_all=1)[4:]
    if comp_obj != ( sorted(bond_types_1),sorted(angle_types_1),sorted(dihedral_types_1)):
        print "ERROR in map_frag: No mapping was discovered that preserved the topology of the original fragment. Check that the two fragments are indeed identical. Exiting..."
        print "{:60s} {:60s}".format("Original","New")
        for count_i,i in enumerate(remapped_atomtypes):
            print "{:<3d} {:60s} {:60s}".format(count_i,i,new_atomtypes[count_i])
        return new_atomtypes,1

    return new_atomtypes,mapping,0        

# # OLD FUNCTION: This function handles remapping the atomtypes of one fragment on to another such that the modes (bonds 
# # angles, dihedrals, etc) in the remapped fragment are the same after remapping. 
# def map_frag(fixed_atomtypes,fixed_adj_mat,remapped_atomtypes,remapped_adj_mat):

#     # A copy of the orginal list of atomtypes (and the original order) is needed for the current algorithm
#     remapped_atomtypes_0 = deepcopy(remapped_atomtypes)

#     # The set of allowed moves is based on the set of link atoms in both the fixed and remapped atomtypes (all link- atoms should be reassignable)
#     movable_ind = [ count_i for count_i,i in enumerate(remapped_atomtypes) if "link-" in i ]
# #    movable_ind = sorted(set( movable_ind + [ count_i for count_i,i in enumerate(fixed_atomtypes) if "link-" in i ] ))

#     # The equivalency of the two topologies is based upon the mode types and number each topology possesses
#     bond_types_0,angle_types_0,dihedral_types_0,one_five_types_0 = Find_modes(remapped_adj_mat,remapped_atomtypes,return_all=1)[4:]
#     comp_obj = ( sorted(bond_types_0),sorted(angle_types_0),sorted(dihedral_types_0) )

#     # The set of allowed moves is based on the set of atoms within "gens" bonds of any link atom in the mappable structure.
# #     gens = 3
# #     seps = graph_seps(remapped_adj_mat)
# #     links_remapped = [ count_i for count_i,i in enumerate(remapped_atomtypes) if "link-" in i ]
# #     movable_ind    = [ count_j for i in links_remapped for count_j,j in enumerate(seps[i]) if j <= gens ]
# #     movable_ind    = sorted(set(movable_ind))
# #    movable_types,movable_ind = zip(*sorted(zip([ remapped_atomtypes[i].split('[')[1].split(']')[0] if 'link-' not in remapped_atomtypes[i]\
# #    else remapped_atomtypes[i].split('link-')[1].split('[')[1].split(']')[0] for i in movable_ind ],movable_ind)))

#     # A flag is used to check if the function fails to identify an equivalent topology
#     flag = 0

#     # permutations is a generator for all unique orderings of the movable indices. 
#     # Thus, if there are four movable indices then permuations will return 4-element tuples comprised of the unique ordering of 0, 1, 2, and 3 
#     # for i in gen_perms(movable_types):
#     for i in permutations(range(len(movable_ind)),len(movable_ind)):

#         # The movable atoms are remapped based on the current permuation. 
#         # NOTE: movable_ind is a list of the atomic movable indices *in the atomtypes lists* 
#         # NOTE: the left side of the equality refers to the original index (count_j) while the right-side is based on the shuffled index (j)
#         for count_j,j in enumerate(i): remapped_atomtypes[movable_ind[count_j]] = remapped_atomtypes_0[movable_ind[j]]

#         # Check if the set of types and number of each mode match the original fragment
#         # NOTE: the fixed_adj_mat is used for determining the connectivity
#         bond_types_1,angle_types_1,dihedral_types_1,one_five_types_1 = Find_modes(fixed_adj_mat,remapped_atomtypes,return_all=1)[4:]
#         if comp_obj == ( sorted(bond_types_1),sorted(angle_types_1),sorted(dihedral_types_1)):
#             flag = 1
#             break

#     # Check for failure to find a matching fragment (this shouldn't happen if a suitable fragment was supplied)
#     if flag == 0:
#         print "ERROR in map_frag: No mapping was discovered that preserved the topology of the original fragment. Check that the two fragments are indeed identical. Exiting..."
#         print "{:60s} {:60s}".format("Original","New")
#         print "movable_ind: {}".format(movable_ind)
#         for count_i,i in enumerate(remapped_atomtypes):
#             print "{:<3d} {:60s} {:60s}".format(count_i,fixed_atomtypes[count_i],i)
#         return remapped_atomtypes_0,1

#     # Successful return
#     else:
#         return remapped_atomtypes,0

# Add hydrogens based upon the supplied atom types. 
# This function is only compatible with TAFFI atom types
# NOTE: Hydrogenation heuristics assume carbon behavior. This isn't usually a problem when the results are refined with transify, but more specific rules should be implemented in the future
def add_hydrogens(geo,adj_mat,atomtypes,elements,saturate=True,retype=True):

    # Initialize the saturation dictionary the first time this function is called
    if not hasattr(add_hydrogens, "sat_dict"):
        add_hydrogens.sat_dict = {  'H':2, 'He':1,\
                                   'Li':1, 'Be':2,                                                                                                                'B':3,     'C':4,     'N':3,     'O':2,     'F':1,    'Ne':1,\
                                   'Na':1, 'Mg':2,                                                                                                               'Al':3,    'Si':4,     'P':3,     'S':2,    'Cl':1,    'Ar':1,\
                                    'K':1, 'Ca':2, 'Sc':None, 'Ti':None,  'V':None, 'Cr':None, 'Mn':None, 'Fe':None, 'Co':None, 'Ni':None, 'Cu':None, 'Zn':None, 'Ga':None, 'Ge':None, 'As':None, 'Se':None, 'Br':1,    'Kr':None,\
                                   'Rb':1, 'Sr':2,  'Y':None, 'Zr':None, 'Nb':None, 'Mo':None, 'Tc':None, 'Ru':None, 'Rh':None, 'Pd':None, 'Ag':None, 'Cd':None, 'In':None, 'Sn':None, 'Sb':None, 'Te':None,  'I':1,    'Xe':None,\
                                   'Cs':1, 'Ba':2, 'La':None, 'Hf':None, 'Ta':None,  'W':None, 'Re':None, 'Os':None, 'Ir':None, 'Pt':None, 'Au':None, 'Hg':None, 'Tl':None, 'Pb':None, 'Bi':None, 'Po':None, 'At':None, 'Rn':None  }

    H_length = 1.1
    N_atoms  = len(geo)
    init_len = len(geo)

    # Loop over the atoms in the geometry
    for count_i,i in enumerate(geo):
        
        # ID undercoordinated atoms
        N_expected = sum(type_adjmat(atomtypes[count_i])[0][0])
        N_current  = sum(adj_mat[count_i])

        # Add hydrogens to undercoordinated atoms
        if N_expected > N_current: 

            # If the saturate option is True, then the atom is fully hydrogenated to be consistent with the corresponding number of bonds in add_hydrogens.sat_dict
            if saturate is True and add_hydrogens.sat_dict[elements[count_i]] is not None and N_expected < add_hydrogens.sat_dict[elements[count_i]]:
                N_expected = add_hydrogens.sat_dict[elements[count_i]]

            old_inds = [ count_j for count_j,j in enumerate(adj_mat[count_i]) if j == 1 ]

            # Protocols for 1 missing hydrogen
            if N_expected - N_current == 1:
                if N_expected == 1:
                    new = i + array([H_length,0.0,0.0])
                elif N_expected == 2:
                    new = -1.0 * normalize(geo[old_inds[0]] - i) * H_length + i + array([random.random(),random.random(),random.random()])*0.01 #random factor added for non-carbon types to relax during FF-opt
                elif N_expected == 3:
                    new = -1.0 * normalize( normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) ) * H_length + i
                elif N_expected == 4:
                    print "WARNING in add_hydrogens: 1:4 (unchecked)"
                    new = -1.0 * normalize( normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) + normalize(geo[old_inds[2]] - i) ) * H_length + i                

                # Update geometry, adj_mat, elements, and atomtypes with one new atoms
                geo = vstack([geo,new])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]
                elements += ["H"]
                tmp = zeros([N_atoms+1,N_atoms+1])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[-1,count_i] = 1
                tmp[count_i,-1] = 1
                adj_mat = tmp                
                N_atoms += 1

            # Protocols for 2 missing hydrogens
            if N_expected - N_current == 2:
                if N_expected == 2:
                    new_1 = i + array([H_length,0.0,0.0])
                    new_2 = i - array([H_length,0.0,0.0])
                elif N_expected == 3:
                    rot_vec = normalize(cross( geo[old_inds[0]] - i, array([random.random(),random.random(),random.random()]) ))
                    new_1 = normalize(axis_rot(geo[old_inds[0]],rot_vec,i,120.0) - i)*H_length + i
                    new_2 = normalize(axis_rot(geo[old_inds[0]],rot_vec,i,240.0) - i)*H_length + i
                elif N_expected == 4:
                    bisector = normalize(geo[old_inds[0]] - i + geo[old_inds[1]] - i) 
                    new_1    = axis_rot(geo[old_inds[0]],bisector,i,90.0)
                    new_2    = axis_rot(geo[old_inds[1]],bisector,i,90.0) 
                    rot_vec  = normalize(cross(new_1-i,new_2-i))
                    angle    = ( 109.5 - acos(dot(normalize(new_1-i),normalize(new_2-i)))*180.0/pi ) / 2.0
                    new_1    = axis_rot(new_1,rot_vec,i,-angle)
                    new_2    = axis_rot(new_2,rot_vec,i,angle)
                    new_1    = -1*H_length*normalize(new_1-i) + i
                    new_2    = -1*H_length*normalize(new_2-i) + i
                    
                # Update geometry, adj_mat, elements, and atomtypes with two new atoms
                geo = vstack([geo,new_1])
                geo = vstack([geo,new_2])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*2
                elements += ["H","H"]
                tmp = zeros([N_atoms+2,N_atoms+2])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2],count_i] = 1
                tmp[count_i,[-1,-2]] = 1
                adj_mat = tmp
                N_atoms += 2

            # Protocols for 3 missing hydrogens
            if N_expected - N_current == 3:
                if N_expected == 3:
                    rot_vec = array([0.0,1.0,0.0])
                    new_1 = i + array([H_length,0.0,0.0])
                    new_2 = axis_rot(new_1,rot_vec,i,120.0)
                    new_3 = axis_rot(new_1,rot_vec,i,240.0)
                if N_expected == 4:
                    rot_vec = normalize(cross( geo[old_inds[0]] - i, array([random.random(),random.random(),random.random()]) ))
                    new_1 = H_length*normalize(axis_rot(geo[old_inds[0]],rot_vec,i,109.5)-i) + i
                    new_2 = axis_rot(new_1,normalize(i-geo[old_inds[0]]),i,120.0)
                    new_3 = axis_rot(new_2,normalize(i-geo[old_inds[0]]),i,120.0)

                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = vstack([geo,new_1])
                geo = vstack([geo,new_2])
                geo = vstack([geo,new_3])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*3
                elements += ["H","H","H"]
                tmp = zeros([N_atoms+3,N_atoms+3])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2,-3],count_i] = 1
                tmp[count_i,[-1,-2,-3]] = 1
                adj_mat = tmp
                N_atoms += 3

            # Protocols for 4 missing hydrogens
            if N_expected - N_current == 4:
                if N_expected == 4:
                    new_1 = i + array([H_length,0.0,0.0])
                    rot_vec = normalize(cross( new_1 - i, array([random.random(),random.random(),random.random()]) ))
                    new_2 = H_length*normalize(axis_rot(new_1,rot_vec,i,109.5)-i) + i
                    new_3 = axis_rot(new_2,normalize(i-new_1),i,120.0)
                    new_4 = axis_rot(new_3,normalize(i-new_1),i,120.0)
                    
                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = vstack([geo,new_1])
                geo = vstack([geo,new_2])
                geo = vstack([geo,new_3])
                geo = vstack([geo,new_4])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*4
                elements += ["H","H","H","H"]
                tmp = zeros([N_atoms+4,N_atoms+4])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2,-3,-4],count_i] = 1
                tmp[count_i,[-1,-2,-3,-4]] = 1
                adj_mat = tmp
                N_atoms += 4

    if retype is True:
        #return geo,id_types(elements,adj_mat,gens=2,geo=geo),elements,adj_mat,range(init_len,len(geo))
        return geo,id_types(elements,adj_mat,gens=2),elements,adj_mat,range(init_len,len(geo)) # id_types v.062520
    else:
        return geo,atomtypes,elements,adj_mat,range(init_len,len(geo))

# Shortcut for normalizing a vector
def normalize(x):
    return x/sum(x**(2.0))**(0.5)

# Recursively generates the permutations of indices corresponding to identical values in types
def gen_perms(types,start=0):
    ind = [ count_i for count_i,i in enumerate(types) if i==types[0] ]
    not_ind = [ count_i for count_i,i in enumerate(types) if i !=types[0] ]
    for i in permutations(range(start,start+len(ind)),len(ind)):
        if len(not_ind) > 0:
            for j in gen_perms([ types[j] for j in not_ind ],start=len(ind)+start):
                yield i+j
        else:
            yield i

# Logger object redirects standard output to a file.
class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/frag_gen.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
            
if __name__ == "__main__":
   main(sys.argv[1:])


