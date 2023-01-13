

import sys,os
from numpy import *
from numpy.linalg import norm

# Add TAFFY Lib to path
from adjacency import *
from write_functions import *
from copy import deepcopy
import random
import itertools

# A simple wrapper script that calculates the mappings for a molecule
def main(argv):

    if len(argv) == 1:
        N_beads = 2
    else:
        N_beads = int(argv[1])

    # Check if the mappings directory exists
    if os.path.isdir("mappings") is True:
        print "mappings already exists, exiting..."
        quit()
    else:
        os.mkdir("mappings")

    # Read in the geometry
    count = 0
    with open(argv[0],'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            if lc == 0:
                geometry = zeros([int(fields[0]),3])
                elements = []
                charges  = zeros(int(fields[0]))
                continue
            if len(fields) == 0:
                continue
            else:
                geometry[count] = array([float(fields[1]),float(fields[2]),float(fields[3])])
                elements       += [fields[0]]
                charges[count]  = float(fields[4])
                count += 1

    # neutralize the monomer
    q_offset = -sum(charges)/float(len(elements))
    for i in range(len(geometry)):
        charges[i] += q_offset
        #print "{} {} {}".format(elements[i],geometry[i],charges[i])

    # Center the geometry about the center of mass
    cm = calc_cm(elements,geometry)    
    geometry = geometry - cm
    xyz_write("mappings/{}".format("geometry.xyz"),elements,geometry)

    # calculate adjmat
    d_0,d_mag_0 = calc_dipole(geometry,charges)

    # Generate the adjacency matrix
    # NOTE: the units are converted back angstroms
    adj_mat = Table_generator(elements,geometry)

    # Generate all of the 2-bead mappings
    print "generating {}-bead mappings...".format(N_beads)
    mappings = split_graph(adj_mat,N_subs=N_beads)

    # Calculate the dipole of each mapping
    d = []
    d_mags = []
    errors = []
    errors_2 = []
    for count_i,i in enumerate(mappings):
        geo_m = array([ calc_cm([ elements[_] for _ in j ],geometry[j,:]) for j in i ])
        q_m   = array([sum(charges[j]) for j in i ])        
        d_tmp,d_mag_tmp = calc_dipole(geo_m,q_m)
        scale = d_mag_0 / d_mag_tmp
        q_m   = q_m * scale
        d_tmp,d_mag_tmp = calc_dipole(geo_m,q_m)
        d += [d_tmp]        
        d_mags += [d_mag_tmp]
        errors += [ sum(abs(d_0-d_tmp)) ]
        errors_2 += [ sum((d_0-d_tmp)**(2.0)) ]
        xyz_write("mappings/{}.xyz".format(str(count_i)),["X","X"],geo_m)

    errors_2,d,d_mags,mappings,errors = zip(*sorted([ (i,d[count_i],d_mags[count_i],mappings[count_i],errors[count_i]) for count_i,i in enumerate(errors_2) ]))

    print "Parse complete, {} total mappings found...".format(N_beads,len(mappings))
    print "Ranking representaions based on dipole moment...\n"
    print "    Target dipole: {}".format(d_0)
    print "    Target |d|:    {}\n".format(d_mag_0)

    print "    {:<20s} {:<80s} {:<20s} {:<20s}".format("Rank","Mapping","sum(abs(d_0-d_map))","sum((d_0-d_map)**(2))")
    for i in range(len(d)):
        print "    {:<20d} {:<80s} {:<20.8f} {:<20.8f}".format(i,mappings[i],errors[i],errors_2[i])
        
    return

# Description: Calculates the center of mass for the molecule/collection of atoms
#              supplied to the geo argument. 
#
# Inputs:      elements (list): Contains a list of strings holding the elemental labels of each atom in geo.
#              geo (Nx3 array): Contains the cartesian coordinates of each atom. The rows are indexed to elements. 
#
# Returns:     cm (1x3 array):  Contains the cartesian center of mass of geo. 
def calc_cm(elements,geo):

    # Define the dictionary of masses on the first call of the function
    if hasattr(calc_cm,'Masses') is False:
        
        calc_cm.Masses = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                          'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                          'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                          'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                          'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                          'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                          'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                          'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Check for missing elements
    if set(elements).issubset(calc_cm.Masses.keys()) == False:
        print "ERROR in calc_cm: Encountered an element that wasn't in the masses dictionary while trying to calculate the center of mass. Exiting..."
        quit()

    # Calculate the center of mass
    else:
        masses = array([ calc_cm.Masses[i] for i in elements ])
        cm = sum(masses[:,None] * geo,axis=0)/sum(masses)

    # Return result
    return cm

# Description: Returns a unique mapping key for use in redundancy checks. The mapping is canonicalized by
#              sorting the sublists with respect to the lowest index atom in each sublist. 
#
# Inputs:      mapping (list of lists): contains the CG mapping as a list of lists, where
#                                       each sublist corresponds to a bead.
#
# Returns:     (str): contains a string representation of the canonical mapping expression.
def return_mapping_key(mapping):
    _,tmp = zip(*sorted([ (sorted(_)[0],sorted(_)) for _ in mapping  ]))
    return str(tmp)    

# Description: Returns all unique ways of splitting a graph into two subgraphs. This is a more complicated
#              problem than just a "choose N" algorithm, since the subgraphs by definition must be connected
#              components. The algorithm here starts by seeding a subgraph search with each (single) node in
#              the adjacency matrix. If the remaining atoms at connected, then this is a valid mapping (the
#              ith atom in one subgraph, and the rest in the other). The atoms that are connected to i are
#              added to i, and used as a seed for the next search. Unique mappings are hashed to avoid 
#              redundant attempts to split the graph. 
#
# Inputs:      adj_mat:  (array-type) the adjacency matrix holding the graph to be split
#              inds:     (list type) holds the indices of the adjacency matrix. This can be useful for
#                        recursive calls on subviews of the adjacency matrix, where the actual indices 
#                        are difference from their index in adj_mat. 
#
# Returns:     mappings: (list of lists) a list of unique mappings. Each mapping is itself a list of two
#                        lists, correpsonding to the indices corresponding to subgraph "a" and "b". For
#                        example, mappings[5][0] would return list of nodes in subgraph "a" of the sixth
#                        representation. The nodes in each subgraph are sorted in ascending order according
#                        to index, the "a" subgraph is determined by which possesses the lowest index node.  
#
def split_graph(adj_mat,inds=None,N_subs=2):

    # Initialize the inds variable
    # else: check the inds list for consistency
    if inds is None:        
        inds   = range(len(adj_mat)) 
    else:
        if len(inds) != len(adj_mat):
            print "ERROR in split_graph: inds variable size does not match adj_mat"
            quit()

    # If there are less than 2 nodes in the graph then there is no way to split it
    if len(adj_mat) < 2:
        return []

    # If one subgraph is requested then the indices of all nodes are simply returned
    if N_subs == 1:
        return [sorted(inds)]

    # Initialize necessary lists and sets
    a_inds   = range(len(adj_mat))
    p_groups = [ [_] for _ in range(len(adj_mat)) ]        # Seed the subgraph search with the individual nodes
    mappings = []                                          # list of the mappings
    kept_keys = set([])                                    # set of the retained mappings 
    attempted_keys = set([])                               # set of the attempted mappings
 
    # Loop over the "potential" subgraphs in "p_groups"
    for i in p_groups:

        # Avoid the case where all nodes are in i (others is empty)
        others = [ _ for _ in a_inds if  _ not in i ]
        if len(others) == 0:
            continue

        # Avoid repeating the test if this representation has already been attempted
        # If continue, else add the current key to the attempted_keys set
        key = return_mapping_key([sorted(i),sorted(others)])
        if key in attempted_keys:
            continue
        else:
            attempted_keys.add(key)
        
        # If the remaining set is still connected (one subgraph) then create the mapping with i,others and add its key to the         
        # kept_keys set. find_subgraphs function returns all connected subgraphs, in this case we only care about the number. 
        if len(find_subgraphs(adj_mat[others,:][:,others])) == 1:
            if key not in kept_keys:
                mappings += [[sorted(i),sorted(others)]]
                kept_keys.add(key)

        # Determine connections with group i atoms
        connections  = set([ count_k for j in i for count_k,k in enumerate(adj_mat[j]) if k == 1 and count_k not in i ])

        # Create the new groups for testing
        for j in connections: p_groups += [ [j]+i ]

    # If N_subs == 2 then only a single split is requested and recursive calls are not required
    if N_subs == 2:
        mappings_final = mappings

        # Converts the final mappings to the absolute indices (based on inds list)
        for count_i,i in enumerate(mappings_final):     
            mappings_final[count_i] = [ [ inds[_] for _ in j ] for j in i ]  

    # Recursive call to handle further decomposition of the graph
    else:
        mappings_final = []
        for count_i,i in enumerate(mappings):

            # Recursive call to split the i[1] (i[0]) subgraph N_subs - 1 additional times. The result is iterated over and added to 
            # combined with the i[1] (i[0]) subgraph to generate valid mappings. The "j" mappings are returned as absolute indicies
            # whereas the i subgraphs need to be remapped
            for j in split_graph(adj_mat[i[0],:][:,i[0]],inds=[ inds[_] for _ in i[0] ],N_subs=N_subs-1):
                mappings_final += [j+[[ inds[_] for _ in i[1] ]]] 
            for j in split_graph(adj_mat[i[1],:][:,i[1]],inds=[ inds[_] for _ in i[1] ],N_subs=N_subs-1):
                mappings_final += [[[ inds[_] for _ in i[0] ]]+j] 

    # Sort the final mappings
    for count_i,i in enumerate(mappings_final):
        mappings_final[count_i] = sorted(mappings_final[count_i])

    # For N>2 subdivisions there is the possibility of redundant mappings. 
    # Redundancies are removed here on the basis of unique mapping keys.
    if N_subs > 2:
        kept_keys = set([])               # set of the retained mappings 
        mappings_red = []                 # non-redundant mappings
        for i in mappings_final:          # loop over the mappings 
            key = return_mapping_key(i)   # determine the mapping key for the ith mapping
            if key in kept_keys:          # skip if the key has already been stored
                continue         
            else:                         
                mappings_red += [i]       # add non-redundant mappping
                kept_keys.add(key)        # add the key to the kept_keys set
        mappings_final = mappings_red     # update mappings_final with the non-redundant set of mappings
            
    return mappings_final

# Description: Returns the subgraphs of an adjacency matrix. The sugraphs are returned as a list of lists, 
#              where each sublist holds the indices of the nodes in each subgraph.
#
# Inputs:      adj_mat (array NxN): Holds the adjacency matrix for the graph whose subgraph structure is being characterized
#
# Returns:     sugraph_list (list of lists): Holds the indices of nodes in the subgraph. The indices are returned based on the 
#                                            values supplied in idx, or by default, their index in the supplied adjacency matrix.
#
def find_subgraphs(adj_mat,idx=None):

    # Initialize list of indices
    if idx is None:
        idx = range(len(adj_mat))

    # Find subgraphs
    subgraph_list = []
    placed_list   = []
    for count_i,i in enumerate(adj_mat):

        # skip if the current node has already been placed in a subgraph
        if count_i in placed_list: continue

        # current holds the subgraph seedlings
        subgraph = [count_i] + [ count_j for count_j,j in enumerate(i) if j == 1 ]

        # new holds the updated list of nodes in the subgraph at the end of each bond search
        new = [ count_j for count_j,j in enumerate(i) if j == 1 ]

        # recursively add nodes to new until no new nodes are found
        while new != []:

            # Initialize list to hold new connnected nodes
            connections = []

            # Iterate over the new nodes and find new connections
            for count_j,j in enumerate(new):
                connections += [ count_k for count_k,k in enumerate(adj_mat[j]) if k == 1 and count_k not in subgraph ] # add new connections not already in the subgraph

            # Update lists
            connections = list(set(connections))
            subgraph += connections
            new = connections

        # Sort nodes in the subgraph (makes them unique)
        subgraph = sorted(subgraph)

        # Add new subgraphs to the master list and add its nodes to the placed_list to avoid redundant searches
        if subgraph not in subgraph_list:
            subgraph_list += [subgraph]
            placed_list += subgraph

    # Relable the subgraphs with the absolute indices supplied by the idx list
    for i in range(len(subgraph_list)):
        subgraph_list[i] = [ idx[j] for j in subgraph_list[i] ]
    
    return subgraph_list

# Description: Calculates the dipole of a collection of point charges.
#              Note for charged molecules the origin for the calculation matters. 
#
# input        geo        : Nx3 array
#              charges    : Nx1 array
#
# returns      dipole     : 1x3 array
#              dipole mag : scalar
#
def calc_dipole(geo,charges,r_0=None):

    # If an origin is supplied, then perform a consistency check and center the geometry about it.
    if r_0 is not None:
        
        # Consistency check
        if len(r_0) != 3:
            print "ERROR in calc_dipole: r_0 must be a 3 element array. Exiting..."
            quit()

        # Center the geometry about the supplied origin
        geo = geo - r_0

    # Center about the centroid by default
    else:
        centroid = mean(geo,axis=0)
        geo = geo - centroid

    # Return the dipole moment and magnitude
    dipole =  sum(geo*charges[:,None],axis=0)
    return dipole,norm(dipole)

# main sentinel
if __name__ == "__main__":
    main(sys.argv[1:])
