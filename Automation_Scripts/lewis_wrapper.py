#!/bin/env python
"""
Created on Thu Oct 12 15:24:44 2017

@author: Stephen Shiring
"""

#import sys,argparse,os,datetime,fnmatch,os,re
import sys,os,argparse
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')

from id_types import *
from file_parsers import *
from adjacency import *

def doLewis(filename):

    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                           'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                           'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                           'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                           'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                           'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                           'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                           'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}
    
    elements, geometry = xyz_parse(filename)
    adj_matrix = Table_generator(elements,geometry)
    atom_types = id_types(elements,adj_matrix)
    
    # Canonicalize by sorting the elements based on hashing (NOTE: range(len(N_atom_types)) is used here rather than "atoms" as in the keep_terminal is True option. 
    Masses = [ mass_dict[elements[i]] for i in range(len(atom_types)) ]
    hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_i,adj_matrix,Masses),i) for count_i,i in enumerate(range(len(atom_types))) ],reverse=True)) ]

    bonding_pref = []

    lone_electrons, bonding_electrons, core_electrons, adj_mat, total_formal_charge, bonding_pref = find_lewis(atom_types, adj_matrix, bonding_pref=bonding_pref, return_pref=True)
    
    return elements, geometry, adj_mat, atom_types, core_electrons, bonding_electrons, lone_electrons, total_formal_charge

# Function that returns the internal modes of a molecule                                                                                                                                                                                   
#                                                                                                                                                                                                                                          
# inputs:       A, Adjacency matrix (numpy array or list of lists)                                                                                                                                                                         
# returns:      dictionary, with lists of 2, 3, and 4-tuples keyed                                                                                                                                                                         
#               to "bonds", "angles", and "dihedrals" (all dihedrals,
#               plus flexible and inflexible dihedrals), respectively                                                                                                                                                                        
#               and corresponding to the same.                                                                                                                                                                                             
def find_modes(A):

    # List comprehension to determine bonds from a loop over the adjacency matrix. Iterates over rows (i) and individual elements                                                                                                          
    # ( elements A[count_i,count_j] = j ) and stores the bond if the element is "1". The count_i < count_j condition avoids                                                                                                                
    # redudant bonds (e.g., (i,j) vs (j,i) ). By convention only the i < j definition is stored.                                                                                                                                           
    bonds      = [ (count_i,count_j) for count_i,i in enumerate(A) for count_j,j in enumerate(i) if count_i < count_j and j != 0 ]

    # List comprehension to determine angles from a loop over the bonds. Note, since there are two bonds in every angle, there will be                                                                                                     
    # redundant angles stored (e.g., (i,j,k) vs (k,j,i) ). By convention only the i < k definition is stored.                                                                                                                              
    angles     = [ (count_j,i[0],i[1]) for i in bonds for count_j,j in enumerate(A[i[0]]) if j != 0 and count_j < i[1] and count_j not in i ]

    # List comprehension to determine dihedrals from a loop over the angles. Note, since there are two angles in every dihedral, there will be                                                                                                
    # redundant dihedrals stored (e.g., (i,j,k,m) vs (m,k,j,i) ). By convention only the i < m definition is stored.    
    # List all the dihedrals in the molecule                                                                                                                   
    dihedrals  = [ (count_j,i[0],i[1],i[2]) for i in angles for count_j,j in enumerate(A[i[0]]) if j != 0 and count_j < i[2] and count_j not in i ]
    
    # Identify the flexible (j-k connected by a single or triple bond) and inflexible (j-k connected by a double bond) dihedrals
    flexible  = [ (count_j,i[0],i[1],i[2]) for i in angles for count_j,j in enumerate(A[i[0]]) if j !=0 and A[i[0],i[1]] != 2 and count_j < i[2] and count_j not in i ]
    inflexible  = [ (count_j,i[0],i[1],i[2]) for i in angles for count_j,j in enumerate(A[i[0]]) if j !=0 and A[i[0],i[1]] == 2 and count_j < i[2] and count_j not in i ]

    return {"bonds":bonds,"angles":angles,"dihedrals":dihedrals,"flexible":flexible,"inflexible":inflexible}

def main(argv):
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument('xyz_list', help = '')

    args = parser.parse_args()
    xyz = args.xyz_list
    
    elements, geometry, adj_mat, atom_types, core_electrons, bonding_electrons, lone_electrons, total_formal_charge = doLewis(xyz)
    
    output = open(xyz[:xyz.find('.')]+'.lewis','w')
    output.write(str(xyz) + '\n')
    output.write('\nElements:\n' + str(elements) + '\n\nGeometry:\n' + str(geometry) + '\n')
    output.write('\nAdjaceny matrix:\n' + str(adj_mat) + '\n')
    output.write('\nAtom types:\n' + str(atom_types) + '\n')
    output.write('\nCore electrons:\n' + str(core_electrons) + '\n')
    output.write('\nBonding electrons:\n' + str(bonding_electrons) + '\n')
    output.write('\nLone electrons:\n' + str(lone_electrons) + '\n')
    output.write('\n\n')
    output.write('{:^10} {:^10} {:^10} {:^10}\n'.format('Element','core e','bonding e', 'lone e'))
    
    for i in range (0,len(elements)):
        #print str(elements[i]) + ' | ' + str(core_electrons[i]) + ' | ' + str(bonding_electrons[i]) + ' | ' + str(lone_electrons[i])
        output.write('{:^10} {:^10d} {:^10d} {:^10d}\n'.format(elements[i], core_electrons[i], bonding_electrons[i], lone_electrons[i]))
    
    output.write('\nmolecular charge: ' + str(total_formal_charge))
    output.write('\n\n')
    
    # Find the modes in A                                                                                                                                                                                                                      
    modes = find_modes(adj_mat)

    # Print results                                                                                                                                                                                                                        
    print("\nBonds discovered ({}):\n".format(len(modes["bonds"])))
    print("\t" + "\t".join([ str(_)+"\n" for _ in modes["bonds"]]))
    
    output.write("\nBonds discovered ({}):\n".format(len(modes["bonds"])))
    output.write("\t" + "\t".join([ str(_)+"\n" for _ in modes["bonds"]]))

    print("\nAngles discovered ({}):\n".format(len(modes["angles"])))
    print("\t" + "\t".join([ str(_)+"\n" for _ in modes["angles"]]))
    
    output.write("\nAngles discovered ({}):\n".format(len(modes["angles"])))
    output.write("\t" + "\t".join([ str(_)+"\n" for _ in modes["angles"]]))

    print("\nDihedrals discovered ({}):\n".format(len(modes["dihedrals"])))
    print("\t" + "\t".join([ str(_)+"\n" for _ in modes["dihedrals"]]))
    
    output.write("\nDihedrals discovered ({}):\n".format(len(modes["dihedrals"])))
    output.write("\t" + "\t".join([ str(_)+"\n" for _ in modes["dihedrals"]]))
    
    print("\nFlexible dihedrals discovered ({}):\n".format(len(modes["flexible"])))
    print("\t" + "\t".join([ str(_)+"\n" for _ in modes["flexible"]]))
    
    output.write("\nFlexible dihedrals discovered ({}):\n".format(len(modes["flexible"])))
    output.write("\t" + "\t".join([ str(_)+"\n" for _ in modes["flexible"]]))
    
    print("\nInflexible Dihedrals discovered ({}):\n".format(len(modes["inflexible"])))
    print("\t" + "\t".join([ str(_)+"\n" for _ in modes["inflexible"]]))
    
    output.write("\nInflexible Dihedrals discovered ({}):\n".format(len(modes["inflexible"])))
    output.write("\t" + "\t".join([ str(_)+"\n" for _ in modes["inflexible"]]))
    
    output.close()
    return

if __name__ == "__main__": 
    main(sys.argv[1:])