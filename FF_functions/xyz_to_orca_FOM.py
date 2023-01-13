#!/bin/env python                                                                                                                                                             

from numpy import *
import sys,argparse,random,os
from scipy.spatial.distance import *
from copy import deepcopy

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from adjacency import *
from file_parsers import *
from id_types import *

def main(argv):

    parser = argparse.ArgumentParser(description='Converts an xyz into an orca geometry optimization input file. Allows for the randomization of atom positions')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_file', help = 'The input file. (currently must be an xyz with the atom types in the fourth column and the ion listed last.')

    #optional arguments
    parser.add_argument('-p', dest='proc_num', default=1,
                        help = 'Sets the number of processors to run the orca job on (default: 1)')

    parser.add_argument('-r', dest='random_factor', default=0.0,
                        help = 'Scaling of the position random number. Everytime an atom is placed, a random number, R, between 1 and 0 is drawn, and the positions moved by '+\
                               'R*scaling factor. By default scaling_fact=0.0 and the positions are copied as in the .xyz file. (default: 0.0)')

    parser.add_argument('-o', dest='output_name', default=[],
                        help = 'Determines the name of the output file. If none is supplied then the default behavior is to use the base of the input filename (default: base_filename)')

    parser.add_argument('-c', dest='constraints', default='',
                        help = 'Determines the constraints to apply to the geometry optimization. This argument expects a space-delimited string. Valid arguments are bond angle and dihedral. '+\
                               'By default, no modes are contrained during the geometry optimization.')

    parser.add_argument('-q', dest='charge', default=0,
                        help = 'Sets the total charge for the calculation. (default: 0)')

    parser.add_argument('-f', dest='functional', default='B3LYP',
                        help = 'Sets the functional for the calculation. (default: B3LYP; other typical options are WB97X-D3, and M062X)')

    parser.add_argument('-b', dest='basis', default=None,
                        help = 'Sets the basis set for the calculation. (default: def2-TZVP for neutral and cationic species, and ma-def2-TZVP for anions)')

    parser.add_argument('-m', dest='multiplicity', default=1,
                        help = 'Sets the multiplicity for the calculation. (default: 1)')

    parser.add_argument('-eps', dest='eps', default=0.0,
                        help = 'Sets the dielectric for the cosmo pcm model. By default no pcm is used.')

    parser.add_argument('--quad', dest='quad_opt', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the orca job will calculate the quadrupole. (default: off)')

    parser.add_argument('--no_freq', dest='no_freq', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the orca job will avoid performing a frequency calculation. Default behavior is to perform the frequency calculation.')

    parser.add_argument('--no_opt', dest='no_opt', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the orca job will avoid performing a geometry optimization. Default behavior is to perform the geometry optimization.')

    parser.add_argument('--no_D3', dest='D3_option', default="D3BJ ", action='store_const', const="",
                        help = 'When this flag is present, the D3-dispersion correction is disabled. By default the dispersion correction is used.')

    parser.add_argument('--RIJCOSX', dest='RIJCOSX_option', default="", action='store_const', const="RIJCOSX ",
                        help = 'When this flag is present, the RIJCOX approximation is enabled. This is usually an excellent approximation and greatly speeds up hybrid calculations. (Default: off)')

    parser.add_argument('--free', dest='free_elements', default=' ', 
                        help = 'This argument expects a list of element labels. When contraints are supplied, any modes involving these elements in a terminal position are left unconstrained. '+\
                               'For example, if dihedrals are being constrained (via -c "dihedral") and --free "H" is supplied, then dihedrals like H-C-C-C, O-C-C-H, and H-C-C-H will be left '+\
                               'unconstrained. (default: none)') 

    # Make relevant inputs lowercase
    args=parser.parse_args()

    # Parse inputs
    args.proc_num = int(args.proc_num)
    args.random_factor = float(args.random_factor)
    if args.output_name != []:
        args.output_name = str(args.output_name)
    else:
        args.output_name = '.'.join([ i for i in args.coord_file.split('.')[:-1]])+'.in'
    args.constraints = args.constraints.split()
    if False in [ i in ["bond","angle","dihedral"] for i in args.constraints ]: print("ERROR: only bond angle and dihedral are valid constraint options. Exiting..."); quit()
    if False not in [ i in args.constraints for i in ["bond","angle","dihedral"] ]: 
        print("ERROR: if all bonds, angles, and dihedrals are constrained then no geometry optimization can be performed. Exiting..."); quit()
    args.free_elements = args.free_elements.split()
    args.charge = int(args.charge)
    args.multiplicity = int(args.multiplicity)
    args.eps = float(args.eps)

    # Check that the input is an .xyz file. 
    if args.coord_file.split('.')[-1] != 'xyz':
        print("ERROR: Check to ensure that the input file is in .xyz format.")
        return
    
    if args.basis == None:
        if args.charge > -0.9:
            args.basis = "def2-TZVP"
        else:
            args.basis = "ma-def2-TZVP"

    # Extract Element list and Coord list from the file
    Elements,Geometry = xyz_parse(args.coord_file)

    # Generate adjacency table
    Adj_mat = Table_generator(Elements,Geometry)
    
    # Parse constrained modes
    Bonds_list = []
    Angles_list = []
    Dihedrals_list = []
    if "bond" in args.constraints:
        Bonds_list = Find_all_bonds(Adj_mat)

        # Remove constraints involving free atoms
        del_list = []
        for count_i,i in enumerate(Bonds_list):
            if Elements[i[0]] in args.free_elements or Elements[i[1]] in args.free_elements:
                del_list += [count_i]
                print("omitting constraint {} ({})".format(i,Elements[i[0]]+'-'+Elements[i[1]]))

        Bonds_list = [ i for count_i,i in enumerate(Bonds_list) if count_i not in del_list ]

    if "angle" in args.constraints:
        Angles_list = Find_all_angles(Adj_mat)

        # Remove constraints involving free atoms
        del_list = []
        for count_i,i in enumerate(Angles_list):
            if Elements[i[0]] in args.free_elements or Elements[i[2]] in args.free_elements: 
                del_list += [count_i]
                print("omitting constraint {} ({})".format(i,Elements[i[0]]+'-'+Elements[i[1]]+'-'+Elements[i[2]]))

        Angles_list = [ i for count_i,i in enumerate(Angles_list) if count_i not in del_list ]

    if "dihedral" in args.constraints:
        Dihedrals_list = Find_all_dihedrals(Adj_mat)

        # Remove constraints involving free atoms
        del_list = []
        for count_i,i in enumerate(Dihedrals_list):
            if Elements[i[0]] in args.free_elements or Elements[i[3]] in args.free_elements: 
                del_list += [count_i]
                print("omitting constraint {} ({})".format(i,Elements[i[0]]+'-'+Elements[i[1]]+'-'+Elements[i[2]]+'-'+Elements[i[3]]))

        Dihedrals_list = [ i for count_i,i in enumerate(Dihedrals_list) if count_i not in del_list ]

    # Use the total number of constraints as a check for the writing of the constraints block in the input file
    N_con = len(Bonds_list) + len(Angles_list) + len(Dihedrals_list)

    # Write input file
    with open(args.output_name,'w') as f:
        f.write("# Geometry optimization for {}\n".format(args.coord_file))
        if args.no_opt == 1:
             f.write("! {} {} TIGHTSCF Grid5 FinalGrid6 {}{}CHELPG\n".format(args.functional,args.basis,args.D3_option,args.RIJCOSX_option))
        else:
             f.write("! {} {} TIGHTSCF TIGHTOPT Grid5 FinalGrid6 Opt {}{}CHELPG\n".format(args.functional,args.basis,args.D3_option,args.RIJCOSX_option))

        # Perform frequency analysis
        if args.no_freq == 0:
            f.write("! Freq\n")

        # Add PCM model
        if args.eps > 0.0:
            f.write("\n%cosmo epsilon {}\nend\n".format(args.eps)) 

        # procs option
        if args.proc_num != 1:
            f.write("\n%pal nprocs {}\nend\n".format(args.proc_num))

        # Add quadrupole calculation
        if args.quad_opt == 1:
            f.write("\n%elprop\nDipole True\nQuadrupole True\nend\n")

        if args.no_opt == 0:
            f.write("\n%geom MaxIter 500\n")

            # Write Geometry constraint(s)    
            if N_con > 0:
                f.write('  Constraints\n')
                for i in Bonds_list:
                    f.write('    {{ B {} {} C }}\n'.format(i[0],i[1]))
                for i in Angles_list:
                    f.write('    {{ A {} {} {} C }}\n'.format(i[0],i[1],i[2]))
                for i in Dihedrals_list:
                    f.write('    {{ D {} {} {} {} C }}\n'.format(i[0],i[1],i[2],i[3]))
                f.write('  end\n')
            f.write('end\n\n')

        # Write job name and geometry
        f.write("\n%base \"geo_opt\"\n\n* xyz {} {}\n".format(args.charge,args.multiplicity))
        for count_i,i in enumerate(Geometry):
            f.write("  {:3s}".format(Elements[count_i]))
            for j in i:
                f.write(" {:< 16.8f}".format((random.random()*2.0-1.0)*args.random_factor+j))
            f.write("\n")
        f.write("*\n")

# This function identifies all bonds based on 1-bond deep walks within the adjacency matrix
def Find_all_bonds(adj_mat):

    bonds = []
    
    # Each row in the adjacency matrix is used to seed a 2-bond deep search
    for count_i,i in enumerate(adj_mat):

        # tmp holds the bonds to the i-th atom
        tmp = [ [count_i,count_j] for count_j,j in enumerate(i) if j == 1 ]

        # Only add new bonds
        bonds += [ j for j in tmp if ( j not in bonds and [j[1],j[0]] not in bonds ) ]

    bond_atoms = [ tuple(i) for i in bonds ]

    return bond_atoms

# This function identifies all angles based on 2-bond deep walks within the adjacency matrix
def Find_all_angles(adj_mat):

    angles = []
    
    # Each row in the adjacency matrix is used to seed a 2-bond deep search
    for count_i,i in enumerate(adj_mat):

        # current holds the angle seedlings
        current = [ [count_i,count_j] for count_j,j in enumerate(i) if j == 1 ]

        # new holds the updated list of angles at the end of each bond search
        new = []
        iteration = 0 

        # recursively add bonds to new until no new bonds are found
        while new != current:
            if iteration == 0: new = deepcopy(current)
            current = deepcopy(new)

            # Iterate over the seedlings and add connections/spawn new seedlings
            for count_j,j in enumerate(current):
                if len(j) == 3: continue                                                                            # avoid adding new elements to full angles
                connections = [ count_k for count_k,k in enumerate(adj_mat[j[-1]]) if k == 1 and count_k not in j ] # add new connections not already in the angle
                for count_k,k in enumerate(connections):
                    if count_k == 0: new[count_j] += [k]                                                            # add the first connection to the existing seedling
                    else: new += [ j + [k] ]                                                                        # add the rest as new elements
            iteration += 1             

        # Only add full angles and new angles to the angles list
        tmp = [ j for j in new if len(j) == 3 ]
        angles += [ j for j in tmp if ( j not in angles and [j[2],j[1],j[0]] not in angles ) ]

    angle_atoms = [ tuple(i) for i in angles ]

    return angle_atoms

# This function identifies all dihedrals based on 3-bond deep walks within the adjacency matrix
def Find_all_dihedrals(adj_mat):

    dihedrals = []
    
    # Each row in the adjacency matrix is used to seed a 3-bond deep search
    for count_i,i in enumerate(adj_mat):

        # current holds the dihedral seedlings
        current = [ [count_i,count_j] for count_j,j in enumerate(i) if j == 1 ]

        # new holds the updated list of dihedrals at the end of each bond search
        new = []
        iteration = 0 

        # recursively add bonds to new until no new bonds are found
        while new != current:
            if iteration == 0: new = deepcopy(current)
            current = deepcopy(new)

            # Iterate over the seedlings and add connections/spawn new seedlings
            for count_j,j in enumerate(current):
                if len(j) == 4: continue                                                                            # avoid adding new elements to full dihedrals
                connections = [ count_k for count_k,k in enumerate(adj_mat[j[-1]]) if k == 1 and count_k not in j ] # add new connections not already in the dihedral
                for count_k,k in enumerate(connections):
                    if count_k == 0: new[count_j] += [k]                                                            # add the first connection to the existing seedling
                    else: new += [ j + [k] ]                                                                        # add the rest as new elements
            iteration += 1             

        # Only add full dihedrals and new dihedrals to the dihedrals list
        tmp = [ j for j in new if len(j) == 4 ]
        dihedrals += [ j for j in tmp if ( j not in dihedrals and [j[3],j[2],j[1],j[0]] not in dihedrals ) ]

    dihedral_atoms = [ tuple(i) for i in dihedrals ]

    return dihedral_atoms

if __name__ == "__main__":
   main(sys.argv[1:])
