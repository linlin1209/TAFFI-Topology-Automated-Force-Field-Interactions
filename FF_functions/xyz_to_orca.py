#!/bin/env python                                                                                                                                                             

from numpy import *
import sys,argparse,subprocess,random,os
from subprocess import PIPE
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

    parser.add_argument('-f', dest='functional', default='wB97X-D3',
                        help = 'Sets the functional for the calculation. (default: wB97X-D3; other typical options are B3LYP, and M062X)')

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
    args=parser.parse_args(argv)

    
    # Create fake option for function input
    option = dict(args.__dict__.items())

    # Check that the input is an .xyz file. 
    if args.coord_file.split('.')[-1] != 'xyz':
        print( "ERROR: Check to ensure that the input file is in .xyz format.")
        return
    
    # Create geo_dict for function input
    geo_dict = {}
    Elements,Geometry = xyz_parse(args.coord_file)
    inchikey = GetInchi(Elements,Geometry)
    geo_dict["elements"] = Elements
    geo_dict["geo"] = Geometry
    
    other_opt = {}
    fun(option,inchikey,geo_dict,other_opt)
    
    


def fun(config,inchikey,geo_dict,other_opt):

    option = miss2default(config,other_opt)
    proc_num = option['proc_num']
    random_factor = option['random_factor']
    constraints = option['constraints']
    charge = option['charge']
    functional = option['functional']
    basis = option['basis']
    multiplicity = option['multiplicity']
    eps = option['eps']
    quad_opt = option['quad_opt']
    no_freq = option['no_freq']
    no_opt = option['no_opt']
    D3_option = option['D3_option']
    RIJCOSX_option = option['RIJCOSX_option']
    free_elements = option['free_elements']

 
    # Parse inputs from config
    # Need to be modified to be able to do the original input
    inchikey = str(inchikey)
    output_name = inchikey

    # Parse inputs
    proc_num = int(proc_num)
    if proc_num > 8: print("ERROR: Too many procs requested, ORCA only supports up to 8 core parallelization. Exiting...");quit()
    random_factor = float(random_factor)
    if output_name != []:
        output_name = str(output_name)
    #else:
    #    output_name = '.'.join([ i for i in coord_file.split('.')[:-1]])+'.in'
    constraints = constraints.split()
    if False in [ i in ["bond","angle","dihedral"] for i in constraints ]: print("ERROR: only bond angle and dihedral are valid constraint options. Exiting..."); quit()
    if False not in [ i in constraints for i in ["bond","angle","dihedral"] ]: 
        print("ERROR: if all bonds, angles, and dihedrals are constrained then no geometry optimization can be performed. Exiting..."); quit()
    free_elements = free_elements.split()
    charge = int(charge)
    multiplicity = int(multiplicity)
    eps = float(eps)

    # Check that the input is an .xyz file. 
    #if coord_file.split('.')[-1] != 'xyz':
    #    print("ERROR: Check to ensure that the input file is in .xyz format.")
    #    return
    
    if basis == None:
        if charge > -0.9:
            basis = "def2-TZVP"
        else:
            basis = "ma-def2-TZVP"

    # Extract Element list and Coord list from the file
    #Elements,Geometry = xyz_parse(coord_file) ##should be modified, directly used geometry and element from frag gen inter
    Elements = geo_dict["elements"]
    Geometry = geo_dict["geo"]
    # the input inchikey is the inter inchikey
    inchi_intra = GetInchi(Elements,Geometry)

    # Generate adjacency table
    Adj_mat = Table_generator(Elements,Geometry)
    
    # Parse constrained modes
    Bonds_list = []
    Angles_list = []
    Dihedrals_list = []
    if "bond" in constraints:
        Bonds_list = Find_all_bonds(Adj_mat)

        # Remove constraints involving free atoms
        del_list = []
        for count_i,i in enumerate(Bonds_list):
            if Elements[i[0]] in free_elements or Elements[i[1]] in free_elements:
                del_list += [count_i]
                print("omitting constraint {} ({})".format(i,Elements[i[0]]+'-'+Elements[i[1]]))

        Bonds_list = [ i for count_i,i in enumerate(Bonds_list) if count_i not in del_list ]

    if "angle" in constraints:
        Angles_list = Find_all_angles(Adj_mat)

        # Remove constraints involving free atoms
        del_list = []
        for count_i,i in enumerate(Angles_list):
            if Elements[i[0]] in free_elements or Elements[i[2]] in free_elements: 
                del_list += [count_i]
                print("omitting constraint {} ({})".format(i,Elements[i[0]]+'-'+Elements[i[1]]+'-'+Elements[i[2]]))

        Angles_list = [ i for count_i,i in enumerate(Angles_list) if count_i not in del_list ]

    if "dihedral" in constraints:
        Dihedrals_list = Find_all_dihedrals(Adj_mat)

        # Remove constraints involving free atoms
        del_list = []
        for count_i,i in enumerate(Dihedrals_list):
            if Elements[i[0]] in free_elements or Elements[i[3]] in free_elements: 
                del_list += [count_i]
                print("omitting constraint {} ({})".format(i,Elements[i[0]]+'-'+Elements[i[1]]+'-'+Elements[i[2]]+'-'+Elements[i[3]]))

        Dihedrals_list = [ i for count_i,i in enumerate(Dihedrals_list) if count_i not in del_list ]

    # Use the total number of constraints as a check for the writing of the constraints block in the input file
    N_con = len(Bonds_list) + len(Angles_list) + len(Dihedrals_list)

    # Make inchi_folder
    if os.path.isdir("{}/Intra/geoopt_{}".format(inchi_intra,inchikey)) is False:
      os.makedirs("{}/Intra/geoopt_{}".format(inchi_intra,inchikey))
    
    # Write input file
    with open("{}/Intra/geoopt_{}/geoopt.in".format(inchi_intra,inchikey),'w') as f:
        f.write("# Geometry optimization for {}\n".format(inchi_intra))
        if proc_num == 1:
            if no_opt == 1:
                f.write("! {} {} TIGHTSCF Grid5 FinalGrid6 {}{}CHELPG\n".format(functional,basis,D3_option,RIJCOSX_option))
            else:
                #f.write("! {} {} TIGHTSCF TIGHTOPT Grid5 FinalGrid6 Opt {}{}CHELPG\n".format(functional,basis,D3_option,RIJCOSX_option))
                f.write("! {} {} TIGHTSCF Opt Grid5 FinalGrid6 {}{}CHELPG\n".format(functional,basis,D3_option,RIJCOSX_option))
        else:
            if no_opt == 1:
                f.write("! {} {} TIGHTSCF Grid5 FinalGrid6 {}{}CHELPG PAL{}\n".format(functional,basis,D3_option,RIJCOSX_option,proc_num))
            else:
                #f.write("! {} {} TIGHTSCF TIGHTOPT Grid5 FinalGrid6 Opt {}{}CHELPG PAL{}\n".format(functional,basis,D3_option,RIJCOSX_option,proc_num))
                f.write("! {} {} TIGHTSCF Opt Grid5 FinalGrid6 {}{}CHELPG PAL{}\n".format(functional,basis,D3_option,RIJCOSX_option,proc_num))

        # Perform frequency analysis
        if no_freq == 0:
            f.write("! Freq\n")

        # Add PCM model
        if eps > 0.0:
            f.write("\n%cosmo epsilon {}\nend\n".format(eps)) 

        # Add quadrupole calculation
        if quad_opt == 1:
            f.write("\n%elprop\nDipole True\nQuadrupole True\nend\n")

        if no_opt == 0:
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
        f.write("\n%base \"geo_opt\"\n\n* xyz {} {}\n".format(charge,multiplicity))
        for count_i,i in enumerate(Geometry):
            f.write("  {:3s}".format(Elements[count_i]))
            for j in i:
                f.write(" {:< 16.8f}".format((random.random()*2.0-1.0)*random_factor+j))
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

def GetInchi(Elements,Geometry): #convert geo2inchi key

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('out.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    # create inchikey
    p = subprocess.Popen("obabel -ixyz out.xyz -oinchikey".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
    inchikey, err = p.communicate()
    inchikey = str(inchikey, 'utf-8').split('\n')[0] ##might be a problem if there are characters utf-8 can't handle
    # remove xyz file
    os.remove("out.xyz")
   
    return inchikey

def miss2default(config,other_opt):

 
    # Create default dict
    options = ['proc_num','random_factor','output_name','constraints','charge','functional','basis','multiplicity','eps','quad_opt','no_freq','no_opt','D3_option','RIJCOSX_option','free_elements']

    defaults = [1, 0.0, [], '', 0, 'B3LYP', None, 1, 0.0, 0, 0, 0, "D3BJ ","",' ']  

    N_position = int (len(options) - len(defaults))

    default = {}
    for count_i,i in enumerate(defaults):
      default[options[N_position + count_i]] = i


    # Combine config and other_opt
    option = {}
    for key in config:
      option[key] = config[key]

    #other_opt's priority > config
    for key in other_opt:
      option[key] = other_opt[key]

    missing = [ i in option for i in options]
    
    
    # set missing option to default
    for count_i,i in enumerate(missing):
      if i is False:
         option[options[count_i]] = default[options[count_i]]

    return option
if __name__ == "__main__":
   main(sys.argv[1:])
