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
from builtins import any  # any has a namespace conflict with numpy. XXX WE REALLY NEED TO CLEAN UP THESE NAMESPACE ISSUES
import json, codecs

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a .xyz file and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    parser.add_argument('-FF', dest='FF_db', default='',
                        help = 'This variable holds the filename(s) of the .db files holding the necessary force-field parameters (default: FF.db)')

    parser.add_argument('-o', dest='outputname', default='frags',
                        help = 'Controls the output folder name for the fragments. (default: frags)')

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'Controls the bond search depth for identifying unique atom types (default: 2)')

    parser.add_argument('-R', dest='recursive_opt', default=False, action='store_const',const=True,
                        help = 'When this flag is present, filenames and wildcards will be searched for recursively.')

    parser.add_argument('-q', dest='q_tot', default=0, 
                        help = 'The total charge on the generated fragments. (default: 0)')

    parser.add_argument('--keep_terminal', dest='keep_terminal_opt', default=False, action='store_const',const=True,
                        help = 'When this flag is present, terminal atoms at the edges of fragments are retained. Default is to perform a full hydrogenation consistent with the hybridization of the atomtype. (default: off)')

    parser.add_argument('--avoid_frags', dest='avoid_frags_opt', default=False, action='store_const',const=True,
                        help = 'When this flag is present, the xyzs are used directly for the parameterization (default: False)')

    # Make relevant inputs lowercase
    args=parser.parse_args()    

    # Convert coords_files and N into lists
    args.FF_db = args.FF_db.split()

    # Parse additional inputs 
    args.gens = int(args.gens)
    args.q_tot = int(args.q_tot)

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
        print("ERROR: Check to ensure that the input file(s) are in .xyz format.")
        return
    if False in [ os.path.isfile(i) for i in args.coord_files ]:
        print("ERROR: Could not find file {}. Exiting...".format(next([ i for i in args.coord_files if os.path.isfile(i) == False ])))
    for i in args.FF_db:
        if os.path.isfile(i) == False:
            print("ERROR: FF file {} couldn't not be found. Exiting...".format(i))
            return

    # Generate base filename for naming output
    #Filename = args.outputname

    # Make directory to hold the output    
    #if os.path.exists(Filename):
    #    print('The desired output folder ({}) already exists. Exiting to avoid overwriting data...'.format(Filename))
    #    return
    #else:
    #    os.makedirs(Filename)
    #    Filename = Filename+'/'+Filename
    sys.stdout = Logger(args.outputname)
    print("PROGRAM CALL: python frag_gen.py {}\n".format(' '.join([ i for i in argv])))

    # Generate a dictionary of the various modes contained in the supplied database files
    FF_db,modes_from_FF = parse_FF_params(args.FF_db)

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
    print("\n"+"*"*120)
    print("* {:^116s} *".format("Parsing modes and generating all-trans optimized conformers for the following files"))
    print("*"*120 + "\n")    
    for file in args.coord_files:        

        # Print the name
        print("  {}".format(file))

        # Initialize subdictionary
        geo_dict[file] = {}

        # Extract Element list and Coord list from the file
        geo_dict[file]["elements"],geo_dict[file]["geo"] = xyz_parse(file)

        # Generate adjacency table
        geo_dict[file]["adj_mat"] = Table_generator(geo_dict[file]["elements"],geo_dict[file]["geo"])

        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        hybridizations = Hybridization_finder(geo_dict[file]["elements"],geo_dict[file]["adj_mat"])

        # Find atom types
        geo_dict[file]["atom_types"] = id_types(geo_dict[file]["elements"],geo_dict[file]["adj_mat"],args.gens,hybridizations,geo_dict[file]["geo"])

        # Canonicalize the geometry
        geo_dict[file]["elements"],geo_dict[file]["geo"],geo_dict[file]["adj_mat"],geo_dict[file]["atom_types"] = canon_geo(geo_dict[file]["elements"],geo_dict[file]["geo"],\
                                                                                                                            geo_dict[file]["adj_mat"],geo_dict[file]["atom_types"])

        # Determine bonding matrix for the compound
        geo_dict[file]["bond_mat"] = find_lewis(geo_dict[file]["atom_types"], geo_dict[file]["adj_mat"], q_tot=args.q_tot, b_mat_only=True,verbose=False)

        # Find modes (returns all modes in canonicalized format)
        bonds_tmp,angles_tmp,dihedrals_tmp,one_fives_tmp,bond_types_tmp,angle_types_tmp,dihedral_types_tmp,one_five_types_tmp = Find_modes(geo_dict[file]["adj_mat"],geo_dict[file]["atom_types"],geo_dict[file]["bond_mat"],return_all=1)

        # Canonicalize modes
        for count_i,i in enumerate(bonds_tmp):
            bond_types_tmp[count_i],bonds_tmp[count_i] = canon_bond(bond_types_tmp[count_i],bonds_tmp[count_i])
        for count_i,i in enumerate(angles_tmp):
            angle_types_tmp[count_i],angles_tmp[count_i] = canon_angle(angle_types_tmp[count_i],angles_tmp[count_i])
        for count_i,i in enumerate(dihedrals_tmp):
            dihedral_types_tmp[count_i],dihedrals_tmp[count_i] = canon_dihedral(dihedral_types_tmp[count_i],dihedrals_tmp[count_i])

        # Prune the bonds angles and dihedrals
        del_ind = set([ count_i for count_i,i in enumerate(bond_types_tmp) if i in FF_db["bonds"] or i in bond_types ])
        bonds_tmp = [ i for count_i,i in enumerate(bonds_tmp) if count_i not in del_ind ]
        bond_types_tmp = [ i for count_i,i in enumerate(bond_types_tmp) if count_i not in del_ind ]
        del_ind = set([ count_i for count_i,i in enumerate(angle_types_tmp) if i in FF_db["angles"] or i in angle_types ])
        angles_tmp = [ i for count_i,i in enumerate(angles_tmp) if count_i not in del_ind ]
        angle_types_tmp = [ i for count_i,i in enumerate(angle_types_tmp) if count_i not in del_ind ]
        del_ind = set([ count_i for count_i,i in enumerate(dihedral_types_tmp) if i in FF_db["dihedrals"] or i in dihedral_types ])
        dihedrals_tmp = [ i for count_i,i in enumerate(dihedrals_tmp) if count_i not in del_ind ]
        dihedral_types_tmp = [ i for count_i,i in enumerate(dihedral_types_tmp) if count_i not in del_ind ]
        del_ind = set([ count_i for count_i,i in enumerate(dihedral_types_tmp) if i in FF_db["dihedrals_harmonic"] or i in dihedral_types ])
        dihedrals_tmp = [ i for count_i,i in enumerate(dihedrals_tmp) if count_i not in del_ind ]
        dihedral_types_tmp = [ i for count_i,i in enumerate(dihedral_types_tmp) if count_i not in del_ind ]

        # Append tmp's to lists
        bonds += bonds_tmp
        bond_types += bond_types_tmp
        bond_files += [ file for j in bonds_tmp ]
        angles += angles_tmp
        angle_types += angle_types_tmp
        angle_files += [ file for j in angles_tmp ]
        dihedrals += dihedrals_tmp
        dihedral_types += dihedral_types_tmp
        dihedral_files += [ file for j in dihedrals_tmp ]

        # If there are ring-dihedrals then the geometry will eventually get unwrapped. 
        # Here the bond, angle, and dihedrals from the unwrapped fragment are added as a new geometry.
        # The following loop is a near identical copy of the above, except that there is a call to Find_modes
        # that handles the generation of the dihedral geometry
        dihedral_types_tmp_0 = deepcopy(dihedral_types_tmp)
        dihedrals_tmp_0 = deepcopy(dihedrals_tmp)
        if args.avoid_frags_opt is False:
            for count_d,d in enumerate(dihedral_types_tmp_0):
                if "R" in d[1] and "R" in d[2]:
                    new = "{}-{}".format(file,count_d)
                    geo_dict[new] = {}

                    # Get the linear analog of the dihedaral (d) containing molecule
                    _,geo_dict[new]["geo"],geo_dict[new]["adj_mat"],geo_dict[new]["elements"],_,_ = \
                        mode_frag(dihedrals_tmp_0[count_d],geo_dict[file]["geo"],geo_dict[file]["adj_mat"],geo_dict[file]["bond_mat"],geo_dict[file]["elements"],\
                          geo_dict[file]["atom_types"],args.gens,args.q_tot,keep_terminal=args.keep_terminal_opt)

                    # Find atom types
                    geo_dict[new]["atom_types"] = id_types(geo_dict[new]["elements"],geo_dict[new]["adj_mat"],args.gens)

                    # Canonicalize the geometry
                    geo_dict[new]["elements"],geo_dict[new]["geo"],geo_dict[new]["adj_mat"],geo_dict[new]["atom_types"] = canon_geo(geo_dict[new]["elements"],geo_dict[new]["geo"],\
                                                                                                                                        geo_dict[new]["adj_mat"],geo_dict[new]["atom_types"])

                    # Determine bonding matrix for the compound
                    geo_dict[new]["bond_mat"] = find_lewis(geo_dict[new]["atom_types"], geo_dict[new]["adj_mat"], q_tot=args.q_tot, b_mat_only=True,verbose=False)                                   

                    # Find modes (returns all modes in canonicalized format)
                    bonds_tmp,angles_tmp,dihedrals_tmp,one_fives_tmp,bond_types_tmp,angle_types_tmp,dihedral_types_tmp,one_five_types_tmp = Find_modes(geo_dict[new]["adj_mat"],geo_dict[new]["atom_types"],\
                                                                                                                                                       geo_dict[new]["bond_mat"],return_all=1)

                    # Canonicalize modes
                    for count_i,i in enumerate(bonds_tmp):
                        bond_types_tmp[count_i],bonds_tmp[count_i] = canon_bond(bond_types_tmp[count_i],bonds_tmp[count_i])
                    for count_i,i in enumerate(angles_tmp):
                        angle_types_tmp[count_i],angles_tmp[count_i] = canon_angle(angle_types_tmp[count_i],angles_tmp[count_i])
                    for count_i,i in enumerate(dihedrals_tmp):
                        dihedral_types_tmp[count_i],dihedrals_tmp[count_i] = canon_dihedral(dihedral_types_tmp[count_i],dihedrals_tmp[count_i])

                    # Prune the bonds angles and dihedrals
                    del_ind = set([ count_i for count_i,i in enumerate(bond_types_tmp) if i in FF_db["bonds"] or i in bond_types ])
                    bonds_tmp = [ i for count_i,i in enumerate(bonds_tmp) if count_i not in del_ind ]
                    bond_types_tmp = [ i for count_i,i in enumerate(bond_types_tmp) if count_i not in del_ind ]
                    del_ind = set([ count_i for count_i,i in enumerate(angle_types_tmp) if i in FF_db["angles"] or i in angle_types ])
                    angles_tmp = [ i for count_i,i in enumerate(angles_tmp) if count_i not in del_ind ]
                    angle_types_tmp = [ i for count_i,i in enumerate(angle_types_tmp) if count_i not in del_ind ]
                    del_ind = set([ count_i for count_i,i in enumerate(dihedral_types_tmp) if i in FF_db["dihedrals"] or i in dihedral_types ])
                    dihedrals_tmp = [ i for count_i,i in enumerate(dihedrals_tmp) if count_i not in del_ind ]
                    dihedral_types_tmp = [ i for count_i,i in enumerate(dihedral_types_tmp) if count_i not in del_ind ]
                    del_ind = set([ count_i for count_i,i in enumerate(dihedral_types_tmp) if i in FF_db["dihedrals_harmonic"] or i in dihedral_types ])
                    dihedrals_tmp = [ i for count_i,i in enumerate(dihedrals_tmp) if count_i not in del_ind ]
                    dihedral_types_tmp = [ i for count_i,i in enumerate(dihedral_types_tmp) if count_i not in del_ind ]

                    # Append tmp's to lists
                    bonds += bonds_tmp
                    bond_types += bond_types_tmp
                    bond_files += [ new for j in bonds_tmp ]
                    angles += angles_tmp
                    angle_types += angle_types_tmp
                    angle_files += [ new for j in angles_tmp ]
                    dihedrals += dihedrals_tmp
                    dihedral_types += dihedral_types_tmp
                    dihedral_files += [ new for j in dihedrals_tmp ]

    # Initialize mass_dict (used for identifying the dihedral among a coincident set that will be explicitly scanned)
    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Collect the first instance of each bond, angle, and dihedral in need of parametrization
    print("\n"+"*"*120)
    print("* {:^116s} *".format("Determining which modes are in need of parametrization"))
    print("*"+"-"*118+"*")
    print("* {:<116s} *".format("NOTE: Dihedrals labled with a (*) are in need of parametrization, but are coincident with other dihedrals"))
    print("* {:<116s} *".format("      and do not need to be independently scanned"))
    print("*"*120 + "\n")    
    print("Modes in need of parametrization:")
    mode_instances = []
    instance_types = []    

    # Bond collection
    print("\n\tBond Types ({}):\n".format(len(set(bond_types))))
    for i in set(bond_types):
        print("\t\t{}".format(i))
        for count_j,j in enumerate(bond_types):
            if i == j:
                mode_instances += [(bond_files[count_j],bonds[count_j])]
                instance_types += [i]
                break

    # Angle collection
    print("\n\tAngle Types ({}):\n".format(len(set(angle_types))))
    for i in set(angle_types):
        print("\t\t{}".format(i))
        for count_j,j in enumerate(angle_types):
            if i == j:
                mode_instances += [(angle_files[count_j],angles[count_j])]
                instance_types += [i]
                break

    # Dihedral collection (longer because it is necessary to identify and remove redundant coincident dihedral scans.
    print("\n\tDihedral Types ({}) (*=implicit):\n".format(len(set(dihedral_types))))
    implicit_dihedrals = []  
    assigned_dihedrals = []
    for i in set(dihedral_types):

        # Skip the dihedral if it has already been assigned
        if i in assigned_dihedrals: continue

        # Find the index for the first instance of dihedral "i" in the dihedral_types list
        ind = next( count_j for count_j,j in enumerate(dihedral_types) if i == j )
        file = dihedral_files[ind]

        # Harmonic dihedrals (TAFFI treats 2-3 atoms connected by a double bond explicitly as harmonic functions. 
        if i[4] == "harmonic":
            assigned_dihedrals += [i]
            mode_instances     += [(dihedral_files[ind],dihedrals[ind])]
            instance_types     += [i]
            print("\t\t{}".format(i))
            continue

        # Gather coincident dihedrals (dihedrals that share the 2-3 bond with this dihedral)
        left_atoms  = [ count_k for count_k,k in enumerate(geo_dict[file]["adj_mat"][dihedrals[ind][1]]) if k == 1 and count_k != dihedrals[ind][2] ]
        right_atoms = [ count_k for count_k,k in enumerate(geo_dict[file]["adj_mat"][dihedrals[ind][2]]) if k == 1 and count_k != dihedrals[ind][1] ]
        coincident_dihedrals = [ tuple([k])+dihedrals[ind][1:3]+tuple([m]) for k in left_atoms for m in right_atoms ]
        coincident_types     = [ tuple([ geo_dict[file]["atom_types"][m] for m in k ] + ["opls"]) for k in coincident_dihedrals ]

        # Canonicalize the dihedrals (probably not necessary but done throughout all scripts just in case somebody does something that relies on it)
        for count_k,k in enumerate(coincident_dihedrals):
            coincident_types[count_k],coincident_dihedrals[count_k] = canon_dihedral(coincident_types[count_k],coincident_dihedrals[count_k])

        # Find the dihedral that will be explicitly scanned (the dihedral with the largest combined mass on the 1 and 4 atoms gets explicitly scanned)
        max_mass = 0.0
        for count_k,k in enumerate(coincident_dihedrals):
            current_mass = mass_dict[geo_dict[file]["elements"][k[0]]] + mass_dict[geo_dict[file]["elements"][k[3]]]
            if current_mass > max_mass:
                max_mass = current_mass
                keep = count_k

        # Add the implicitly treated dihedrals to the implicit_dihedrals list and update the mode_instances and instance_types lists with the explicit dihedral
        implicit_dihedrals += [ k for k in coincident_types if ( k != coincident_types[keep] and k != coincident_types[keep][::-1] ) ]                                
        mode_instances += [(dihedral_files[ind],coincident_dihedrals[keep])]   
        instance_types += [coincident_types[keep]]
        assigned_dihedrals += coincident_types
        print("\t\t{}".format(coincident_types[keep]))

    # Print the dihedrals that are being implicitly scanned
    for i in set(implicit_dihedrals):
        print("\t\t{} *".format(i))
    print("\n\tTotal modes being explicitly parametrized: {}".format(len(mode_instances)))

    # Generate fragments for each mode being parametrized
    print("\n"+"*"*120)
    print("* {:^116s} *".format("Beginning fragment generation"))
    print("*"+"-"*118+"*")
    print("* {:<116s} *".format("NOTE: Each mode will be parametrized using the smallest molecular fragment that is consistent with its"))
    print("* {:<116s} *".format("      topological definition"))
    print("*"*120 + "\n")
    modes_tot      = []
    modetypes_tot  = []
    atomtypes_tot  = []
    geos_tot       = []
    elements_tot   = []
    adj_mat_tot    = []
    bond_mat_tot   = []
    hash_lists_tot = []    
    for count_i,i in enumerate(mode_instances):

        # Use the original xyz
        if args.avoid_frags_opt:
            mode_ind        = i[1]
            frag_geo        = geo_dict[i[0]]["geo"]
            frag_adj_mat    = geo_dict[i[0]]["adj_mat"]
            frag_bond_mat   = geo_dict[i[0]]["bond_mat"]
            frag_elements   = geo_dict[i[0]]["elements"]
            frag_atom_types = geo_dict[i[0]]["atom_types"]
            frag_masses     = [ mass_dict[frag_elements[j]] for j in range(len(frag_atom_types)) ]
            frag_hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_k,frag_adj_mat,frag_masses),k) for count_k,k in enumerate(range(len(frag_atom_types))) ],reverse=True)) ]            

            # Update lists/arrays based on the sorted atoms
            mode_ind        = tuple([ atoms.index(j) for j in mode_ind ])
            frag_geo        = frag_geo[atoms]
            frag_adj_mat    = frag_adj_mat[atoms]
            frag_adj_mat    = frag_adj_mat[:,atoms]
            frag_bond_mat   = [ j[atoms] for j in frag_bond_mat ]
            frag_bond_mat   = [ j[:,atoms] for j in frag_bond_mat ]
            frag_elements   = [ frag_elements[j] for j in atoms ]
            frag_atom_types = [ frag_atom_types[j] for j in atoms ]

        # Generate canonical fragment for mode i
        else:
            mode_ind,frag_geo,frag_adj_mat,frag_elements,frag_atom_types,frag_hash_list = \
            mode_frag(i[1],geo_dict[i[0]]["geo"],geo_dict[i[0]]["adj_mat"],geo_dict[i[0]]["bond_mat"],geo_dict[i[0]]["elements"],\
                      geo_dict[i[0]]["atom_types"],args.gens,args.q_tot,keep_terminal=args.keep_terminal_opt)

            # Determine bonding matrix for the compound
            frag_bond_mat = find_lewis(frag_atom_types, frag_adj_mat, q_tot=args.q_tot, b_mat_only=True,verbose=False)

        # Build total lists
        modes_tot      += [mode_ind]
        modetypes_tot  += [instance_types[count_i]]
        atomtypes_tot  += [frag_atom_types]
        geos_tot       += [frag_geo]
        elements_tot   += [frag_elements]
        adj_mat_tot    += [frag_adj_mat]
        bond_mat_tot   += [frag_bond_mat]
        hash_lists_tot += [frag_hash_list]

    # Save one instance of each unique fragment and a frag-*.modes file to keep track of what modes to parametrize for each fragment
    geo_dict_out = {}  
    
    run_list = []
    N_frag = 0
    for count_i,i in enumerate(hash_lists_tot):

        # Skip already assigned fragments
        if count_i in run_list: continue        

        # Find equivalent fragments based on a comparison of the hash lists
        match_idx = []
        for count_j,j in enumerate(hash_lists_tot):
            if i == j:

                # Add the current match, and remap the atomtypes and mode indices
                match_idx += [count_j]
                atomtypes_tot[count_j],mapping,status = map_frag(atomtypes_tot[count_i],adj_mat_tot[count_i],bond_mat_tot[count_i],i,atomtypes_tot[count_j],adj_mat_tot[count_j],bond_mat_tot[count_j],j)
                modes_tot[count_j] = tuple([ mapping[k] for k in modes_tot[count_j] ])

                # Print problem fragment(s) in the case of failure status being returned by map_frag
                if status == 1:
                    write_xyz("incomp_fixed",elements_tot[count_i],geos_tot[count_i],atomtypes_tot[count_i])                    
                    write_xyz("incomp_mapped",elements_tot[count_j],geos_tot[count_j],atomtypes_tot[count_j])                    
                    quit()

        # Save the final transified geometry
        opt_geo = transify(geos_tot[match_idx[0]],adj_mat_tot[match_idx[0]],elements=elements_tot[match_idx[0]],opt_terminals=True,opt_final=True)
        #write_xyz("{}/frag-{}".format(args.outputname,N_frag),elements_tot[match_idx[0]],opt_geo,atomtypes_tot[match_idx[0]])
        
        #inchikey = GetInchi(elements_tot[match_idx[0]],opt_geo)
        #geo_dict_out[inchikey] = {}
        #geo_dict_out[inchikey]["elements"] = elements_tot[match_idx[0]]
        #geo_dict_out[inchikey]["geo"] = opt_geo
        #geo_dict_out[inchikey]["atom_types"] = atomtypes_tot[match_idx[0]]
        print("\n\tModes asssociated with fragment {:<12s}".format(str(N_frag)+":"))

        # Assemble the dictionary of modes being scanned and print diagnostics for the user
        # NOTE: dihedrals are treated differently from bonds and angles because in some cases two dihedrals have identical fragments and simply
        #       differ in the atomtypes/VDW/charges used when fitting the potential. 
        scanned_modes = {}
        placed = []
        for j in [ k for k in match_idx if len(modes_tot[k]) == 2]:
            scanned_modes[modetypes_tot[j]] = { "modes":[modes_tot[j]], "atomtypes":[atomtypes_tot[j]] }
            print("\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
        for j in [ k for k in match_idx if len(modes_tot[k]) == 3]:
            scanned_modes[modetypes_tot[j]] = { "modes":[modes_tot[j]], "atomtypes":[atomtypes_tot[j]] }
            print("\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
        for j in [ k for k in match_idx if len(modes_tot[k]) == 4]:
            print("\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
            if j in placed: continue
            scanned_modes[modetypes_tot[j]] = { "modes":[], "atomtypes":[] }
            
            # If j is a harmonic dihedral then redundancy is assessed on the basis of all four mode 
            # indices matching (this is necessary beacuse of the way that harmonic dihedrals are scanned)
            if 2 in [ n[modes_tot[j][1],modes_tot[j][2]] for n in bond_mat_tot[j] ]:

                # Iterate over all dihedrals that have identical atomic indices in their definition
                for m in [ n for n in match_idx if len(modes_tot[n]) == 4]:

                    # Save dihedrals that have indentical atomic indices
                    if modes_tot[j] == modes_tot[m] or modes_tot[j] == modes_tot[m][::-1]:
                        placed += [m]
                        scanned_modes[modetypes_tot[j]]["modes"]     += [modes_tot[m]]
                        scanned_modes[modetypes_tot[j]]["atomtypes"] += [atomtypes_tot[m]]

            # If j is a fourier dihedral then redundancy is assessed on the basic of the 2-3 atoms in the dihedral
            else:

                # Iterate over all dihedrals that have identical fragments
                for m in [ n for n in match_idx if len(modes_tot[n]) == 4]:

                    # Save dihedrals that have identical 2-3 atoms in their definition
                    if modes_tot[j][1:3] == modes_tot[m][1:3] or modes_tot[j][1:3] == modes_tot[m][1:3][::-1]:
                        placed += [m]
                        scanned_modes[modetypes_tot[j]]["modes"]     += [modes_tot[m]]
                        scanned_modes[modetypes_tot[j]]["atomtypes"] += [atomtypes_tot[m]]
        
        # Write the frag-*.modes file associated with this fragment

         # Make inchi_folder
        #for file in args.coord_files:        
        #    inchikey_in = GetInchi(geo_dict[file]["elements"],geo_dict[file]["geo"])
        #if os.path.isdir("{}".format(inchikey)) is False:
        #    os.makedirs("{}".format(inchikey))
        ## Write the frag-*.modes file associated with this fragment
        #write_modelist("{}/out_{}.modes".format(inchikey,inchikey_in),scanned_modes,bond_mat_tot[count_i])
        #write_xyz("{}/out".format(inchikey),elements_tot[match_idx[0]],opt_geo,atomtypes_tot[match_idx[0]]) 
        

        write_modelist("{}/frag-{}.modes".format(args.outputname,N_frag),scanned_modes,bond_mat_tot[count_i])

        # Update the list of modes that have been associated with a fragment
        run_list += match_idx        
        N_frag += 1

    # Create submission script for geometry optimization(s)
    #with open("{}/mass_geoopt.sh".format(args.outputname),'w') as f:
    #    f.write("#!/bin/bash\n"+\
    #            "\n"+\
    #            "# generate qc input files\n"+\
    #            "for i in *.xyz; do\n"+\
    #            "    name=$( echo ${i} | awk -F '.' '{print $1 }')\n"+\
    #            "    python {}".format('/'.join(os.path.abspath(__file__).split('/')[:-1])+"/xyz_to_orca.py")+" ${i} -p 8 -r 0.0 --no_freq "+"-q {}".format(args.q_tot)+" -o ${name}_geoopt.in\n"+\
    #            "    mkdir ${name}_geoopt\n"+\
    #            "    mv ${name}_geoopt.in ${name}_geoopt/.\n"+\
    #            "done\n")

    print("\n"+"*"*120)
    print("* {:^116s} *".format("Fragment generation complete!"))
    print("*"*120 + "\n")
    geo_dump = geo_dict_out
    for key in geo_dump:
      geo_dump[key]['geo'] = geo_dump[key]['geo'].tolist()
    
    if os.path.isdir(args.outputname) is False:
      print("Error: output folder {}  doesn't exist".format(args.outputname))
      quit()
    json.dump(geo_dump,codecs.open(args.outputname+'/out_intra.json', 'w', encoding='utf-8'),indent=10)
    return geo_dict_out 

# Wrapper function for write commands for *modes files
def write_modelist(name,modes,bond_mat=[]):
    with open(name,'w') as f:
        for count_i,i in enumerate(modes.keys()):
            if len(i) == 2:
                modetype='bond'
            elif len(i) == 3:
                modetype='angle'
            elif len(modes[i]["modes"][0]) == 4 and 2 in [ j[modes[i]["modes"][0][1],modes[i]["modes"][0][2]] for j in bond_mat ]:
                modetype='harmonic_dihedral'
            elif len(modes[i]["modes"][0]) == 4:
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
def mode_frag(M,Geometry,Adj_mat,Bond_mat,Elements,Atom_types,gens,q_tot=0,keep_terminal=False,keep_rings=False,keep_types=True):

    # Initialize mass_dict (used for identifying the dihedral among a coincident set that will be explicitly scanned)
    if hasattr(mode_frag,"mass_dict") is False:
        mode_frag.mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                               'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                               'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                               'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                               'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                               'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                               'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                               'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Check the consistency of the supplied arguments
    if gens < 1: print("ERROR in mode_frag: gens variable must be an integer >= 1. Exiting..."); quit();
    if keep_terminal not in [ True, False]: print("ERROR in mode_frag: keep_terminal argument must be set to a bool. Exiting..."); quit();
    if keep_rings not in [ True, False]: print("ERROR in mode_frag: keep_rings argument must be set to a bool. Exiting..."); quit();

    # Generate the fragment geometry/attribute lists corresponding to the mode
    N_M,N_Geometry,N_Adj_mat,tmp = mode_geo(M,Geometry,Adj_mat,gens,dup=[Elements,Atom_types])
    N_Elements,N_Atom_types = tmp

    # Preservation protocol depends on the mode type the fragment was generated for.
    fixed_bonds = []
    if len(N_M) == 1:
        loop_list = N_M
    elif len(N_M) == 2:
        loop_list = N_M

        # Updated fixed bonds based on the resonance structure with the largest bond order
        bmat_ind = [ i[M[0],M[1]] for i in Bond_mat ]
        bmat_ind = bmat_ind.index(max(bmat_ind))
        fixed_bonds += [(N_M[0],N_M[1],int(Bond_mat[bmat_ind][M[0],M[1]]))]
    elif len(N_M) == 3:
        loop_list = [N_M[1]]

        # Updated fixed bonds based on the resonance structure with the largest bond order of the bonds involved in the bend
        bmat_ind = [ i[M[0],M[1]]+i[M[1],M[2]] for i in Bond_mat ]
        bmat_ind = bmat_ind.index(max(bmat_ind))
        fixed_bonds += [(N_M[0],N_M[1],int(Bond_mat[bmat_ind][M[0],M[1]]))]
        fixed_bonds += [(N_M[1],N_M[2],int(Bond_mat[bmat_ind][M[1],M[2]]))]
    elif len(N_M) == 4:
        if sum([ N_Adj_mat[N_M[0],j] for j in N_M[1:] ]) == 3: 
            print("WARNING: USING IMPROPER CRITERIA ON {} CHECK THE RESULT".format([ N_Atom_types[i] for i in N_M ]))
            loop_list = N_M[0]
        else:
            loop_list = N_M[1:3]

            # Updated fixed bonds based on the resonance structure with the largest bond order
            bmat_ind = [ i[M[1],M[2]] for i in Bond_mat ]
            bmat_ind = bmat_ind.index(max(bmat_ind))            
            fixed_bonds += [(N_M[1],N_M[2],int(Bond_mat[bmat_ind][M[1],M[2]]))]
    else:
        print("ERROR in mode_frag: Protocol doesn't exist for modes involving more than 4 atoms.")
        quit()

    # Include the atoms in the mode and connected atoms within the preserve list.
    preserve = []
    gs = graph_seps(N_Adj_mat)
    for i in loop_list:
        preserve += [ count_i for count_i,i in enumerate(gs[i]) if i < gens ]
    preserve = set(preserve)

    # Perform Hydrogenation    
    N_Geometry,tmp_Atom_types,N_Elements,N_Adj_mat,added_idx = add_hydrogens(N_Geometry,N_Adj_mat,deepcopy(N_Atom_types),N_Elements,q_tot=q_tot,preserve=preserve,fixed_bonds=fixed_bonds)

    # Update the link types
    N_Atom_types += [ 'link-'+tmp_Atom_types[i] for i in added_idx ]

    # Canonicalize by sorting the elements based on hashing (NOTE: range(len(N_atom_types)) is used here rather than "atoms" as in the keep_terminal is True option. 
    Masses = [ mode_frag.mass_dict[N_Elements[i]] for i in range(len(N_Atom_types)) ]
    hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_i,N_Adj_mat,Masses),i) for count_i,i in enumerate(range(len(N_Atom_types))) ],reverse=True)) ]

    # Update lists/arrays based on the sorted atoms
    N_M          = tuple([ atoms.index(i) for i in N_M ])
    N_Geometry   = N_Geometry[atoms]
    N_Adj_mat    = N_Adj_mat[atoms]
    N_Adj_mat    = N_Adj_mat[:,atoms]
    N_Elements   = [ N_Elements[i] for i in atoms ]
    N_Atom_types = [ N_Atom_types[i] for i in atoms ]

    # If keep_types is False then the atomtypes are recalculated here
    if keep_types is False:
        N_Atom_types = id_types(N_Elements,N_Adj_mat,gens=gens,geo=N_Geometry)

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
                        FF_dict["dihedrals"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
                    elif fields[5] == "harmonic":
                        FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ] 
                    elif fields[5] == "quadratic":
                        FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(fields[6]),float(fields[7]) ]
                if fields[0].lower() == "vdw":    
                    FF_dict["vdw"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                    FF_dict["vdw"][(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]

    return FF_dict,modes_from_FF

# # This function handles remapping the atomtypes of one fragment on to another such that the modes (bonds
# # angles, dihedrals, etc) in the remapped fragment are the same after remapping. 
def map_frag(fixed_atomtypes,fixed_adj_mat,fixed_bond_mats,fixed_hash,remapped_atomtypes,remapped_adj_mat,remapped_bond_mats,remapped_hash):
    
    # The equivalency of the two topologies is based upon the mode types and number each topology possesses
    bond_types_0,angle_types_0,dihedral_types_0,one_five_types_0 = Find_modes(remapped_adj_mat,remapped_atomtypes,remapped_bond_mats,return_all=1)[4:]
    comp_obj = ( sorted(bond_types_0),sorted(angle_types_0),sorted(dihedral_types_0) )

    # Seed the mapping by placing the first atom in the mappable set onto the first instance 
    # of its hash type in the fixed set.
    #first = len(remapped_hash)-1
    first = 0
    mapping   = {first:fixed_hash.index(remapped_hash[first])}  # fixed indices keyed to remapped indices
    R_mapping = {fixed_hash.index(remapped_hash[first]):first}  # remapped indices keyed to fixed indices
    remapped_place_me = set(range(len(remapped_atomtypes)))     # remapped atom indices that need to be placed
    fixed_place_me = set(range(len(remapped_atomtypes)))        # fixed atom indices that haven't found a match in the remapped set
    remapped_place_me.remove(first)                             # remove the starting atom index
    fixed_place_me.remove(mapping[first])                       # remove the starting atom index

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
    # NOTE: the fixed_bond_mats is used for determining dihedral types
    bond_types_1,angle_types_1,dihedral_types_1,one_five_types_1 = Find_modes(fixed_adj_mat,new_atomtypes,fixed_bond_mats,return_all=1)[4:]
    if comp_obj != ( sorted(bond_types_1),sorted(angle_types_1),sorted(dihedral_types_1)):
        print("ERROR in map_frag: No mapping was discovered that preserved the topology of the original fragment. Check that the two fragments are indeed identical. Exiting...")
        print("{:60s} {:60s}".format("Original","New"))
        for count_i,i in enumerate(remapped_atomtypes):
            print("{:<3d} {:60s} {:60s}".format(count_i,i,new_atomtypes[count_i]))
        return new_atomtypes,mapping,1
    
    return new_atomtypes,mapping,0        

# Add hydrogens based upon the supplied atom types. 
# This function is only compatible with TAFFI atom types
# NOTE: Hydrogenation heuristics for geometry assume carbon behavior. This isn't usually a problem when the results are refined with transify, but more specific rules should be implemented in the future
def add_hydrogens(geo,adj_mat,atomtypes,elements,q_tot=0,preserve=[],saturate=True,retype=True,fixed_bonds=[]):

    # Initialize the saturation dictionary the first time this function is called
    if not hasattr(add_hydrogens, "sat_dict"):
        add_hydrogens.sat_dict = {  'H':1, 'He':1,\
                                   'Li':1, 'Be':2,                                                                                                                'B':3,     'C':4,     'N':3,     'O':2,     'F':1,    'Ne':1,\
                                   'Na':1, 'Mg':2,                                                                                                               'Al':3,    'Si':4,     'P':3,     'S':2,    'Cl':1,    'Ar':1,\
                                    'K':1, 'Ca':2, 'Sc':None, 'Ti':None,  'V':None, 'Cr':None, 'Mn':None, 'Fe':None, 'Co':None, 'Ni':None, 'Cu':None, 'Zn':None, 'Ga':None, 'Ge':None, 'As':None, 'Se':None, 'Br':1,    'Kr':None,\
                                   'Rb':1, 'Sr':2,  'Y':None, 'Zr':None, 'Nb':None, 'Mo':None, 'Tc':None, 'Ru':None, 'Rh':None, 'Pd':None, 'Ag':None, 'Cd':None, 'In':None, 'Sn':None, 'Sb':None, 'Te':None,  'I':1,    'Xe':None,\
                                   'Cs':1, 'Ba':2, 'La':None, 'Hf':None, 'Ta':None,  'W':None, 'Re':None, 'Os':None, 'Ir':None, 'Pt':None, 'Au':None, 'Hg':None, 'Tl':None, 'Pb':None, 'Bi':None, 'Po':None, 'At':None, 'Rn':None  }

        add_hydrogens.lone_e = {    'H':0, 'He':2,\
                                   'Li':0, 'Be':2,                                                                                                                'B':0,     'C':0,     'N':2,     'O':4,     'F':6,    'Ne':8,\
                                   'Na':0, 'Mg':2,                                                                                                               'Al':0,    'Si':0,     'P':2,     'S':4,    'Cl':6,    'Ar':8,\
                                    'K':0, 'Ca':2, 'Sc':None, 'Ti':None,  'V':None, 'Cr':None, 'Mn':None, 'Fe':None, 'Co':None, 'Ni':None, 'Cu':None, 'Zn':None, 'Ga':None, 'Ge':0,    'As':3,    'Se':4,    'Br':6,    'Kr':None,\
                                   'Rb':0, 'Sr':2,  'Y':None, 'Zr':None, 'Nb':None, 'Mo':None, 'Tc':None, 'Ru':None, 'Rh':None, 'Pd':None, 'Ag':None, 'Cd':None, 'In':None, 'Sn':None, 'Sb':None, 'Te':None,  'I':6,    'Xe':None,\
                                   'Cs':0, 'Ba':2, 'La':None, 'Hf':None, 'Ta':None,  'W':None, 'Re':None, 'Os':None, 'Ir':None, 'Pt':None, 'Au':None, 'Hg':None, 'Tl':None, 'Pb':None, 'Bi':None, 'Po':None, 'At':None, 'Rn':None  }

        add_hydrogens.frag = 0

    # Intermediate scalars
    H_length = 1.1
    N_atoms  = len(geo)
    init_len = len(geo)

    # If the user specifies a set of atoms to preserve as is, then
    # then bonding_pref entry is set to full saturation.
    if preserve != []:
        bonding_pref = [ (i,add_hydrogens.sat_dict[elements[i]]) for i in preserve ]
    else:
        bonding_pref = None

    # Get the lewis structure
    lone_electrons,bonding_electrons,core_electrons,bonding_pref = check_lewis(atomtypes,adj_mat,q_tot=q_tot,bonding_pref=bonding_pref,return_pref=True,fixed_bonds=fixed_bonds)

    # Update the preserved atoms (check_lewis will extend this list if there are special groups (e.g., nitro moieties) that need to be conserved for the sake of the lewis structure)
    preserve = set([ i[0] for i in bonding_pref ])

    # Loop over the atoms in the geometry
    for count_i,i in enumerate(geo):
        
        # ID undercoordinated atoms
        if count_i in preserve:
            continue
        elif add_hydrogens.sat_dict[elements[count_i]] is not None:
            B_expected = add_hydrogens.sat_dict[elements[count_i]]
        else:
            print("ERROR in add_hydrogens: could not determine the number of hydrogens to add to {}. Exiting...".format(elements[count_i]))
            quit()
        B_current  = bonding_electrons[count_i]

        # Determine the number of nuclei that are attached and expected.
        N_current   = sum(adj_mat[count_i])
        N_expected = N_current + (B_expected - B_current)

        # Add hydrogens to undercoordinated atoms
        if N_expected > N_current:
            
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
                    print("WARNING in add_hydrogens: 1:4 (unchecked)")
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
            # ISSUE, NEW ALGORITHM IS BASED ON BONDED ATOMS NOT BONDED CENTERS
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
        return geo,id_types(elements,adj_mat,gens=2,geo=geo),elements,adj_mat,list(range(init_len,len(geo)))
    else:
        return geo,atomtypes,elements,adj_mat,list(range(init_len,len(geo)))

# Shortcut for normalizing a vector
def normalize(x):
    return x/sum(x**(2.0))**(0.5)

# Recursively generates the permutations of indices corresponding to identical values in types
def gen_perms(types,start=0):
    ind = [ count_i for count_i,i in enumerate(types) if i==types[0] ]
    not_ind = [ count_i for count_i,i in enumerate(types) if i !=types[0] ]
    for i in permutations(list(range(start,start+len(ind))),len(ind)):
        if len(not_ind) > 0:
            for j in gen_perms([ types[j] for j in not_ind ],start=len(ind)+start):
                yield i+j
        else:
            yield i

def GetInchi(Elements,Geometry): #convert geo2inchi key

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('inchi.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    # create inchikey
    p = subprocess.Popen("obabel -ixyz inchi.xyz -oinchikey".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
    inchikey, err = p.communicate()
    inchikey = str(inchikey, 'utf-8').split('\n')[0] ##might be a problem if there are characters utf-8 can't handle
    # remove xyz file
    os.remove("inchi.xyz")
   
    return inchikey
# Logger object redirects standard output to a file.
class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/frag_gen.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

            
if __name__ == "__main__":
   main(sys.argv[1:])


