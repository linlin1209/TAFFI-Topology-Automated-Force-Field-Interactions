#!/bin/env python                                                                                                                                                             

import sys,argparse,subprocess,os,time,math
from subprocess import PIPE
import numpy as np
from math import sqrt,sin,cos,tan,factorial
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
import json, codecs
from builtins import any #python 3

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in .xyz file(s) and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode. '+\
                                                 'Model compounds are placed in individual folders named by inchi-key, and the generated file, dependency.txt, contains the dependency graph.' )

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    parser.add_argument('-FF', dest='FF_db', default='',
                        help = 'This variable holds the filename(s) of the .db files holding existing force-field parameters to be avoided in the fragment generation (default: FF.db)')

    parser.add_argument('-xyz_avoid', dest='xyz_avoid', default='',
                        help = 'This variable holds the space-delimited filename(s) of xyz files. Any atomtypes for which these xyz files are canonical parametrization fragments '+\
                               'for are avoided when generating fragments (default: None)')

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
    args=parser.parse_args(argv)    

    sys.stdout = Logger()
    print ("PROGRAM CALL: python frag_gen_inter.py {}\n".format(' '.join([ i for i in argv])))

    # Create a config dictionary based on command line options. This is passed to fun() which is an importable main.
    # XXX Reformat config definition calls. 7/18/20
    config = {} ##create a fake config 
    config['xyz'] = args.coord_files.split()
    config['ff'] = args.FF_db.split()
    config['gens'] = args.gens 
    config['xyz_avoid'] = args.xyz_avoid.split()
    config['outputname'] = args.outputname
    config['gens'] = int(args.gens)
    config['keep_terminal_opt'] = args.keep_terminal_opt
    config['q_tot'] = int(args.q_tot)
    config['avoid_frags_opt'] = args.avoid_frags_opt
    config['recursive_opt'] = args.recursive_opt

    # Generate the model compounds
    geo_dict_out = fun(config)
      
    return geo_dict_out

# Importable main function. Does the heavy lifting of generating the model compounds and dependency tree.
def fun(config):

    # Find files
    config['xyz'] = find_files(config['xyz'],recursive=config['recursive_opt'])
   
    # Perform some consistency checks
    if False in [ i.split('.')[-1] == 'xyz' for i in config['xyz'] ]:
        print( "ERROR: Check to ensure that the input file(s) are in .xyz format.")
        return
    if False in [ os.path.isfile(i) for i in config['xyz'] ]:
        print( "ERROR: Could not find file {}. Exiting...".format(next([ i for i in config['xyz'] if os.path.isfile(i) == False ])))
    for i in config['ff']:
        if os.path.isfile(i) == False:
            print ("ERROR: FF file {} couldn't not be found. Exiting...".format(i))
            return
    for i in config['xyz_avoid']:
        if os.path.isfile(i) == False:
            print ("ERROR: xyz_avoid file {} couldn't not be found. Exiting...".format(i))
            return

    # Generate a dictionary of the various modes contained in the supplied database files
    FF_db,modes_from_FF = parse_FF_params(config['ff'])
 
    # Generate model comopounds.
    mc_dict, mode_instances, instance_types = model_compound_generator(config,FF_db)

    # Generate dependency lists and assign generations. These are assigned to geo_dict, so there is no return value. 
    assign_generations(mc_dict,mode_instances)
    
    # Map redundant model compounds, write mode files associated with each model compound,
    # and write the transified geometries to dedicated folders (based on inchi key).
    geo_dict_out = generate_final(config,mc_dict,mode_instances,instance_types)

    # Does anything actually use this? If not, then remove and simplify generate_final. 7/17/20
    return geo_dict_out

# Description: recursively generates model compounds for all modes present in config['xyz'] but not in the supplied 
#              force-field files (FF_db). A new model compound is generated for each mode that is discovered. This 
#              potentially results in many equivalent model compounds being generated, however redundant models are 
#              later removed by the generate_final function. 
# Returns:
# geo_dict:       dictionary containing the  model compounds and their properties,
# mode_instances: a list of tuples holding the key in geo_dict of the correponding compound and the atom
#                 indices in the model compound for the mode (e.g., (5,(1,0,2)) would be an angle 
#                 corresponding to atoms 1, 0, and 2 in model compound geo_dict[5]. Later, assign_generations
#                 inserts the generation for parameterization into the second element.
# instance_types: a list of mode types, indexed to mode_instances. Each mode is specified as a tuple of the 
#                 atomtypes in the mode. Dihedrals have an extra term corresponding to the potential type.
def model_compound_generator(config,FF_db):

    # Generate a dictionary of the various atomtypes whose canonical fragments match the supplied xyz_avoid files
    xyz_avoid_set = parse_avoided_types(config["xyz_avoid"],config["gens"])

    # Initialize dictionary and lists for holding modes and geometry information
    geo_dict       = {}
    atoms          = []
    atom_types     = []
    atom_keys      = []
    bonds          = []
    bond_types     = []
    bond_keys      = []
    angles         = []
    angle_types    = []
    angle_keys     = []
    dihedrals      = []
    dihedral_types = []
    dihedral_keys  = []

    # Loop over the input structures, generate the all trans-minimized conformer,  and parse the modes
    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Parsing atomtypes and modes in need of parameterization for the following files"))
    print( "*"*120 + "\n")    
    for i in config["xyz"]: print( "  {}".format(i))        

    # Initialize geo_dict with the supplied coord_file information
    for count_i,i in enumerate(config["xyz"]):

        # Initialize subdictionary
        geo_dict[count_i] = {}

        # Extract elements, coords, and charge from file
        geo_dict[count_i]["elements"],geo_dict[count_i]["geo"],geo_dict[count_i]["q_tot"] = xyz_parse(i,q_opt=True)

        # Generate adjacency table
        geo_dict[count_i]["adj_mat"] = Table_generator(geo_dict[count_i]["elements"],geo_dict[count_i]["geo"])

        # Canonicalize the geometry
        geo_dict[count_i]["elements"],geo_dict[count_i]["geo"],geo_dict[count_i]["adj_mat"],geo_dict[count_i]["atom_types"] = canon_geo(geo_dict[count_i]["elements"],geo_dict[count_i]["geo"],\
                                                                                                                                        geo_dict[count_i]["adj_mat"],[0]*len(geo_dict[count_i]["elements"]))

        # Determine bonding matrix for the compound
        lone_electrons,bonding_electrons,core_electrons,geo_dict[count_i]["bond_mat"],geo_dict[count_i]["fc"] =\
            find_lewis(geo_dict[count_i]["elements"],geo_dict[count_i]["adj_mat"],q_tot=geo_dict[count_i]["q_tot"],return_pref=False,verbose=False,return_FC=True)
        
        # Find radicals for atom type assignment
        keep_lones = [ [ count_j for count_j,j in enumerate(lone_electron) if j % 2 != 0] for lone_electron in lone_electrons] 

        # Find atom_types
        geo_dict[count_i]["atom_types"] = id_types(geo_dict[count_i]["elements"],geo_dict[count_i]["adj_mat"],gens=config["gens"],fc=geo_dict[count_i]["fc"],keep_lone=keep_lones,return_index=False)

    # Find atom types in need of parameterization
    # XXX Another headache of using depth 2 angles and dihedrals is that we have to generate model compounds for the atom types first
    #     then generate the model compounds for the rest of the modes. 
    geo_dict_0 = deepcopy(geo_dict)
    geo_keys = list(geo_dict.keys())
    mc_dict = {}
    for i in geo_keys:

        # Find first instance of each atom and mode in need of parameterization
        # Generate their model compounds and propagate the information to the atom/bond/angle/dihedral lists
        # mc_dict holds the model compound information and gets updated by the helper function. 
        # For atoms
        for count_j,j in enumerate(geo_dict[i]["atom_types"]):
            if j not in atom_types:
                # Generate model
                update_mc_dict([count_j],j,geo_dict,mc_dict,config,geo_keys,i,atoms,atom_types,atom_keys)
        
        del geo_dict[i]

    # Determine atom type dependencies between model compounds
    keys = mc_dict.keys()
    for i in keys:
        tmp = []
        for j in keys:
            if any([ k in mc_dict[j]["min_types"] for k in mc_dict[i]["dep_types"]]):
                tmp += [j]
        mc_dict[i]["dep_keys"] = tmp

    # Assign the generation to the model compounds for atom types (at_mc)
    # Description: Algorithm starts all calculations out at the same generation and moves individual 
    # fragments to later generations based on their dependences. Break occurs once no generation update occur
    for i in keys: mc_dict[i]["gen"] = 1
    flag = True
    count = 0
    while flag is True:
        flag = False
        count += 1

        # Loop over each fragment and if it is in the same generation as one of its dependencies, increment its generation upward
        for i in keys:
            tmp_max = [ mc_dict[j]["gen"] for j in mc_dict[i]["dep_keys"] ]
            if len(tmp_max) > 0:
                tmp_max = max(tmp_max)
                if mc_dict[i]["gen"] <= tmp_max:
                    mc_dict[i]["gen"] = tmp_max + 1
                    flag = True

    # 1) Restore geo_dict to its original input structures,
    # 2) add a "gen" element that comes after the atom type model compounds, 
    # 3) add atomtype model compounds (current mc_dict) to restored geo_dict. 
    max_gen = max([ mc_dict[i]["gen"] for i in mc_dict.keys() ])
    geo_dict = geo_dict_0                                        # restore original dictionary of input structures.
    for i in geo_dict.keys(): geo_dict[i]["gen"] = max_gen + 1   # add gen argument
    for i in mc_dict.keys():                                     # add atomtype model compounds to geo_dict
        geo_dict[len(geo_dict.keys())] = mc_dict[i]

    # Find atoms and modes in need of parameterization starting from lowest generation up
    # This ordering is required to avoid angles and dihedrals being parameterized to parent-dependent atoms at positions more than two bonds away.
    # Model compounds for parameterization are added to geo_dict at each iteration
    # The loop only terminates when no new model compounds are generated
    atoms          = []
    atom_types     = []
    atom_keys      = []
    _,geo_keys = map(list,zip(*sorted([(geo_dict[i]["gen"],i) for i in geo_dict.keys() ])))
    mc_dict = {} # Holds the actual model compounds information 
    for i in geo_keys:

        # Find modes (returns all modes in canonicalized format)
        bonds_tmp,angles_tmp,dihedrals_tmp,one_fives_tmp,bond_types_tmp,angle_types_tmp,dihedral_types_tmp,one_five_types_tmp = Find_modes(geo_dict[i]["adj_mat"],\
                                                                                                                                           geo_dict[i]["atom_types"],\
                                                                                                                                           geo_dict[i]["bond_mat"],return_all=1)

        # Find first instance of each atom and mode in need of parameterization
        # Generate their model compounds and propagate the information to the atom/bond/angle/dihedral lists
        # mc_dict holds the model compound information and gets updated by the helper function. 
        # For atoms
        for count_j,j in enumerate(geo_dict[i]["atom_types"]):
            if j not in atom_types and j not in xyz_avoid_set and ( (j,j) not in FF_db["vdw"].keys() or j not in FF_db["charges"].keys() ):                
                # Generate model
                update_mc_dict([count_j],j,geo_dict,mc_dict,config,geo_keys,i,atoms,atom_types,atom_keys)

        # For bonds
        for count_j,j in enumerate(bond_types_tmp):
            if j not in bond_types and j not in FF_db["bonds"]:
                # Generate model
                update_mc_dict(bonds_tmp[count_j],j,geo_dict,mc_dict,config,geo_keys,i,bonds,bond_types,bond_keys)

        # For angles
        for count_j,j in enumerate(angle_types_tmp):
            if j not in angle_types and j not in FF_db["angles"]:
                # Generate model
                update_mc_dict(angles_tmp[count_j],j,geo_dict,mc_dict,config,geo_keys,i,angles,angle_types,angle_keys)

        # For dihedrals
        for count_j,j in enumerate(dihedral_types_tmp):
            if j not in dihedral_types and j not in FF_db["dihedrals"]:
                # Generate model
                update_mc_dict(dihedrals_tmp[count_j],j,geo_dict,mc_dict,config,geo_keys,i,dihedrals,dihedral_types,dihedral_keys)

        # Remove the current compound after all of its model compounds have been generated and added to mc_dict
        del geo_dict[i]

    # Define mass dictionary for determining which dihedrals are explicitly scanned
    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Collect mode_instances and instance_types
    mode_instances = [] # list of ( geo_key, mode ) tuples, where mode is a tuple of the atom indices of the mode
    instance_types = [] # list of  mode_types indexed to mode_instances
    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Force-field terms in need of parameterization"))
    print( "*"*120 + "\n")    
    if len(atom_types)>0: 
        print("    Atomtypes ({}):\n\n\t".format(len(set(atom_types)))+"\n\t".join(atom_types)+"\n")
        mode_instances += [ (atom_keys[count_i],i) for count_i,i in enumerate(atoms) ]
        instance_types += [ (i,) for i in atom_types ]
    if len(bond_types)>0: 
        print("    Bonds ({}):\n\n\t".format(len(set(bond_types)))+"\n\t".join([ "{:<40s} {:<40s}".format(*_) for _ in  bond_types])+"\n")
        mode_instances += [ (bond_keys[count_i],i) for count_i,i in enumerate(bonds) ]
        instance_types += bond_types        
    if len(angle_types)>0: 
        print("    Angles ({}):\n\n\t".format(len(set(angle_types)))+"\n\t".join([ "{:<40s} {:<40s} {:<40s}".format(*_) for _ in  angle_types])+"\n")
        mode_instances += [ (angle_keys[count_i],i) for count_i,i in enumerate(angles) ]
        instance_types += angle_types                
    if len(dihedral_types)>0: 
        print("    Dihedrals ({}) (*=implicit):\n".format(len(set(dihedral_types))))

        implicit_dihedrals = []  
        assigned_dihedrals = []
        for count_i,i in enumerate(dihedral_types):

            # Skip the dihedral if it has already been assigned
            if i in assigned_dihedrals: continue

            # Harmonic dihedrals (TAFFI treats 2-3 atoms connected by a double bond explicitly as harmonic functions. 
            if i[4] == "harmonic":
                assigned_dihedrals += [i]
                mode_instances     += [(dihedral_keys[count_i],dihedrals[count_i])]
                instance_types     += [i]
                print("\t{:<40s} {:<40s} {:<40s} {:<40s} {:<10s}".format(*i))
                continue

            # Gather coincident dihedrals (dihedrals that share the 2-3 bond with this dihedral)
            left_atoms  = [ count_k for count_k,k in enumerate(mc_dict[dihedral_keys[count_i]]["adj_mat"][dihedrals[count_i][1]]) if k == 1 and count_k != dihedrals[count_i][2] ]
            right_atoms = [ count_k for count_k,k in enumerate(mc_dict[dihedral_keys[count_i]]["adj_mat"][dihedrals[count_i][2]]) if k == 1 and count_k != dihedrals[count_i][1] ]
            coincident_dihedrals = [ tuple([k])+dihedrals[count_i][1:3]+tuple([m]) for k in left_atoms for m in right_atoms ] 
            coincident_types     = [ tuple([ mc_dict[dihedral_keys[count_i]]["atom_types"][m] for m in k ] + ["opls"]) for k in coincident_dihedrals ]

            # Canonicalize the dihedrals 
            for count_k,k in enumerate(coincident_dihedrals):
                coincident_types[count_k],coincident_dihedrals[count_k] = canon_dihedral(coincident_types[count_k],coincident_dihedrals[count_k])

            # Find the dihedral that will be explicitly scanned (the dihedral with the largest combined mass on the 1 and 4 atoms gets explicitly scanned)
            max_mass = 0.0
            for count_k,k in enumerate(coincident_dihedrals):
                current_mass = mass_dict[mc_dict[dihedral_keys[count_i]]["elements"][k[0]]] + mass_dict[mc_dict[dihedral_keys[count_i]]["elements"][k[3]]]
                if current_mass > max_mass:
                    max_mass = current_mass
                    keep = count_k

            # Add the implicitly treated dihedrals to the implicit_dihedrals list and update the mode_instances and instance_types lists with the explicit dihedral
            implicit_dihedrals += [ k for k in coincident_types if ( k != coincident_types[keep] and k != coincident_types[keep][::-1] ) ]                                
            mode_instances += [(dihedral_keys[count_i],coincident_dihedrals[keep])]   
            instance_types += [coincident_types[keep]]
            assigned_dihedrals += coincident_types
            print( "\t{:<40s} {:<40s} {:<40s} {:<40s} {:<10s}".format(*coincident_types[keep]))

        # Print the dihedrals that are being implicitly scanned
        for i in set(implicit_dihedrals):
            print( "\t{:<40s} {:<40s} {:<40s} {:<40s} {:<10s} *".format(*i))
    print( "\n\tTotal modes being explicitly parametrized: {}".format(len(mode_instances)))

    return mc_dict, mode_instances, instance_types

# Description: helper function for model_compound_generator that updates the model compound dictionary
#              based on the supplied mode and the parent compound information in geo_dict. The logic is 
#              that the model compound is generated for each mode and added to mc_dict. This ends up 
#              producing a lot of redundant model compounds, but these are later removed. The algorithm is
#              simpler this way rather than trying to generate the minimal set directly.  
def update_mc_dict(M,M_type,geo_dict,mc_dict,config,geo_keys,current,M_list,M_types,M_keys):

    # Generate model
    geo_id = max(geo_dict.keys()) + 1
    geo_keys += [geo_id]    
    geo_dict[geo_id] = {}
    geo_dict[geo_id]["q_tot"] = geo_dict[current]["q_tot"] ### XXX Needs to be updated to propagate the charge of ionic model compounds. maybe this needs to take place within mode_frag 7/10/20
    N_M,geo_dict[geo_id]["geo"],geo_dict[geo_id]["adj_mat"],geo_dict[geo_id]["elements"],geo_dict[geo_id]["atom_types"],geo_dict[geo_id]["hash_list"] = \
        mode_frag(M,geo_dict[current]["geo"],geo_dict[current]["adj_mat"],geo_dict[current]["bond_mat"],geo_dict[current]["elements"],geo_dict[current]["atom_types"],\
                  config["gens"],geo_dict[current]["q_tot"],keep_terminal=config["keep_terminal_opt"],keep_types=True,force_linear=True)

    # Determine bonding matrix for the model compound
    lone_electrons,bonding_electrons,core_electrons,geo_dict[geo_id]["bond_mat"],geo_dict[geo_id]["fc"] =\
            find_lewis(geo_dict[geo_id]["elements"],geo_dict[geo_id]["adj_mat"],q_tot=geo_dict[geo_id]["q_tot"],return_pref=False,verbose=False,return_FC=True)

    # Find radicals for atom typ assignment
    keep_lones = [ [ count_i for count_i,i in enumerate(lone_electron) if i % 2 != 0] for lone_electron in lone_electrons] 

    # Save model compound and original atom types to mc_dict
    mc_id = len(mc_dict.keys())
    mc_dict[mc_id] = { _:geo_dict[geo_id][_] for _ in geo_dict[geo_id].keys() }
    M_list += [N_M]
    M_types += [M_type]
    M_keys += [mc_id]   

    # Determine true atomtypes for the model compound (note, we are using keep_types = True for mode_frag). geo_dict needs these for
    # correct mode assignment in subsequent iterations. And atomtype model compounds also need these in mc_dict.
    geo_dict[geo_id]["atom_types"] = id_types(geo_dict[geo_id]["elements"],geo_dict[geo_id]["adj_mat"],gens=config["gens"],fc=geo_dict[geo_id]["fc"],keep_lone=keep_lones,return_index=False)

    # Because angles and dihedrals use d=2 model compounds, calculated with respect to the central atoms, they are parameterized with the parent compound atom types
    # however, bond and atom type model compounds need to use the actual atomtypes corresponding to the model compound. XXX We will revisit this in the future by potentially
    # using depth=1 atom types for angle and dihedral fitting.
    if len(M) in [1,2]:
        mc_dict[mc_id]["atom_types"] = deepcopy(geo_dict[geo_id]["atom_types"])
        
    # Find the atom types for which this is a model compound (empty if this is a model for an intramolecular mode) 
    if len(M) == 1:
        mc_dict[mc_id]["min_types"] = set([ j for j in mc_dict[mc_id]["atom_types"] if minimal_structure(j,elements=mc_dict[mc_id]["elements"],adj_mat=mc_dict[mc_id]["adj_mat"],\
                                                                                                         atomtypes=mc_dict[mc_id]["atom_types"],gens=config["gens"]) is True ])
    else:
        mc_dict[mc_id]["min_types"] = set([])

    # Find the atom types that this mode depends on
    if len(M) == 1:
        mc_dict[mc_id]["dep_types"] = set([ j.replace("link-","") for j in mc_dict[mc_id]["atom_types"] if minimal_structure(j,elements=mc_dict[mc_id]["elements"],adj_mat=mc_dict[mc_id]["adj_mat"],\
                                                                                                                             atomtypes=mc_dict[mc_id]["atom_types"],gens=config["gens"]) is False ])
    else:
        mc_dict[mc_id]["dep_types"] = set([ j.replace("link-","") for j in mc_dict[mc_id]["atom_types"] ])
    return

# Description: determines the dependencies among the model compounds based on the atom types. 
#              "dep_keys" holds a list of model compound indices that each model compound depends on
#              "gen" holds the generation that the model compound belongs to. The generation info 
#              also added to each mode_instance as the second element of the tuple.
# Returns:     None, data is added to the supplied mc_dict dictionary. 
# XXX If this breaks, we may need to Revise in future to sort dependence based on atoms, then bonds, then angles, then dihedrals.  
def assign_generations(mc_dict,mode_instances):

    # Find the model compounds corresponding to atom types
    at_mc = [ i[0] for i in mode_instances if len(i[1]) == 1 ]
    others = [ _ for _ in mc_dict.keys() if _ not in at_mc ]

    # Determine atom type dependencies between model compounds
    keys = mc_dict.keys()
    for i in keys:
        tmp = []
        for j in keys:
            if any([ k in mc_dict[j]["min_types"] for k in mc_dict[i]["dep_types"]]):
                tmp += [j]
        mc_dict[i]["dep_keys"] = tmp

    # Assign the generation to the model compounds for atom types (at_mc)
    # Description: Algorithm starts all calculations out at the same generation and moves individual 
    # fragments to later generations based on their dependences. Break occurs once no generation update occur
    for i in at_mc: mc_dict[i]["gen"] = 1
    flag = True
    count = 0
    while flag is True:
        flag = False
        count += 1

        # Loop over each fragment and if it is in the same generation as one of its dependencies, increment its generation upward
        for i in at_mc:
            tmp_max = [ mc_dict[j]["gen"] for j in mc_dict[i]["dep_keys"] ]
            if len(tmp_max) > 0:
                tmp_max = max(tmp_max)
                if mc_dict[i]["gen"] <= tmp_max:
                    mc_dict[i]["gen"] = tmp_max + 1
                    flag = True

        # Check for an unreasonable number of generations
        if count >= 1000:
            print("There appears to be a mutual dependence among model compounds. assign_generations failed after 1000 attempts.")
            quit()

    # Assign the generation to the model compounds for modes (others)
    for i in others:
        mc_dict[i]["gen"] = max([ mc_dict[j]["gen"] for j in mc_dict[i]["dep_keys"] ])

    # Add generation information to each mode_instance as the second element
    for count_i,i in enumerate(mode_instances):
        mode_instances[count_i] = (i[0],mc_dict[i[0]]["gen"],i[1])

    return



# Description: this function handles mapping redundant model compounds, writing the mode files, and outputting
#              the json dictionary. Mapping is necessary because a model compound is generated for every mode
#              by the model_compoound_generator function. Many of these model compounds are equivalent, and the
#              map_frag function here takes care of determining equivalency and making sure that the mode indices
#              are aligned. 
#
# Returns:     mc_dict_out: a dictionary with the final model compound geometry atomtypes etc.
def generate_final(config,mc_dict,mode_instances,instance_types):

    #### Instead of output a lot of folders and xyz file, now we simply save all the ouput into a dict: mc_dict_out
    #### similar to one used in frag_gen, but now mc_dict_out[file], the filename is now inchikey
    #### which contains: elements,geometry,atom_types,gen ; so that we can directly input it into frag_gen.py
    #### Main difference: has generation
    _,keys = zip(*sorted([ (mc_dict[i]["gen"],i) for i in mc_dict.keys() ]))

    mc_dict_out = {} 

    # Save one instance of each unique fragment and a frag-*.modes file to keep track of what modes to parametrize for each fragment
    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Beginning canonical fragment generation"))
    print( "*"+"-"*118+"*")
    print( "* {:<116s} *".format("NOTE: Each atomtype and bond are parametrized using the smallest molecular fragment that is consistent with its"))
    print( "* {:<116s} *".format("      topological definition. For angles and dihedrals, the smallest molecular fragment is used that is consistent"))
    print( "* {:<116s} *".format("      with their central atom, and bond, respectively."))
    print( "*"*120 + "\n")
    run_list = []
    N_frag = 0
    for i in keys:

        # Skip already assigned fragments
        if i in run_list: continue        

        # Find equivalent fragments based on a comparison of the hash lists
        match_idx = []
        for j in keys:
            if mc_dict[i]["hash_list"] == mc_dict[j]["hash_list"]:

                # Add the current match, and remap the atomtypes and mode indices
                # Note: map_frag does the heavy lifting of mapping the atom indices in j to i (graph isomorphism problem!) the idea is that multiple modes use the same model compound
                # so only one is subjected to explicit calculations, but the distinct atomtypes for fitting each mode are retained in teh *.modes file.
                match_idx += [j]
                mc_dict[j]["atom_types"],mapping,status = map_frag(mc_dict[i]["atom_types"],mc_dict[i]["adj_mat"],mc_dict[i]["bond_mat"],mc_dict[i]["hash_list"],\
                                                                    mc_dict[j]["atom_types"],mc_dict[j]["adj_mat"],mc_dict[j]["bond_mat"],mc_dict[j]["hash_list"])
                mode_instances = [ (m[0],m[1],tuple([ mapping[k] for k in m[2] ])) if m[0] == j else m for m in mode_instances ]

                # Print problem fragment(s) in the case of failure status being returned by map_frag
                if status == 1:
                    write_xyz("incomp_fixed",mc_dict[i]["elements"],mc_dict[i]["geo"],mc_dict[i]["atom_types"])                    
                    write_xyz("incomp_mapped",mc_dict[j]["elements"],mc_dict[j]["geo"],mc_dict[j]["atom_types"])
                    quit()

                # add this model compound key to run list because it is being assigned here
                run_list += [j]

        # Generate inchi key for this model compound
        opt_geo = transify(mc_dict[i]["geo"],mc_dict[i]["adj_mat"],elements=mc_dict[i]["elements"],opt_terminals=True,opt_final=True)
        inchikey = GetInchi(mc_dict[i]["elements"],opt_geo)

        # mode_inds
        mode_inds = [ count_j for count_j,j in enumerate(mode_instances) if j[0] in match_idx ]
        
        # Assemble the dictionary of modes being scanned and print diagnostics for the user
        # NOTE: dihedrals are treated differently from bonds and angles because in some cases two dihedrals have identical fragments and simply
        #       differ in the atomtypes/VDW/charges used when fitting the potential. 
        print( "\n\tModes associated with model compound ({:<s}) (inchi: {}) (gen in leading parentheses):".format(str(N_frag),inchikey))
        scanned_modes = {}
        placed = []
        for j in [ k for k in mode_inds if len(mode_instances[k][2]) == 1]:
            scanned_modes[instance_types[j]] = { "gens":[mode_instances[j][1]], "modes":[mode_instances[j][2]], "atomtypes":[mc_dict[mode_instances[j][0]]["atom_types"]] }
            print( "\t\t({})    {}".format(mode_instances[j][1]," ".join([ "{:40s}".format(m) for m in instance_types[j] ])))
        for j in [ k for k in mode_inds if len(mode_instances[k][2]) == 2]:
            scanned_modes[instance_types[j]] = { "gens":[mode_instances[j][1]], "modes":[mode_instances[j][2]], "atomtypes":[mc_dict[mode_instances[j][0]]["atom_types"]] }
            print( "\t\t({})    {}".format(mode_instances[j][1]," ".join([ "{:40s}".format(m) for m in instance_types[j] ])))
        for j in [ k for k in mode_inds if len(mode_instances[k][2]) == 3]:
            scanned_modes[instance_types[j]] = { "gens":[mode_instances[j][1]], "modes":[mode_instances[j][2]], "atomtypes":[mc_dict[mode_instances[j][0]]["atom_types"]] }
            print( "\t\t({})    {}".format(mode_instances[j][1]," ".join([ "{:40s}".format(m) for m in instance_types[j] ])))
        for j in [ k for k in mode_inds if len(mode_instances[k][2]) == 4]:
            print( "\t\t({})    {}".format(mode_instances[j][1]," ".join([ "{:40s}".format(m) for m in instance_types[j] ])))
            if j in placed: continue
            scanned_modes[instance_types[j]] = { "gens":[], "modes":[], "atomtypes":[] }

            # If j is a harmonic dihedral then redundancy is assessed on the basis of all four mode 
            # indices matching (this is necessary beacuse of the way that harmonic dihedrals are scanned)
            if 2 in [ n[mode_instances[j][2][1],mode_instances[j][2][2]] for n in mc_dict[mode_instances[j][0]]["bond_mat"] ]:

                # Iterate over all dihedrals that have identical atomic indices in their definition
                for m in [ n for n in mode_inds if len(mode_instances[n][2]) == 4]:

                    # Save dihedrals that have indentical atomic indices
                    if mode_instances[j][2] == mode_instances[m][2] or mode_instances[j][2] == mode_instances[m][2][::-1]:
                        placed += [m]
                        scanned_modes[instance_types[j]]["gens"]      += [mode_instances[m][1]]
                        scanned_modes[instance_types[j]]["modes"]     += [mode_instances[m][2]]
                        scanned_modes[instance_types[j]]["atomtypes"] += [mc_dict[mode_instances[m][0]]["atom_types"]]

            # If j is a fourier dihedral then redundancy is assessed on the basic of the 2-3 atoms in the dihedral
            else:

                # Iterate over all dihedrals that have identical fragments
                for m in [ n for n in mode_inds if len(mode_instances[n][2]) == 4]:

                    # Save dihedrals that have identical 2-3 atoms in their definition
                    if mode_instances[j][2][1:3] == mode_instances[m][2][1:3] or mode_instances[j][2][1:3] == mode_instances[m][2][1:3][::-1]:
                        placed += [m]
                        scanned_modes[instance_types[j]]["gens"]      += [mode_instances[m][1]]
                        scanned_modes[instance_types[j]]["modes"]     += [mode_instances[m][2]]
                        scanned_modes[instance_types[j]]["atomtypes"] += [mc_dict[mode_instances[m][0]]["atom_types"]]

        # Find real atomtypes for this model compound (needed for subsequent parsers)
        lone_electrons,bonding_electrons,core_electrons,bond_mat,fc =\
            find_lewis(mc_dict[i]["elements"],mc_dict[i]["adj_mat"],q_tot=mc_dict[i]["q_tot"],return_pref=False,verbose=False,return_FC=True)
        keep_lones = [ [ count_j for count_j,j in enumerate(lone_electron) if j % 2 != 0] for lone_electron in lone_electrons] 
        model_atomtypes = id_types(mc_dict[i]["elements"],mc_dict[i]["adj_mat"],gens=config["gens"],fc=mc_dict[i]["fc"],keep_lone=keep_lones,return_index=False)        

        # If mc_dict_out is not needed for the return, then these objects should just be passed to the json dictionary and not go through this intermediate. 7/17/20
        mc_dict_out[inchikey] = {}
        mc_dict_out[inchikey]["gen"] = sorted(set([ k for j in scanned_modes.keys() for k in scanned_modes[j]["gens"] ]))
        mc_dict_out[inchikey]["atom_types"] = model_atomtypes
        mc_dict_out[inchikey]["q"] = mc_dict[i]["q_tot"]

        # Make inchi_folder
        if os.path.isdir("{}".format(inchikey)) is False:
            os.makedirs("{}".format(inchikey))

        # Write the frag-*.modes file associated with this fragment
        write_modelist("{}/{}.modes".format(inchikey,inchikey),scanned_modes,mc_dict[i]["bond_mat"])
        write_xyz("{}/{}".format(inchikey,inchikey),mc_dict[i]["elements"],opt_geo,model_atomtypes)

        # Update the list of modes that have been associated with a fragment
        run_list += match_idx        
        N_frag += 1

    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Fragment generation complete!"))
    print( "*"*120 + "\n")

    # dump all.json for the taffi driver to utilize for running the parameterization jobs
    flag_keywords = ["geoopt","modescan"]
    all = { _:False for _ in flag_keywords }
    all["gens"] = {i:{ _ : False for _ in ["initial_fit","charges","after_charge","vdw","final_fit"] } for i in set( k for j in mc_dict_out.keys() for k in mc_dict_out[j]["gen"] ) }
    all["mc"] = { i:{"gen":mc_dict_out[i]["gen"], "atom_types":mc_dict_out[i]["atom_types"], "q":mc_dict_out[i]["q"] } for i in mc_dict_out.keys() }
    json.dump(all,codecs.open('all.json', 'w', encoding='utf-8'),indent=10)    
    
    return mc_dict_out



# Wrapper function for write commands for *modes files
def write_modelist(name,modes,bond_mat=[]):
    with open(name,'w') as f:
        for count_i,i in enumerate(modes.keys()):
            if len(i) == 1:
                modetype='atom'                
            elif len(i) == 2:
                modetype='bond'
            elif len(i) == 3:
                modetype='angle'
            elif len(modes[i]["modes"][0]) == 4 and 2 in [ j[modes[i]["modes"][0][1],modes[i]["modes"][0][2]] for j in bond_mat ]:
                modetype='harmonic_dihedral'
            elif len(modes[i]["modes"][0]) == 4:
                modetype='dihedral'                
            f.write("\n{} start\n".format(modetype))
            f.write('{}\n'.format(len(modes[i]["modes"])))
            f.write('{}\n'.format(" ".join([ "{:<60s}".format(str(j)) for j in modes[i]["gens"] ])))
            f.write("{}\n".format(" ".join([ "{:<60s}".format("_".join([ str(k) for k in j])) for j in modes[i]["modes"] ])))
            for count_j,j in enumerate(modes[i]["atomtypes"][0]):
                f.write('{}\n'.format(" ".join([ "{:<60s}".format(modes[i]["atomtypes"][k][count_j]) for k in range(len(modes[i]["atomtypes"])) ])))
            f.write("{} end\n".format(modetype))
    return 


# Returns the geometry and relevant property lists corresponding to the
# smallest fragment that is consistent with the mode being parametrized.
# If the 1 or 4 atom is a ring, then the entire ring is included. 
# XXX Doesn't return the charge state of the new model compound 7/10/20
def mode_frag(M,Geometry,Adj_mat,Bond_mat,Elements,Atom_types,gens,q_tot=0,fc=[],keep_lone=[],keep_terminal=False,keep_rings=False,keep_types=True,force_linear=False,return_FC=False):
    
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
    if gens < 1: print( "ERROR in mode_frag: gens variable must be an integer >= 1. Exiting..."); quit();
    if keep_terminal not in [ True, False]: print( "ERROR in mode_frag: keep_terminal argument must be set to a bool. Exiting..."); quit();
    if keep_rings not in [ True, False]: print( "ERROR in mode_frag: keep_rings argument must be set to a bool. Exiting..."); quit();

    # Generate the fragment geometry/attribute lists corresponding to the mode
    N_M,N_Geometry,N_Adj_mat,N_Bond_mat,tmp = new_mode_geo(M,Geometry,Adj_mat,Bond_mat,gens,dup=[Elements,Atom_types],force_linear=force_linear)
    N_Elements,N_Atom_types = tmp

    # determine loop list and fixed_bonds
    # loops over resonance structures and only fixes bonds if consistent bonding occurs across structures.
    # loop_list contains the atoms in the mode whose bond structure is fixed (if consistent across resonance structures)
    total_fixed_bonds = []
    for lb,bond_mat in enumerate(N_Bond_mat):

        fixed_bonds = []
        if len(M) == 1:
            loop_list = N_M

        elif len(M) == 2:
            loop_list = N_M
            fixed_bonds += [(N_M[0],N_M[1],int(bond_mat[N_M[0],N_M[1]]))]

        elif len(N_M) == 3:
            loop_list = [N_M[1]]
            fixed_bonds += [(N_M[0],N_M[1],int(bond_mat[N_M[0],N_M[1]]))] # XXX I think these should be removed, redundant with loop_list? 7/10/20
            fixed_bonds += [(N_M[1],N_M[2],int(bond_mat[N_M[1],N_M[2]]))]

        elif len(N_M) == 4:
            if sum([ N_Adj_mat[N_M[0],j] for j in N_M[1:] ]) == 3: 
                print("WARNING: USING IMPROPER DIHEDRAL CRITERIA ON {} CHECK THE RESULT".format([ N_Atom_types[i] for i in N_M ]))
                loop_list = N_M[0]
            else:
                loop_list = N_M[1:3]
                fixed_bonds += [(N_M[1],N_M[2],int(bond_mat[N_M[1],N_M[2]]))]
        else:
            print("ERROR in mode_frag: Protocol doesn't exist for modes involving more than 4 atoms.")
            quit() 

        # Include the atoms in the mode and connected atoms within the preserve list.    
        for i in loop_list:     
            fixed_bonds += [(i,j,int(k)) for j,k in enumerate(bond_mat[i]) if k > 1]

        total_fixed_bonds += [fixed_bonds]

    # only if all fixed_bonds are same, take it as fixed_bonds; else, fixed_bonds=[]
    if len(list(map(list,set(map(tuple,total_fixed_bonds))))) > 1: 
        fixed_bonds = []
    else:
        fixed_bonds = total_fixed_bonds[0] 

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
        N_Atom_types = id_types(N_Elements,N_Adj_mat,gens=gens)

    return N_M,N_Geometry,N_Adj_mat,N_Elements,N_Atom_types,hash_list

# Description: This is a simple implementation of the Dijkstra algorithm for 
#              finding the backbone of a polymer 
def Dijkstra(Adj_mat):

    # Remove terminal sites (sites with only a single length 2
    # self walk). Continue until all of the terminal structure 
    # has been removed from the topology.
    Adj_trimmed = np.copy(Adj_mat)

    # Initialize Distances, Previous, and Visited lists    
    Distances = np.array([0] + [100000]*(len(Adj_mat)-1)) # Holds shortest distance to origin from each site
    Previous = np.array([0]*len(Adj_mat)) # Holds the previous site on the short distance to origin
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

# # Adds parameters from the FF file(s) to the FF_dict.
# def parse_FF_params(FF_files,FF_dict={"masses":{},"charges":{},"bonds":{},"angles":{},"dihedrals":{},"dihedrals_harmonic":{},"vdw":{}}):
                   
#     modes_from_FF = []
#     for i in FF_files:
#         with open(i,'r') as f:
#             for lines in f:
#                 fields = lines.split()
#                 if len(fields) == 0: continue
#                 if fields[0].lower() == "atom":   FF_dict["masses"][fields[1]] = float(fields[3])
#                 if fields[0].lower() == "charge": FF_dict["charges"][fields[1]] = float(fields[2])
#                 if fields[0].lower() == "bond":   
#                     modes_from_FF += [(fields[1],fields[2])]
#                     modes_from_FF += [(fields[2],fields[1])]
#                     FF_dict["bonds"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
#                 if fields[0].lower() == "angle":
#                     modes_from_FF += [(fields[1],fields[2],fields[3])]
#                     modes_from_FF += [(fields[3],fields[2],fields[1])]
#                     FF_dict["angles"][(fields[1],fields[2],fields[3])] = [fields[4],float(fields[5]),float(fields[6])]
#                 if fields[0].lower() in ["dihedral","torsion"]: 
#                     modes_from_FF += [(fields[1],fields[2],fields[3],fields[4])]
#                     modes_from_FF += [(fields[4],fields[3],fields[2],fields[1])]
#                     if fields[5] == "opls":       
#                         FF_dict["dihedrals"][(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
#                     elif fields[5] == "harmonic":
#                         FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ] 
#                     elif fields[5] == "quadratic":
#                         FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),float(fields[7]) ]
#                 if fields[0].lower() == "vdw":    
#                     FF_dict["vdw"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
#                     FF_dict["vdw"][(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]

#     return FF_dict,modes_from_FF

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
        print( "ERROR in map_frag: No mapping was discovered that preserved the topology of the original fragment. Check that the two fragments are indeed identical. Exiting...")
        print( "{:60s} {:60s}".format("Original","New"))
        for count_i,i in enumerate(remapped_atomtypes):
            print( "{:<3d} {:60s} {:60s}".format(count_i,i,new_atomtypes[count_i]))
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
        bonding_pref = []

    # Get the lewis structure
    #lone_electrons,bonding_electrons,core_electrons,bonding_pref = check_lewis(atomtypes,adj_mat,q_tot=q_tot,bonding_pref=bonding_pref,return_pref=True,fixed_bonds=fixed_bonds)
    # since the input doesn't contain fc infor yet, set fc all be zero
    fc = [0] * len(elements)
    lone_electrons,bonding_electrons,core_electrons,bonding_pref = frag_find_lewis(elements,adj_mat,q_tot=q_tot,fc_0=fc,keep_lone=[],fixed_bonds=fixed_bonds,bonding_pref=bonding_pref,return_pref=True,check_lewis_flag=True)

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
            print( "ERROR in add_hydrogens: could not determine the number of hydrogens to add to {}. Exiting...".format(elements[count_i]))
            quit()
        B_current  = bonding_electrons[count_i]

        # Determine the number of nuclei that are attached and expected.
        N_current   = np.sum(adj_mat[count_i])
        N_expected = N_current + (B_expected - B_current)

        # Add hydrogens to undercoordinated atoms
        if N_expected > N_current:

            old_inds = [ count_j for count_j,j in enumerate(adj_mat[count_i]) if j == 1 ]

            # Protocols for 1 missing hydrogen
            if N_expected - N_current == 1:
                if N_expected == 1:
                    new = i + np.array([H_length,0.0,0.0])
                elif N_expected == 2:
                    new = -1.0 * normalize(geo[old_inds[0]] - i) * H_length + i + np.array([random.random(),random.random(),random.random()])*0.01 #random factor added for non-carbon types to relax during FF-opt
                elif N_expected == 3:
                    new = -1.0 * normalize( normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) ) * H_length + i
                elif N_expected == 4:
                    print( "WARNING in add_hydrogens: 1:4 (unchecked)")
                    new = -1.0 * normalize( normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) + normalize(geo[old_inds[2]] - i) ) * H_length + i                

                # Update geometry, adj_mat, elements, and atomtypes with one new atoms
                geo = np.vstack([geo,new])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]
                elements += ["H"]
                tmp = np.zeros([N_atoms+1,N_atoms+1])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[-1,count_i] = 1
                tmp[count_i,-1] = 1
                adj_mat = tmp                
                N_atoms += 1

            # Protocols for 2 missing hydrogens
            # ISSUE, NEW ALGORITHM IS BASED ON BONDED ATOMS NOT BONDED CENTERS
            if N_expected - N_current == 2:
                if N_expected == 2:
                    new_1 = i + np.array([H_length,0.0,0.0])
                    new_2 = i - np.array([H_length,0.0,0.0])
                elif N_expected == 3:
                    rot_vec = normalize(np.cross( geo[old_inds[0]] - i, np.array([random.random(),random.random(),random.random()]) ))
                    new_1 = normalize(axis_rot(geo[old_inds[0]],rot_vec,i,120.0) - i)*H_length + i
                    new_2 = normalize(axis_rot(geo[old_inds[0]],rot_vec,i,240.0) - i)*H_length + i
                elif N_expected == 4:
                    bisector = normalize(geo[old_inds[0]] - i + geo[old_inds[1]] - i) 
                    new_1    = axis_rot(geo[old_inds[0]],bisector,i,90.0)
                    new_2    = axis_rot(geo[old_inds[1]],bisector,i,90.0) 
                    rot_vec  = normalize(np.cross(new_1-i,new_2-i))
                    angle    = ( 109.5 - acos(np.dot(normalize(new_1-i),normalize(new_2-i)))*180.0/pi ) / 2.0
                    new_1    = axis_rot(new_1,rot_vec,i,-angle)
                    new_2    = axis_rot(new_2,rot_vec,i,angle)
                    new_1    = -1*H_length*normalize(new_1-i) + i
                    new_2    = -1*H_length*normalize(new_2-i) + i
                    
                # Update geometry, adj_mat, elements, and atomtypes with two new atoms
                geo = np.vstack([geo,new_1])
                geo = np.vstack([geo,new_2])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*2
                elements += ["H","H"]
                tmp = np.zeros([N_atoms+2,N_atoms+2])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2],count_i] = 1
                tmp[count_i,[-1,-2]] = 1
                adj_mat = tmp
                N_atoms += 2

            # Protocols for 3 missing hydrogens
            if N_expected - N_current == 3:
                if N_expected == 3:
                    rot_vec = np.array([0.0,1.0,0.0])
                    new_1 = i + np.array([H_length,0.0,0.0])
                    new_2 = axis_rot(new_1,rot_vec,i,120.0)
                    new_3 = axis_rot(new_1,rot_vec,i,240.0)
                if N_expected == 4:
                    rot_vec = normalize(np.cross( geo[old_inds[0]] - i, np.array([random.random(),random.random(),random.random()]) ))
                    new_1 = H_length*normalize(axis_rot(geo[old_inds[0]],rot_vec,i,109.5)-i) + i
                    new_2 = axis_rot(new_1,normalize(i-geo[old_inds[0]]),i,120.0)
                    new_3 = axis_rot(new_2,normalize(i-geo[old_inds[0]]),i,120.0)

                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = np.vstack([geo,new_1])
                geo = np.vstack([geo,new_2])
                geo = np.vstack([geo,new_3])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*3
                elements += ["H","H","H"]
                tmp = np.zeros([N_atoms+3,N_atoms+3])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2,-3],count_i] = 1
                tmp[count_i,[-1,-2,-3]] = 1
                adj_mat = tmp
                N_atoms += 3

            # Protocols for 4 missing hydrogens
            if N_expected - N_current == 4:
                if N_expected == 4:
                    new_1 = i + np.array([H_length,0.0,0.0])
                    rot_vec = normalize(np.cross( new_1 - i, np.array([random.random(),random.random(),random.random()]) ))
                    new_2 = H_length*normalize(axis_rot(new_1,rot_vec,i,109.5)-i) + i
                    new_3 = axis_rot(new_2,normalize(i-new_1),i,120.0)
                    new_4 = axis_rot(new_3,normalize(i-new_1),i,120.0)
                    
                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = np.vstack([geo,new_1])
                geo = np.vstack([geo,new_2])
                geo = np.vstack([geo,new_3])
                geo = np.vstack([geo,new_4])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*4
                elements += ["H","H","H","H"]
                tmp = np.zeros([N_atoms+4,N_atoms+4])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2,-3,-4],count_i] = 1
                tmp[count_i,[-1,-2,-3,-4]] = 1
                adj_mat = tmp
                N_atoms += 4

    if retype is True:
        atom_types=id_types(elements,adj_mat,gens=2)
        return geo,atom_types,elements,adj_mat,range(init_len,len(geo))
    else:
        return geo,atomtypes,elements,adj_mat,range(init_len,len(geo))

# Shortcut for normalizing a vector
def normalize(x):
    return x/np.sum(x**(2.0))**(0.5)

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

# Description:   Returns the list of atomtypes whose canonical fragments are included in
#                the xyz_avoid geometries. 
#
# Inputs:        xyz_avoid:     A list of xyz filenames
#                gens:          The bond depth for determining atom uniqueness
#
# Outputs:       xyz_avoid_set: A set of atomtypes whose canonical fragments are 
#                               contained in xyz_avoid.
def parse_avoided_types(xyz_avoid,gens):

    xyz_avoid_set = []

    for files in xyz_avoid:

        # Extract Element list and Coord list from the file
        elements,geo = xyz_parse(files)

        # Generate adjacency table
        adj_mat = Table_generator(elements,geo)

        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        hybridizations = Hybridization_finder(elements,adj_mat)

        # Find atom types
        #atom_types = id_types(elements,adj_mat,gens,hybridizations,geo)
        atom_types = id_types(elements,adj_mat,gens)

        for types in set(atom_types):
            
            if minimal_structure(types,geo,elements,adj_mat=adj_mat,atomtypes=atom_types,gens=gens) is True:
                xyz_avoid_set += [types]
                
    return set(xyz_avoid_set)

# Description:   Checks if the supplied geometry corresponds to the minimal structure of the molecule
# 
# Inputs:        atomtype:      The taffi atomtype being checked
#                geo:           Geometry of the molecule
#                elements:      elements, indexed to the geometry 
#                adj_mat:       adj_mat, indexed to the geometry (optional)
#                atomtypes:     atomtypes, indexed to the geometry (optional)
#                gens:          number of generations for determining atomtypes (optional, only used if atomtypes are not supplied)
# 
# Outputs:       Boolean:       True if geo is the minimal structure for the atomtype, False if not.
def minimal_structure(atomtype,geo=None,elements=None,adj_mat=None,atomtypes=None,gens=2):

    # If required find the atomtypes for the geometry
    if atomtypes is None or adj_mat is None:
        if  geo is None or ( len(elements) != len(geo) ):
            print( "ERROR in minimal_structure: While trying to automatically assign atomtypes, the elements argument must have dimensions equal to geo. Exiting...")
            quit()

        # Generate the adjacency matrix
        # NOTE: the units are converted back angstroms
        adj_mat = Table_generator(elements,geo)

        # Generate the atomtypes
        #atomtypes = id_types(elements,adj_mat,gens,Hybridization_finder(elements,adj_mat),geo)
        atomtypes = id_types(elements,adj_mat,gens)

    # Check if this is a ring type, if not and if there are rings
    # in this geometry then it is not a minimal structure. 
    if "R" not in atomtype:
        if True in [ "R" in i for i in atomtypes ]:
            return False
        
    # Check minimal conditions
    count = 0
    for count_i,i in enumerate(atomtypes):

        # If the current atomtype matches the atomtype being searched for then proceed with minimal geo check
        if i == atomtype:
            count += 1

            # Initialize lists for holding indices in the structure within "gens" bonds of the seed atom (count_i)
            keep_list = [count_i]
            new_list  = [count_i]
            
            # Carry out a "gens" bond deep search
            for j in range(gens):

                # Id atoms in the next generation
                tmp_new_list = []                
                for k in new_list:
                    tmp_new_list += [ count_m for count_m,m in enumerate(adj_mat[k]) if m == 1 and count_m not in keep_list ]

                # Update lists
                tmp_new_list = list(set(tmp_new_list))
                if len(tmp_new_list) > 0:
                    keep_list += tmp_new_list
                new_list = tmp_new_list
            
            # Check for the minimal condition
            keep_list = set(keep_list)
            if False in [ elements[j] == "H" for j in range(len(elements)) if j not in keep_list ]:
                minimal_flag = False
            else:
                minimal_flag = True
        
    return minimal_flag

# Logger object redirects standard output to a file.
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("frag_gen_all.log", "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass
            
if __name__ == "__main__":
   main(sys.argv[1:])


