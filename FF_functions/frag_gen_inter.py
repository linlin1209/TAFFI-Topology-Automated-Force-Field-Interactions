#!/bin/env python                                                                                                                                                             

import sys,argparse,subprocess,os,time,math
from subprocess import PIPE
import numpy as np#from numpy import *
from math import sqrt,sin,cos,tan,factorial
#from scipy import *
#from numpy.linalg import *
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
#from __builtin__ import any  # any has a namespace conflict with numpy. XXX WE REALLY NEED TO CLEAN UP THESE NAMESPACE ISSUES
from builtins import any #python 3

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a .xyz file and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode.')

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

    # Create a fake config for function input
    config = {} ##create a fake config 
    config['ff'] = args.FF_db
    config['gens'] = args.gens 
    
    # Convert coords_files and N into lists
    #args.FF_db = args.FF_db.split()
    #args.xyz_avoid = args.xyz_avoid.split()

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
    
    config['xyz'] = args.coord_files




    # Perform some consistency checks
    #if False in [ i.split('.')[-1] == 'xyz' for i in args.coord_files ]:
    #    print( "ERROR: Check to ensure that the input file(s) are in .xyz format.")
    #    return
    #if False in [ os.path.isfile(i) for i in args.coord_files ]:
    #    print( "ERROR: Could not find file {}. Exiting...".format(next([ i for i in args.coord_files if os.path.isfile(i) == False ])))
    #for i in args.FF_db:
    #    if os.path.isfile(i) == False:
    #        print ("ERROR: FF file {} couldn't not be found. Exiting...".format(i))
    #        return
    #for i in args.xyz_avoid:
    #    if os.path.isfile(i) == False:
    #        print ("ERROR: xyz_avoid file {} couldn't not be found. Exiting...".format(i))
    #        return

    # Generate base filename for naming output
    #Filename = args.outputname

    # Make directory to hold the output    
    #if os.path.exists(Filename):
    #    print ('The desired output folder ({}) already exists. Exiting to avoid overwriting data...'.format(Filename))
    #    return
    #else:
    #    os.makedirs(Filename)
    #    Filename = Filename+'/'+Filename
    #    sys.stdout = Logger(args.outputname)
    #    print ("PROGRAM CALL: python frag_gen_inter.py {}\n".format(' '.join([ i for i in argv])))

    # Generate a dictionary of the various modes contained in the supplied database files
    #FF_db,modes_from_FF = parse_FF_params(args.FF_db)

    # Generate a dictionary of the various atomtypes whose canonical fragments match the supplied xyz_avoid files
    #xyz_avoid_set = parse_avoided_types(args.xyz_avoid,args.gens)
      
    geo_dict_out = fun(config)
      
    return geo_dict_out



def fun(config):

    #parser = argparse.ArgumentParser(description='Reads in a .xyz file and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode.')

    #required (positional) arguments                                                                                                  
    #parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    #parser.add_argument('-FF', dest='FF_db', default='',
    #                    help = 'This variable holds the filename(s) of the .db files holding existing force-field parameters to be avoided in the fragment generation (default: FF.db)')

    #parser.add_argument('-xyz_avoid', dest='xyz_avoid', default='',
    #                    help = 'This variable holds the space-delimited filename(s) of xyz files. Any atomtypes for which these xyz files are canonical parametrization fragments '+\
    #                           'for are avoided when generating fragments (default: None)')

    #parser.add_argument('-o', dest='outputname', default='frags',
    #                    help = 'Controls the output folder name for the fragments. (default: frags)')

    #parser.add_argument('-gens', dest='gens', default=2,
    #                    help = 'Controls the bond search depth for identifying unique atom types (default: 2)')

    #parser.add_argument('-R', dest='recursive_opt', default=False, action='store_const',const=True,
    #                    help = 'When this flag is present, filenames and wildcards will be searched for recursively.')

    #parser.add_argument('-q', dest='q_tot', default=0, 
    #                    help = 'The total charge on the generated fragments. (default: 0)')

    #parser.add_argument('--keep_terminal', dest='keep_terminal_opt', default=False, action='store_const',const=True,
    #                    help = 'When this flag is present, terminal atoms at the edges of fragments are retained. Default is to perform a full hydrogenation consistent with the hybridization of the atomtype. (default: off)')

    #parser.add_argument('--avoid_frags', dest='avoid_frags_opt', default=False, action='store_const',const=True,
    #                    help = 'When this flag is present, the xyzs are used directly for the parameterization (default: False)')

    # Make relevant inputs lowercase
    #parser.parse_)    
    
   
    #default setting
    FF_db = ''
    xyz_avoid = ''
    outputname = 'frags'
    gens = 2   
    recursive_opt = False
    q_tot = 0
    keep_terminal_opt = False
    avoid_frags_opt = False
    
    #assign from config file
    FF_db = config['ff']
    gens = config['gens']
    outputname = 'Intermolecular_Modes' ##need to be deleted, we don't output any folder now
    
    # Convert coords_files and N into lists
    FF_db = FF_db.split()
    xyz_avoid = xyz_avoid.split()

    # Parse additional inputs 
    gens = int(gens)
    q_tot = int(q_tot)

    # Find files
    #coord_files = coord_files.split()
    #wild_card_files  = [ i for i in coord_files if "*" in i ]
    #if recursive_opt == True:
    #    coord_files = [ i for i in coord_files if "*" not in i ]
    #    coord_files += [ dp+"/"+files for i in wild_card_files for dp,dn,fn in os.walk('.') for files in fn if fnmatch(files,i) ] 
    #else:
    #    coord_files = [ i for i in coord_files if "*" not in i ]
    #    for i in wild_card_files:
    #        path = '/'.join(i.split('/')[:-1])
    #        if len(path) == 0:
    #            path = "."
    #        coord_files += [ path+"/"+files for files in os.listdir(path) if fnmatch(files,i) ] 

    # Handle "./" condition
    #for count_i,i in enumerate(coord_files):
    #    if i[:2] == "./": coord_files[count_i] = coord_files[count_i][2:]    
    #coord_files = list(set(coord_files))
    coord_files = config['xyz'] 
   
    # Perform some consistency checks
    if False in [ i.split('.')[-1] == 'xyz' for i in coord_files ]:
        print( "ERROR: Check to ensure that the input file(s) are in .xyz format.")
        return
    if False in [ os.path.isfile(i) for i in coord_files ]:
        print( "ERROR: Could not find file {}. Exiting...".format(next([ i for i in coord_files if os.path.isfile(i) == False ])))
    for i in FF_db:
        if os.path.isfile(i) == False:
            print ("ERROR: FF file {} couldn't not be found. Exiting...".format(i))
            return
    for i in xyz_avoid:
        if os.path.isfile(i) == False:
            print ("ERROR: xyz_avoid file {} couldn't not be found. Exiting...".format(i))
            return

    # Generate base filename for naming output
    #Filename = outputname

    # Make directory to hold the output    
    #if os.path.exists(Filename):
    #    print ('The desired output folder ({}) already exists. Exiting to avoid overwriting data...'.format(Filename))
    #    return
    #else:
    #    os.makedirs(Filename)
    #    Filename = Filename+'/'+Filename
    #    sys.stdout = Logger(outputname)
    #    print ("PROGRAM CALL: python frag_gen_inter.py {}\n".format(' '.join([ i for i in argv])))

    # Generate a dictionary of the various modes contained in the supplied database files
    FF_db,modes_from_FF = parse_FF_params(FF_db)

    # Generate a dictionary of the various atomtypes whose canonical fragments match the supplied xyz_avoid files
    xyz_avoid_set = parse_avoided_types(xyz_avoid,gens)

    # Initialize dictionary and lists for holding modes and geometry information
    geo_dict      = {}
    fit_atomtypes = []
    fit_atominds  = []
    fit_keys      = []

    # Loop over the input structures, generate the all trans-minimized conformer,  and parse the modes
    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Parsing atomtypes in need of LJ/partial charge parametrization for the following files"))
    print( "*"*120 + "\n")    
    for i in coord_files: print( "  {}".format(i))        

    # Populate geo_dict with supplied coord_file information
    for count_i,i in enumerate(coord_files):

        # Initialize subdictionary
        geo_dict[count_i] = {}

        # Extract Element list and Coord list from the file i
        geo_dict[count_i]["elements"],geo_dict[count_i]["geo"] = xyz_parse(i)
        
        # Extract total charge
        #geo_dict[count_i]["q_tot"] = parse_q(i)
        geo_dict[count_i]["q_tot"] = 0

        # Generate adjacency table
        geo_dict[count_i]["adj_mat"] = Table_generator(geo_dict[count_i]["elements"],geo_dict[count_i]["geo"])

        # Canonicalize the geometry
        geo_dict[count_i]["elements"],geo_dict[count_i]["geo"],geo_dict[count_i]["adj_mat"],geo_dict[count_i]["atom_types"] = canon_geo(geo_dict[count_i]["elements"],geo_dict[count_i]["geo"],\
                                                                                                                                        geo_dict[count_i]["adj_mat"],[0]*len(geo_dict[count_i]["elements"]))

        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        hybridizations = Hybridization_finder(geo_dict[count_i]["elements"],geo_dict[count_i]["adj_mat"])

        # Determine bonding matrix for the compound
        lone_electrons,bonding_electrons,core_electrons,geo_dict[count_i]["bond_mat"],geo_dict[count_i]["fc"] =\
            find_lewis(geo_dict[count_i]["elements"],geo_dict[count_i]["adj_mat"],q_tot=geo_dict[count_i]["q_tot"],return_pref=False,verbose=False,return_FC=True)
        #geo_dict[count_i]["bond_mat"] = find_lewis(geo_dict[count_i]["atom_types"], geo_dict[count_i]["adj_mat"], q_tot=q_tot,  b_mat_only=True,verbose=False)
        
        # Find atom types
        #geo_dict[count_i]["atom_types"] = id_types(geo_dict[count_i]["elements"],geo_dict[count_i]["adj_mat"],gens,hybridizations,geo_dict[count_i]["geo"])
        keep_lones = [ [ count_j for count_j,j in enumerate(lone_electron) if j % 2 != 0] for lone_electron in lone_electrons] 
        geo_dict[count_i]["atom_types"] = id_types(geo_dict[count_i]["elements"],geo_dict[count_i]["adj_mat"],gens=gens,fc=geo_dict[count_i]["fc"],keep_lone=keep_lones,return_index=False)

        # Find atomtypes in need of VDW and/or charge parametrization
        tmp_atomtypes = [ (count_j,j) for count_j,j in enumerate(geo_dict[count_i]["atom_types"]) \
                          if j not in fit_atomtypes and j not in xyz_avoid_set \
                          and ( (j,j) not in FF_db["vdw"].keys() or j not in FF_db["charges"].keys() ) ]

        # Remove duplicates
        keep_ind = []
        for j in set([ k[1] for k in tmp_atomtypes ]):
            keep_ind += [next( count_m for count_m,m in enumerate(tmp_atomtypes) if m[1] == j )]

        # Update the master lists
        fit_atomtypes += [ tmp_atomtypes[j][1] for j in keep_ind ]
        fit_atominds  += [ tmp_atomtypes[j][0] for j in keep_ind ]
        fit_keys += [ count_i for j in keep_ind ]

        # Add the linear analogs of any ring atom types
        if avoid_frags_opt is False:
            if True in [ "R" in j for j in geo_dict[count_i]["atom_types"] ]:

                # Assemble a list of ring atom indices whose linear mode_frag will be generated then whose atomtypes will be included in the list. 
                modes = [ [count_j] for count_j,j in enumerate(geo_dict[count_i]["atom_types"]) if "R" in j ]

                for count_j,j in enumerate(modes):

                    # Get the linear analog of the molecule about the ring atoms or the ring bonds (j)
                    geo_dict["{}-{}-lin".format(count_i,count_j)] = {}
                    _,geo_dict["{}-{}-lin".format(count_i,count_j)]["geo"],geo_dict["{}-{}-lin".format(count_i,count_j)]["adj_mat"],geo_dict["{}-{}-lin".format(count_i,count_j)]["elements"],\
                    geo_dict["{}-{}-lin".format(count_i,count_j)]["atom_types"],_ = \
                    mode_frag(j,geo_dict[count_i]["geo"],geo_dict[count_i]["adj_mat"],geo_dict[count_i]["bond_mat"],geo_dict[count_i]["elements"],geo_dict[count_i]["atom_types"],\
                              gens,q_tot,keep_terminal=keep_terminal_opt,keep_types=False,force_linear=True)
                    
                    # Canonicalize the geometry
                    geo_dict["{}-{}-lin".format(count_i,count_j)]["elements"],geo_dict["{}-{}-lin".format(count_i,count_j)]["geo"],\
                    geo_dict["{}-{}-lin".format(count_i,count_j)]["adj_mat"],geo_dict["{}-{}-lin".format(count_i,count_j)]["atom_types"] = canon_geo(geo_dict["{}-{}-lin".format(count_i,count_j)]["elements"],\
                                                                                                                                                     geo_dict["{}-{}-lin".format(count_i,count_j)]["geo"],\
                                                                                                                                                     geo_dict["{}-{}-lin".format(count_i,count_j)]["adj_mat"],\
                                                                                                                                                     geo_dict["{}-{}-lin".format(count_i,count_j)]["atom_types"])

                    # Determine the bonding matrix for the linear analog
                    geo_dict["{}-{}-lin".format(count_i,count_j)]["bond_mat"] = find_lewis(geo_dict["{}-{}-lin".format(count_i,count_j)]["elements"], \
                                                                                           geo_dict["{}-{}-lin".format(count_i,count_j)]["adj_mat"], q_tot=q_tot,  b_mat_only=True,verbose=False)                

                    # Find atomtypes in need of VDW and/or charge parametrization
                    tmp_atomtypes = [ (count_k,k) for count_k,k in enumerate(geo_dict["{}-{}-lin".format(count_i,count_j)]["atom_types"]) \
                                      if k not in fit_atomtypes and k not in xyz_avoid_set \
                                      and ( (k,k) not in FF_db["vdw"].keys() or k not in FF_db["charges"].keys() ) ]

                    # Remove duplicates
                    keep_ind = []
                    for n in set([ k[1] for k in tmp_atomtypes ]):
                        keep_ind += [next( count_m for count_m,m in enumerate(tmp_atomtypes) if m[1] == n )]

                    # Update the master lists
                    fit_atomtypes += [ tmp_atomtypes[k][1] for k in keep_ind ]
                    fit_atominds  += [ tmp_atomtypes[k][0] for k in keep_ind ]
                    fit_keys += [ "{}-{}-lin".format(count_i,count_j) for k in keep_ind ]

    # Define mass dictionary for avoid_frags == True case
    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Generate fragments for each atomtype being parametrized with additional fragments being
    # appended to handle any missing dependent atom types.
    modes_tot      = []
    modetypes_tot  = []
    atomtypes_tot  = []
    geos_tot       = []
    elements_tot   = []
    adj_mat_tot    = []
    bond_mat_tot   = []
    hash_lists_tot = []    
    min_types_tot  = []
    dep_types_tot  = []
    for count_i,i in enumerate(fit_atomtypes):

        # Use the original xyz
        if avoid_frags_opt:
            mode_ind        = [fit_atominds[count_i]]
            frag_geo        = geo_dict[fit_keys[count_i]]["geo"]
            frag_adj_mat    = geo_dict[fit_keys[count_i]]["adj_mat"]
            frag_bond_mat   = geo_dict[fit_keys[count_i]]["bond_mat"]
            frag_elements   = geo_dict[fit_keys[count_i]]["elements"]
            frag_atom_types = geo_dict[fit_keys[count_i]]["atom_types"]
            frag_masses     = [ mass_dict[frag_elements[j]] for j in range(len(frag_atom_types)) ]
            frag_hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_k,frag_adj_mat,frag_masses),k) for count_k,k in enumerate(range(len(frag_atom_types))) ],reverse=True)) ]            

            # Update lists/arrays based on the sorted atoms
            mode_ind        = tuple([ atoms.index(j) for j in mode_ind ])
            frag_geo        = frag_geo[atoms]
            frag_adj_mat    = frag_adj_mat[atoms]
            frag_adj_mat    = frag_adj_mat[:,atoms]
            frag_elements   = [ frag_elements[j] for j in atoms ]
            frag_atom_types = [ frag_atom_types[j] for j in atoms ]

        # Generate canonical fragment for atomtype i
        else:
            mode_ind,frag_geo,frag_adj_mat,frag_elements,frag_atom_types,frag_hash_list = \
            mode_frag([fit_atominds[count_i]],geo_dict[fit_keys[count_i]]["geo"],geo_dict[fit_keys[count_i]]["adj_mat"],geo_dict[fit_keys[count_i]]["bond_mat"],geo_dict[fit_keys[count_i]]["elements"],\
                      geo_dict[fit_keys[count_i]]["atom_types"],gens,q_tot,keep_terminal=keep_terminal_opt,keep_types=False)

            # Determine bonding matrix for the compound
            #frag_bond_mat = find_lewis(frag_atom_types, frag_adj_mat, q_tot=q_tot, b_mat_only=True,verbose=False)
            frag_bond_mat = find_lewis(frag_elements, frag_adj_mat, q_tot=q_tot, b_mat_only=True,verbose=False)

        # Build total lists        
        modes_tot      += [mode_ind]  
        modetypes_tot  += [(i,)]
        atomtypes_tot  += [frag_atom_types]
        geos_tot       += [frag_geo]
        elements_tot   += [frag_elements]
        adj_mat_tot    += [frag_adj_mat]
        bond_mat_tot   += [frag_bond_mat]
        hash_lists_tot += [frag_hash_list]
        min_types_tot  += [set([ j for j in frag_atom_types if minimal_structure(j,frag_geo,frag_elements,adj_mat=frag_adj_mat,atomtypes=frag_atom_types,gens=gens) is True ])]
        dep_types_tot  += [set([ j for j in frag_atom_types if minimal_structure(j,frag_geo,frag_elements,adj_mat=frag_adj_mat,atomtypes=frag_atom_types,gens=gens) is False ])]

        # Add the current fragment to geo_dict
        new_key = len(geo_dict.keys())
        geo_dict[new_key] = {}
        geo_dict[new_key]["elements"]   = frag_elements
        geo_dict[new_key]["geo"]        = frag_geo
        geo_dict[new_key]["adj_mat"]    = frag_adj_mat
        geo_dict[new_key]["bond_mat"]   = frag_bond_mat
        geo_dict[new_key]["atom_types"] = frag_atom_types

        # Find atomtypes in need of VDW and/or charge parametrization in the new fragment
        tmp_atomtypes = [ (count_j,j) for count_j,j in enumerate(geo_dict[new_key]["atom_types"]) \
                          if j not in fit_atomtypes and j not in xyz_avoid_set \
                          and ( (j,j) not in FF_db["vdw"].keys() or j not in FF_db["charges"].keys() ) ]

        # Remove duplicates
        keep_ind = []
        for j in set([ k[1] for k in tmp_atomtypes ]):
            keep_ind += [next( count_m for count_m,m in enumerate(tmp_atomtypes) if m[1] == j )]

        # Update the master lists
        fit_atomtypes += [ tmp_atomtypes[j][1] for j in keep_ind ]
        fit_atominds  += [ tmp_atomtypes[j][0] for j in keep_ind ]
        fit_keys += [ new_key for j in keep_ind ]

    mode_instances = []
    instance_types = []    

    # Atom type collection
    print( "\n\tAtomtypes in need of parametrization ({}):\n".format(len(set(fit_atomtypes))))
    for i in set(fit_atomtypes):
        print( "\t\t{}".format(i))

    # Assign dependence lists
    dep_lists_tot = []
    for count_i,i in enumerate(dep_types_tot):
        tmp = []
        for count_j,j in enumerate(min_types_tot):
            if any([ k in j for k in i ]):
                tmp += [count_j]
        dep_lists_tot += [tmp]
    
    # Assign the generation to each fragment
    # Description: Algorithm starts all calculations out at the same generation and moves individual 
    # fragments to later generations based on their dependences. Break occurs once no generation update occur
    flag = True
    gen_list_tot = [ 1 for i in hash_lists_tot ]
    while flag is True:
        flag = False

        # Loop over each fragment's generation and check if it needs to be updated  
        for count_i,i in enumerate(gen_list_tot):

            # Check the largest generation for a dependence
            tmp_max = [ gen_list_tot[k] for k in dep_lists_tot[count_i] ]            
            if len(tmp_max) > 0:
                tmp_max = max(tmp_max)            
                if i <= tmp_max:
                    gen_list_tot[count_i] = tmp_max + 1
                    flag = True

    # Make gen-* folders
    #gen_list = sorted(set(gen_list_tot))
    #for i in gen_list:
    #    if os.path.isdir("{}/gen-{}".format(outputname,i)) is False:
    #        os.makedirs("{}/gen-{}".format(outputname,i))
    
    #### Instead of output a lot of folders and xyz file, now we simply save all the ouput into a dict: geo_dict_out
    #### similar to one used in frag_gen, but now geo_dict_out[file], the filename is now inchikey
    #### which contains: elements,geometry,atom_types,gen ; so that we can directly input it into frag_gen.py
    #### Main difference: has generation
    geo_dict_out = {} 

    # Save one instance of each unique fragment and a frag-*.modes file to keep track of what modes to parametrize for each fragment
    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Beginning canonical fragment generation"))
    print( "*"+"-"*118+"*")
    print( "* {:<116s} *".format("NOTE: Each atomtype will be parametrized using the smallest molecular fragment that is consistent with its"))
    print( "* {:<116s} *".format("      topological definition"))
    print( "*"*120 + "\n")
    run_list = []
    N_frag = 0
    for count_i,i in enumerate(hash_lists_tot): #loop over generation

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
                ##need modified
                if status == 1:
                    write_xyz("incomp_fixed",elements_tot[count_i],geos_tot[count_i],atomtypes_tot[count_i])                    
                    write_xyz("incomp_mapped",elements_tot[count_j],geos_tot[count_j],atomtypes_tot[count_j])                    
                    quit()

        # Make directory if necessary
        #if os.path.isdir("{}/gen-{}".format(outputname,gen_list_tot[count_i])) is False:
        #    os.makedirs("{}/gen-{}".format(outputname,gen_list_tot[count_i]))
        
        # Save the final transified geometry in dict
        opt_geo = transify(geos_tot[match_idx[0]],adj_mat_tot[match_idx[0]],elements=elements_tot[match_idx[0]],opt_terminals=True,opt_final=True)
        #write_xyz("{}/gen-{}/frag-{}".format(outputname,gen_list_tot[count_i],N_frag),elements_tot[match_idx[0]],opt_geo,atomtypes_tot[match_idx[0]])
        inchikey = GetInchi(elements_tot[match_idx[0]],opt_geo)
        geo_dict_out[inchikey] = {}
        geo_dict_out[inchikey]["gen"] = gen_list_tot[count_i]
        geo_dict_out[inchikey]["elements"] = elements_tot[match_idx[0]]
        geo_dict_out[inchikey]["geo"] = opt_geo
        geo_dict_out[inchikey]["atom_types"] = atomtypes_tot[match_idx[0]]
        
        # Make inchi directory to hold the output    
        if os.path.isdir("{}".format(inchikey)) is False:
            os.makedirs(inchikey)
        sys.stdout = Logger(inchikey) ##Logger have to be put after the inchikey created, but this way we lose some output, may need modified (or not?) since we don't use these files
        write_xyz("{}/out_inter".format(inchikey),elements_tot[match_idx[0]],opt_geo,atomtypes_tot[match_idx[0]])
        
        
        print( "\n\tLennard-Jones/Partial charges asssociated with fragment {:<s} (inchi: {})(gen: {}):".format(str(N_frag),inchikey,gen_list_tot[count_i]))

        # Assemble the dictionary of modes being scanned and print diagnostics for the user
        # NOTE: dihedrals are treated differently from bonds and angles because in some cases two dihedrals have identical fragments and simply
        #       differ in the atomtypes/VDW/charges used when fitting the potential. 
        scanned_modes = {}
        placed = []
        for j in [ k for k in match_idx if len(modes_tot[k]) == 1]:
            scanned_modes[modetypes_tot[j]] = { "modes":[modes_tot[j]], "atomtypes":[atomtypes_tot[j]] }
            print( "\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
        for j in [ k for k in match_idx if len(modes_tot[k]) == 2]:
            scanned_modes[modetypes_tot[j]] = { "modes":[modes_tot[j]], "atomtypes":[atomtypes_tot[j]] }
            print( "\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
        for j in [ k for k in match_idx if len(modes_tot[k]) == 3]:
            scanned_modes[modetypes_tot[j]] = { "modes":[modes_tot[j]], "atomtypes":[atomtypes_tot[j]] }
            print( "\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
        for j in [ k for k in match_idx if len(modes_tot[k]) == 4]:
            print( "\t\t{}".format(" ".join([ "{:40s}".format(m) for m in modetypes_tot[j] ])))     
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
        #write_modelist("{}/gen-{}/frag-{}.modes".format(outputname,gen_list_tot[count_i],N_frag),scanned_modes,bond_mat_tot[count_i])

        # Update the list of modes that have been associated with a fragment
        run_list += match_idx        
        N_frag += 1

    # Create submission scripts for geometry optimization(s)
    #for i in gen_list:
    #    with open("{}/gen-{}/mass_geoopt.sh".format(outputname,i),'w') as f:
    #        f.write("#!/bin/bash\n"+\
    #                "\n"+\
    #                "# generate qc input files\n"+\
    #                "for i in *.xyz; do\n"+\
    #                "    name=$( echo ${i} | awk -F '.' '{print $1 }')\n"+\
    #                "    python {}".format('/'.join(os.path.abspath(__file__).split('/')[:-1])+"/xyz_to_orca.py")+" ${i} -p 8 -r 0.0 --no_freq "+"-q {}".format(q_tot)+" -o ${name}_geoopt.in\n"+\
    #                "    mkdir ${name}_geoopt\n"+\
    #                "    mv ${name}_geoopt.in ${name}_geoopt/.\n"+\
    #                "done\n")

    #for i in gen_list:
    #    with open("{}/gen-{}/mass_intramolecular.sh".format(outputname,i),'w') as f:
    #        f.write("#!/bin/bash\n"+\
    #                "\n"+\
    #                "# generate intramolecular parametrization files\n"+\
    #                "python {}".format('/'.join(os.path.abspath(__file__).split('/')[:-1])+"/frag_gen.py")+" '*.xyz' -o Intramolecular_Modes \n")

    print( "\n"+"*"*120)
    print( "* {:^116s} *".format("Fragment generation complete!"))
    print( "*"*120 + "\n")
    #geo_dict_out= {k: v for k, v in sorted(geo_dict_out.items(), key=lambda item: item[1]["gen"])}
    # json can't dump np array, so create one that is json dumpable
    geo_dump = geo_dict_out
    for key in geo_dump:
      geo_dump[key]['geo'] = geo_dump[key]['geo'].tolist()
    
    json.dump(geo_dump,codecs.open('out_inter.json', 'w', encoding='utf-8'),indent=10)
    #obj_text = codecs.open('./out_inter.json', 'r').read()
    #jload = json.loads(obj_text)
    #for key in jload:
    #  jload[key]['geo'] = np.np.array(jload[key]['geo'])
    #print(jload)
    
    
    return geo_dict_out

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
            f.write("{}\n".format(" ".join([ "{:<60s}".format("_".join([ str(k) for k in j])) for j in modes[i]["modes"] ])))
            for count_j,j in enumerate(modes[i]["atomtypes"][0]):
                f.write('{}\n'.format(" ".join([ "{:<60s}".format(modes[i]["atomtypes"][k][count_j]) for k in range(len(modes[i]["atomtypes"])) ])))
            f.write("{} end\n".format(modetype))
    return 


# Returns the geometry and relevant property lists corresponding to the
# smallest fragment that is consistent with the mode being parametrized.
# If the 1 or 4 atom is a ring, then the entire ring is included. 
# TO DO: Doesn't return the charge state of the new model compound 7/10/20
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
    #N_M,N_Geometry,N_Adj_mat,tmp = mode_geo(M,Geometry,Adj_mat,gens,dup=[Elements,Atom_types],force_linear=force_linear)
    #N_Elements,N_Atom_types = tmp
    N_M,N_Geometry,N_Adj_mat,N_Bond_mat,tmp = new_mode_geo(M,Geometry,Adj_mat,Bond_mat,gens,dup=[Elements,Atom_types],force_linear=force_linear)
    N_Elements,N_Atom_types = tmp

    # determine loop list and fixed_bonds
    total_fixed_bonds = []
    for lb,bond_mat in enumerate(N_Bond_mat):
        fixed_bonds = []
        if len(N_M) == 1:
            loop_list = N_M

        elif len(N_M) == 2:
            loop_list = N_M

            # Updated fixed bonds based on the resonance structure with the largest bond order
            bmat_ind = [ i[M[0],M[1]] for i in bond_mat ]
            bmat_ind = bmat_ind.index(max(bmat_ind))
            fixed_bonds += [(N_M[0],N_M[1],int(bond_mat[bmat_ind][M[0],M[1]]))]

        elif len(N_M) == 3:
            loop_list = [N_M[1]]
            # Updated fixed bonds based on the resonance structure with the largest bond order of the bonds involved in the bend
            bmat_ind = [ i[M[0],M[1]]+i[M[1],M[2]] for i in bond_mat ]
            bmat_ind = bmat_ind.index(max(bmat_ind))
            fixed_bonds += [(N_M[0],N_M[1],int(bond_mat[bmat_ind][M[0],M[1]]))]
            fixed_bonds += [(N_M[1],N_M[2],int(bond_mat[bmat_ind][M[1],M[2]]))]

        elif len(N_M) == 4:
            if sum([ N_Adj_mat[N_M[0],j] for j in N_M[1:] ]) == 3: 
                print("WARNING: USING IMPROPER CRITERIA ON {} CHECK THE RESULT".format([ N_Atom_types[i] for i in N_M ]))
                loop_list = N_M[0]
            else:
                loop_list = N_M[1:3]

            # Updated fixed bonds based on the resonance structure with the largest bond order
            bmat_ind = [ i[M[1],M[2]] for i in bond_mat ]
            bmat_ind = bmat_ind.index(max(bmat_ind))            
            fixed_bonds += [(N_M[1],N_M[2],int(bond_mat[bmat_ind][M[1],M[2]]))]
        else:
            print("ERROR in mode_frag: Protocol doesn't exist for modes involving more than 4 atoms.")
            quit() 

        # Include the atoms in the mode and connected atoms within the preserve list.    
        for i in loop_list:     
            fixed_bonds += [(i,j,int(k)) for j,k in enumerate(bond_mat[i]) if k > 1]

        total_fixed_bonds += [fixed_bonds]

    # only if all fixed_bonds are same, take it as fixed_bonds; else, fixed_bonds=[]
    if len(list(map(list,set(map(tuple,total_fixed_bonds))))) > 1: fixed_bonds = []
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
        #N_Atom_types = id_types(N_Elements,N_Adj_mat,gens=gens,geo=N_Geometry)
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
def minimal_structure(atomtype,geo,elements,adj_mat=None,atomtypes=None,gens=2):

    # If required find the atomtypes for the geometry
    if atomtypes is None or adj_mat is None:
        if len(elements) != len(geo):
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
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/frag_gen_inter.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass
            
if __name__ == "__main__":
   main(sys.argv[1:])


