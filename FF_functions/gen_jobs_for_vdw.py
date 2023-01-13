#!/bin/env python                                                                                                                                                             

import sys,argparse,os,time,ast,random,re,fnmatch,matplotlib
matplotlib.use('Agg') # Needed for cluster image generation 
from pylab import *
from numpy import *
from numpy.linalg import norm
from scipy.spatial.distance import cdist
from shutil import move,copyfile
from copy import deepcopy

def main(argv):
    
    parser = argparse.ArgumentParser(description='Extracts pair-wise configurations from a LAMMPS trajectory and generates an input for an orca binding energy calculation.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('traj_files', help = 'A lammps trajectory file (usually outputted from a vdw_self_gen.py job). Alternatively, if a folder is supplied, the program will attempt to automatically carry out the vdw cycling protocol. '+\
                                             'In the latter case, the program will search for a configs folder and md-* folders to automatically identify which trajectories are in need of sampling.')

    #optional arguments
    parser.add_argument('-map', dest="map_files", default="", help = 'The map file holding the correspondence between the atoms in the data file and the atom types in the VDW file '+\
                                                                     '(usually outputted from vdw_self_gen.py). By default, the map file will be initialized with base name of the folder '
                                                                     'containing the traj_file. In the case of folder-based vdw cycling, the map files are initialized based on the md-* folders that are present.')

    parser.add_argument('-start', dest='frame_start', default=0,
                        help = 'Frame to start parsing from in each trajectory (first=0; default: 0)')

    parser.add_argument('-end', dest='frame_end', default=0,
                        help = 'Frame to stop parsing from the trajectory (zero-indexing, like the start option; 0 defaults to last frame;  default: last)')

    parser.add_argument('-every', dest='every', default=1,
                        help = 'Every N frames will be parsed, where N is this variable (i.e. if mod(frame_num,every) == 0 then it is parsed) (parse all = 1; default: 1)')

    parser.add_argument('-N', dest='N_configs', default=100,
                        help = 'Holds the number of configurations to be generated (for each pair type) before exiting (default: 100)')

    parser.add_argument('-q_A', dest='q_A', default = 0,
                        help = 'Holds the charge of molecules in the A_list. (elementary charge; default: 0)')

    parser.add_argument('-q_B', dest='q_B', default = 0,
                        help = 'Holds the charge of molecules in the B_list. (elementary charge; default: 0)')

    parser.add_argument('-m_A', dest='m_A', default = 1,
                        help = 'Holds the multiplicity of molecules in the A_list. (default: 1)')

    parser.add_argument('-m_B', dest='m_B', default = 1,
                        help = 'Holds the multiplicity of molecules in the B_list. (default: 1)')

    parser.add_argument('-A_list', dest='A_list', default=[],
                        help = 'Holds a list of mol types (see --print_mol), based on the trajectory file, included for pair-wise configurations. For example, if x is present in both A_list and B_list, then x-x type '+\
                               'configurations will be generated. If x is present only in either A_list or B_list, then no x-x type configurations will be generated. If x is present '+\
                               'in neither list, then x is treated as a spectator and no x-containing configurations will be generated. Entries can be added using the : mini-language for '+\
                               'ranges of molecule indices, where 2:4 means 2,3,4 and 2:2:8 means 2,4,6,8 (by default both lists are populated with all mol ids)')

    parser.add_argument('-B_list', dest='B_list', default=[],
                        help = 'Holds a list of mol types (see --print_mol), based on the trajectory file, included for pair-wise configurations. For example, if x is present in both A_list and B_list, then x-x type '+\
                               'configurations will be generated. If x is present only in either A_list or B_list, then no x-x type configurations will be generated. If x is present '+\
                               'in neither list, then x is treated as a spectator and no x-containing configurations will be generated. Entries can be added using the : mini-language for '+\
                               'ranges of molecule indices, where 2:4 means 2,3,4 and 2:2:8 means 2,4,6,8 (by default both lists are populated with all mol ids)')

    parser.add_argument('-r_min_scale',dest='r_min_scale', default=1.0,
                        help = 'Configurations are drawn based on if the pair-wise separations lie within a threshold. This threshold is based on the lj r_min of the pair potential multiplied by this scaling factor. '+\
                               'An attempt is made to parse the pair potentials from the *.in.settings file. If a settings file is not discovered in the job path or parameters are missing, then the program defaults to UFF parameters. '+\
                               'The separation threshold can be made larger by setting this value > 1.0, or more restrictive by setting '+\
                               'this threshold to < 1.0. (default: 1.0)')
    
    parser.add_argument('-p', dest='procs', default=8,
                        help = 'holds the number of processors to run each orca job on (default: 8)')

    parser.add_argument('-o', dest='filename', default="vdw_orca_jobs",
                        help = 'holds the name of the folder where the logfile and output data are stored. (default: vdw_orca_jobs)')

    parser.add_argument('--write_unwrapped', dest='write_unwrapped', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the unwrapped coordinates of each parsed frame are saved. Useful for debugging. (default: off)')

    parser.add_argument('-QC_types', dest='QC_types', default='dft mp2',
                        help = 'This variable controls whether dft, mp2, or both types of jobs are initialized. The variables expects a case-insensitive '+\
                               'space-delimited string (default: "dft mp2")')

    parser.add_argument('-FF', dest='FF_db', default='',
                        help = 'This variable holds the filename of a parameters database. When supplied, the program will avoid drawing configurations for any '+\
                               'vdw pairs that already exist in the database. Multiple files can be supplied as a space-delimited string. (default: None)')

    parser.add_argument('--print_mol', dest='print_mol', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the unique molecule info is parsed from the simulation data and printed to screen. The molecule ids printed here can then be ' +\
                               'used for setting the A and B list arguments. (default: off)')
    
    parser.add_argument('-mode', dest='o_mode', default=0,
                        help = 'This variable selects the mode which determines which output files are generated. 0 = ORCA. 1 = SAPT. (default: 0)')
    
    parser.add_argument('-SAPT_name', dest='SAPT_name', default='complex',
                        help = 'This variable sets the name of the complex for a SAPT calculation.')
    
    parser.add_argument('-SAPT_metod', dest='SAPT_method', default='ssapt0',
                        help = 'This variable sets the SAPT method to use. ssapt0, sapt0, sapt2 (default: ssapt0).')
    
    parser.add_argument('-f', dest='functional', default='wB97X-D3',
                        help = 'Sets the functional for the calculation. (default: wB97X-D3; other typical options are B3LYP, and M062X)')

    parser.add_argument('--no_D3', dest='D3_option', default="D3BJ ", action='store_const', const="",
                        help = 'When this flag is present, the D3-dispersion correction is disabled. By default the dispersion correction is used.')

    parser.add_argument('--remove_cross', dest='remove_cross', default=False, action='store_const', const=True,
                        help = 'When this flag is present, the cross-pairs are not parsed for configurations.')

    # Assign stdin arguments to variables
    args=parser.parse_args(argv)
    args.traj_files  = args.traj_files.split()
    args.map_files   = args.map_files.split()
    args.FF_db       = args.FF_db.split()
    args.frame_start = int(args.frame_start)
    args.frame_end   = int(args.frame_end)
    args.every       = int(args.every)
    args.filename    = str(args.filename)
    args.N_configs   = int(args.N_configs)
    args.q_A         = int(float(args.q_A))
    args.m_A         = int(float(args.m_A))
    args.q_B         = int(float(args.q_B))
    args.m_B         = int(float(args.m_B))
    args.procs       = int(args.procs)
    args.r_min_scale = float(args.r_min_scale)
    args.QC_types    = [ str(i).lower() for i in args.QC_types.split() ]
    args.o_mode      = int(args.o_mode)
    if len(args.QC_types) == 0: print("ERROR: -QC_types requires at least one argument ('dft' 'mp2' or 'dft mp2'). Exiting..."); quit()
    if len([ i for i in args.QC_types if i not in ['dft','mp2']]) > 0: print("ERROR: -QC_types only accepts 'mp2' or 'dft' as arguments. Exiting..."); quit()
    if args.A_list != []:
        args.A_list = parse_mol_list(args.A_list)
    if args.B_list != []:
        args.B_list = parse_mol_list(args.B_list)

    # If a directory is supplied then try and gather information on which cycle to start the simulation at. 
    if os.path.isdir(args.traj_files[0]):

        # Collect the md-* subfolders that are present in the supplied directory
        md_folders = [d for dp, dn, filenames in os.walk(args.traj_files[0]) for d in dn if fnmatch.fnmatch(d,"md-*") ]                    

        # If no md-* subfolders are present, print an error and exit
        if len(md_folders) == 0:
            print("ERROR: No viable molecular dynamics jobs were discovered in this folder. The script expects the folder to possess md-* labled subfolders. Exiting...")
            quit()

        # If a configs subfolder is present, read its subfolders to determine which cycles have already been sampled
        if os.path.isdir(args.traj_files[0]+"/configs") == True:
            try: parsed_cycles = sorted(set([ int(i[0].split('/')[-1].split('-')[0]) for i in os.walk(args.traj_files[0]+'/configs') if '-' in i[0].split('/')[-1] ]))                
            except: 
                print("ERROR: Cycle protocol failed, the script could not determine the cycle number from the contents of {}/configs. Exiting...".format(args.traj_files[0]))
                quit()

        # Else, start sampling at cycle 1 and create configs subfolder
        else:
            parsed_cycles = []
            if args.print_mol == 0:
                os.makedirs(args.traj_files[0]+"/configs")

        # Initialize output folder name, list of lammps trajectories and map files for drawing samples from, and a list of the configuration folder labels
        args.filename   = args.traj_files[0]
        md_folders      = [ i for i in md_folders if int(i.split('-')[-1]) not in parsed_cycles ]
        main_folder     = args.traj_files[0]
        args.traj_files = [ main_folder+'/{}/equil.lammpstrj'.format(i) for i in md_folders ]
        args.map_files  = [ main_folder+'/{}/{}.map'.format(i,i) for i in md_folders ]
        label_list      = [ int(i.split('-')[-1]) for i in md_folders ]

        # Check that the input trajectory exists
        for i in args.traj_files:
            if os.path.isfile(i) == False:
                print("ERROR: Couldn't find {}. Double-check the name and location of the trajectory file. Exiting...".format(i))
                quit()

        # Check that the map files exist
        for i in args.map_files:
            if os.path.isfile(i) == False:
                print("ERROR: Couldn't find {}. Double-check the name and location of the map file. Exiting...".format(i))
                quit()

        # Quit if there were no trajectories discovered that need to be parsed.
        if len(args.traj_files) == 0:
            print("\nNo trajectories were found that are in need of parsing. Exiting...")
            quit()

        # Start logger if args.print_mol == 0
        if args.print_mol == 0:
            sys.stdout = Logger(args.filename)
            print("PROGRAM CALL: python gen_jobs_for_vdw.py {}\n".format(' '.join([ i for i in argv])))

    # Else, use non-cycle direct protocol, a new folder is created and the configurations are sampled in a one-off fashion
    else:

        # Check that the input trajectory exists
        for i in args.traj_files:
            if os.path.isfile(i) == False:
                print("ERROR: Couldn't find {}. Double-check the name and location of the trajectory file. Exiting...".format(i))
                quit()

        # If a map file was not supplied then automatically generate a list of files basd on the root names of the trajectory folder
        if len(args.map_files) == 0:
            args.map_files += [ "/".join(i.split('/')[:-1])+'/'+i.split('/')[-2]+'.map' for i in args.traj_files ]

        # Else, pad the map files list with the last map file until there is a map file for each trajectory file
        else:
            while len(args.map_files) < len(args.traj_files):
                args.map_files += [args.map_files[-1]]                                  

        # Check that the map files exist
        for i in args.map_files:
            if os.path.isfile(i) == False:
                print("ERROR: Couldn't find {}. Double-check the name and location of the map file. Exiting...".format(i))
                quit()

        # Check to ensure that data won't be overwritten
        if os.path.exists(args.filename):
            print('ERROR: A folder named {} already exists in this directory. Exiting to avoid overwriting data...'.format(args.filename))
            quit()
        elif args.print_mol == 0:
            os.makedirs(args.filename)
            os.makedirs(args.filename+'/configs/')     
            sys.stdout = Logger(args.filename)
            print("PROGRAM CALL: python gen_jobs_for_vdw.py {}\n".format(' '.join([ i for i in argv])))

        # Initialize label list
        label_list = [ count_i+1 for count_i,i in enumerate(args.traj_files) ]

    # Read in parameters from database
    if len(args.FF_db) > 0:
        for i in args.FF_db:
            if os.path.isfile(i) == False:
                print("ERROR: Couldn't find {}. Double-check the name and location of the supplied database file. Exiting...".format(i))
                quit()
        FF_Data = get_FF_data(args.FF_db)

    # If multiple trajectories are being processed, check that all trajectories have matching molecule data
    if len(args.map_files) > 0:
        base_molecules, mol_id2type = parse_molecules(args.map_files[0])
        for count_t,t in enumerate(args.map_files):            
            tmp_molecules,tmp_mol_id2type = parse_molecules(t)
            if "{}".format(tmp_molecules) != "{}".format(base_molecules):
                print("ERROR: When parsing multiple trajectories, the script requires the same molecular composition in each trajectory. Exiting...")
                quit()

    # Parse the unique molecules from the simulation data
    for m in args.map_files:        
        molecules,mol_id2type = parse_molecules(m,verbose=1)

    # If the user just wanted to see the molecule info then the program is exited here.
    if args.print_mol == 1:
        quit()

    # Auto-populate A_list 
    if args.A_list == []:
        args.A_list = [ str(i) for i in list(mol_id2type.keys()) ]
    else:
        tmp = []
        for i in args.A_list:            
            if i > len(molecules):
                print("ERROR: mol number {} was in A_list, but not found in the list of unique molecules. Exiting...".format(i))
                quit()
            else:
                tmp += [ str(j) for j in list(mol_id2type.keys()) if mol_id2type[j] == i ]
        args.A_list = tmp

    # Auto-populate B_list 
    if args.B_list == []:
        args.B_list = [ str(i) for i in list(mol_id2type.keys()) ]
    else:
        tmp = []
        for i in args.B_list:
            if i > len(molecules):
                print("ERROR: mol number {} was in B_list, but not found in the list of unique molecules. Exiting...".format(i))
                quit()
            else:
                tmp += [ str(j) for j in list(mol_id2type.keys()) if mol_id2type[j] == i ]
        args.B_list = tmp

    # Loop over the trajectories that need to be sampled
    for count_t,t in enumerate(args.traj_files):

        # Check that the files are properly formatted
        status,types,mols,num_atoms = check_dump_format(t)

        if status==1:
            problems = [ t for count_i,i in enumerate(status) if i == 1 ]
            print("{}".format("*"*80))
            print("* {:76s} *".format("ERROR: the parser requires *sorted* *atom* dump types with the mol id in the 8th"))
            print("* {:76s} *".format("       column. Scripts exist to sort the dump files during post-processing."))
            print("{}".format("*"*80))        
            print("\nThe following file(s) had formatting errors:")
            for i in problems:
                print("\t{}".format(i))
            print("")
            quit()

        # Print diagnostic
        print("\n{}".format("*"*104))
        print("* {:^100s} *".format("Sampling trajectory {}".format(t)))
        print("{}\n".format("*"*104))

        # Look up the mapping between the datafile and the mapfile
        map2data_dict,data2map_dict,data2elem_dict,data2charge_dict = Generate_Map(t,args.map_files[count_t])

        # Generate adjacency matrix
        adj_mat =  Read_Adj_Mat(args.map_files[count_t])

        # Check for a *.in.settings file to read in VDW parameters for the md job. These parameters will be used as the initial guess for generating 
        # r_thresh and the initial guess in the eventual fitting procedure (these parameters are written to the output folder as "initial_guess.db" 
        initial_params = parse_pair_params(args.map_files[count_t],data2map_dict)

        # Generate pair_threshold dictionary. This holds the r_min thresholds for keeping each pair
        # that is used for sampling configurations and the configuration counts for each pair. These values are held in 
        # subdictionaries ("r_thresh" and "configurations"), keyed to the pair type (an ordered tuple of
        # the two atomtypes in the pair). For example, pair_dict[([1[6[1][1][1]]],[6[1][1][1][1]])]["configurations"] holds the
        # number of configurations satisfying r_thresh that have been sampled and pair_dict[([1[6[1][1][1]]],[6[1][1][1][1]])]["r_thresh"]
        # holds the r_threshold for keeping a configuration involving that pair. 
        pair_dict,initial_guess_dict = Generate_Pair_Info(t,args.A_list,args.B_list,data2map_dict,initial_params,r_min_scale=args.r_min_scale)

        # XXX remove all cross-pairs if indicated. 7/25/20
        if args.remove_cross == True:
            pair_dict_new = {x:pair_dict[x] for x in pair_dict.keys() if x[0]==x[1]}
            pair_dict = pair_dict_new
        
        # If a FF_dict file is supplied check for pairs that already have parameters
        if args.FF_db:
            del_keys = [ i for i in list(pair_dict.keys()) if i in list(FF_Data['vdws'].keys()) ]
            if len(del_keys) > 0:
                print("\n{}".format("*"*104))
                print("* {:^100s} *".format("The Following Pair Type(s) Are Already in the FF and will not be Sampled"))
                print("*{}*".format("-"*102))
                print("* {:<100s} *".format("Pair_type"))
                print("{}".format("*"*104))
                for i in del_keys:
                    print("  {:<80s}".format(str(i)))
                    pair_dict.pop(i,None)        
                    initial_guess_dict.pop(i,None)

        # Print summary of pair types being sampled for configurations
        print("\n{}".format("*"*104))
        print("* {:^100s} *".format("Summary of Pair Type(s) Being Sampled and Parsing Criteria"))
        print("*{}*".format("-"*102))
        print("* {:<70s}  {:<19s} {:<8s} *".format("Pair_type","r_thresh (ang)","Source"))
        print("{}".format("*"*104))
        for i in pair_dict:
            if i in list(initial_params.keys()):
                print("  {:<70s} {:< 20.4f} {:<9s}".format(str(i),pair_dict[i]["r_thresh"],str(initial_params[i][0])))        
            else:
                print("  {:<70s} {:< 20.4f} {:<9s}".format(str(i),pair_dict[i]["r_thresh"],"UFF"))        

        # Write initial parameter dictionary                
        with open(args.filename+'/initial_params.db','w') as f:
            f.write("\n# VDW definitions for initial guess\n#\n{:10s} {:60s} {:60s} {:11s} {:20s} {:20s}\n".format("#","Atom_type","Atom_type","Potential","eps (kcal/mol)","sigma (ang)"))        
            for i in pair_dict:
                f.write("{:10s} {:60s} {:60s} {:10s} {}\n".format("vdw",str(i[0]),str(i[1]),initial_guess_dict[i][0]," ".join([ "{:< 20.6f}".format(j) for j in initial_guess_dict[i][1:]])))        

        #############################################################
        # Parse the pair-wise configurations from LAMMPS trajectory #
        #############################################################
        print("\n{}".format("*"*104))
        print("* {:^100s} *".format("Parsing configurations from LAMMPS trajectory "))
        print("{}\n".format("*"*104))

        # Get the number of frames in each trajectory
        num_frames = get_num_frames(t)

        # Appropriately set the termination frame for each trajectory
        if args.frame_end == 0 or num_frames-1 < args.frame_end: # option for parsing all frames in all trajectories, or if less frames than the requested are present parse what's there. 
                args.frame_end = num_frames-1                    # set last frame to the total number of frames (-1 for zero indexing)

        # Parse configurations
        config_count = 0
        complete_flag = 0
        print("{:45} Configurations_Drawn".format('Frame_diagnostic'))
        print(range(args.frame_start,args.frame_end+1,args.every))
        for frame in range(args.frame_start,args.frame_end+1,args.every):

            # Print diagnostic
            #print("parsing configurations from frame {:<10}: {}".format(frame,config_count))
            print("start parsing configurations from frame {:<10}".format(frame))

            # Grab frame data (num_atoms allows more efficient initialization of the parser arrays)
            ids,types,geo,mol,box = get_frame(frame,t,num_atoms)

            # Unwrap molecules
            geo = unwrap(geo,adj_mat,mol,box)

            # Write unwrapped coords
            if args.write_unwrapped == 1:
                with open(args.filename+'/{}-unwrapped_frames.xyz'.format(label_list[count_t]),'a') as f:
                    f.write("{}\n\n".format(len(geo)))
                    for count_i,i in enumerate(geo):
                        f.write(" {:20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n".format(data2elem_dict[types[count_i]],i[0],i[1],i[2]))

            # Grab relevant pairs and 
            red_list = []
            tmp_count = 0
            for count_i,i in enumerate(args.A_list):
                for count_j,j in enumerate(args.B_list):

                    # avoid the same molecule interacting with itself, or duplicate pairs
                    if i == j or (i,j) in red_list:
                        continue
                    if False not in [ pair_dict[k]["configurations"] >= args.N_configs for k in list(pair_dict.keys()) ]:
                    # We stop when enough "atom type" pairs configurations are extracted which means mol pairs is gonna be less or more than the threshold configurations we set for atom pairs
                        break

                    # add the current pair to the red_list to avoid future duplicate configurations
                    red_list += [(i,j),(j,i)]

                    # Intialize subgeometries for each molecule
                    A_inds     = [ count_k for count_k,k in enumerate(mol) if k == i ] # atoms of mol i in A_list
                    B_inds     = [ count_k for count_k,k in enumerate(mol) if k == j ] # atoms of mol j in B_list
                    A_geo      = geo[A_inds,:]
                    B_geo      = geo[B_inds,:]
                    seps       = cdist(A_geo,B_geo) # all the atom pairs distance for mol pair (i,j)

                    # Check if any of the pairs are close enough to include the configuration
                    # A molecular pair was retained if at least one pair of atoms
                    # from each molecule resided within r_thresh
                    keep_flag = 0
                    for count_k,k in enumerate(seps):
                        for count_m,m in enumerate(k):

                            # Grab pair type
                            A_type = data2map_dict[types[A_inds[count_k]]]
                            B_type = data2map_dict[types[B_inds[count_m]]]

                            if A_type >= B_type:
                                pair = (A_type,B_type)
                            else:
                                pair = (B_type,A_type)

                            if pair in list(pair_dict.keys()) and m <= pair_dict[pair]["r_thresh"] and pair_dict[pair]["configurations"] < args.N_configs:
                                keep_flag = 1
                                pair_dict[pair]["configurations"] += 1  

                    # if the keep flag is set, generate orca input files and *.xyz files
                    if keep_flag == 1:

                        # Generate element lists
                        A_elem = [ data2elem_dict[types[m]] for m in A_inds ]
                        B_elem = [ data2elem_dict[types[m]] for m in B_inds ]

                        # Generate type lists
                        A_types = [ data2map_dict[types[m]] for m in A_inds ]
                        B_types = [ data2map_dict[types[m]] for m in B_inds ]
                        
                        
                        # Generate charge lists
                        A_charges = [ data2charge_dict[types[m]] for m in A_inds ] 
                        B_charges = [ data2charge_dict[types[m]] for m in B_inds ] 

                        # Generate xyz with geometry, types, mol_ids, and charges
                        gen_xyz(args.filename+'/configs/{}-{}/'.format(str(label_list[count_t]),str(config_count)),'{}-{}'.format(str(label_list[count_t]),str(config_count)),A_elem,B_elem,A_geo,B_geo,A_charges,B_charges,A_types,B_types)

                        if args.o_mode == 1:
                            # Generate SAPT input file
                            gen_SAPT(args.filename+'/configs/{}-{}/'.format(str(label_list[count_t]),str(config_count)),'{}-{}'.format(str(label_list[count_t]),str(config_count)),A_elem,A_types,A_geo,B_elem,B_types,B_geo,args.q_A,\
                                 args.m_A,args.q_B,args.m_B,args.SAPT_method,args.SAPT_name)
                            
                        else:
                            # Generate orca input file
                            gen_orca(args.filename+'/configs/{}-{}/'.format(str(label_list[count_t]),str(config_count)),'{}-{}'.format(str(label_list[count_t]),str(config_count)),A_elem,A_types,A_geo,B_elem,B_types,B_geo,args.functional,args.D3_option,args.q_A,\
                                 args.m_A,args.q_B,args.m_B,args.procs,args.QC_types)

                        # Iterate configuration counter
                        config_count += 1
            # Print diagnostic
            print("parsing configurations from frame {:<10}: {}".format(frame,config_count))

            # Break if enough configurations have been sampled
            if False not in [ pair_dict[k]["configurations"] >= args.N_configs for k in list(pair_dict.keys()) ]:
                complete_flag = 1
                break

        # Print diagnostic
        print("\nTotal configurations parsed: {}".format(config_count))

        # Print information about the missing pairs
        #fixed_pairs = []
        if False in [ pair_dict[k]["configurations"] >= args.N_configs for k in list(pair_dict.keys()) ]:
            print("\n{}".format("*"*104))
            print("* {:^100s} *".format("The parser Ran out of trajectory data before finding enough useful"))
            print("* {:^100s} *".format("configurations for all pair types"))
            print("*"*104)
            #print("\nThe following pair types will be added to fixed_params.db and held constant during fitting:\n")
            #for k in sorted(pair_dict.keys()):
            #    if pair_dict[k]["configurations"] < args.N_configs:
            #        fixed_pairs += [k]
            #        print("\t{}".format(k))

        # Write initial parameter dictionary
        #with open(args.filename+'/fixed_params.db','w') as f:
        #    f.write("\n# VDW definitions for inadequately sampled pairs\n#\n{:10s} {:60s} {:60s} {:11s} {:20s} {:20s}\n".format("#","Atom_type","Atom_type","Potential","eps (kcal/mol)","sigma (ang)"))        
        #    for i in fixed_pairs:
        #        f.write("{:10s} {:60s} {:60s} {:10s} {}\n".format("vdw",str(i[0]),str(i[1]),initial_guess_dict[i][0]," ".join([ "{:< 20.6f}".format(j) for j in initial_guess_dict[i][1:]])))        

    # Print completion banner
    print("\n{}".format("*"*104))
    print("* {:^100s} *".format("Parse Complete! The configurations are waiting in {}".format(args.filename)))
    print("{}\n".format("*"*104))

    return

# # Unwrap the pbc geometry of a frame
# def unwrap(coords,adj_mat,mol,box):

#     # Calculate box lengths
#     L = array([box[1]-box[0],box[3]-box[2],box[5]-box[4]])

#     # Cycle through atoms and unwrap long bonds. Only the count_i atom is moved in any step
#     # and once an atom has been moved based on a bond, both bonding atoms are left immobile. 
#     # The cycle is continued until the geometry is looped over without any translations being
#     # executed.
#     while 1:
#         placed = []
#         for count_i,i in enumerate(coords):
#             if count_i in placed: continue
#             bond_ind = [ count_j for count_j,j in enumerate(adj_mat[count_i]) if j == 1 ]
#             flag = 0 
#             for count_j,j in enumerate(bond_ind):

#                 # Calculate bond displacment
#                 delta = i-coords[j]

#                 # Apply periodic boundary conditions
#                 if delta[0] >  L[0]/2.0: coords[count_i,0] -= L[0]; placed += [count_i,j]; flag = 1
#                 if delta[0] < -L[0]/2.0: coords[count_i,0] += L[0]; placed += [count_i,j]; flag = 1
#                 if delta[1] >  L[1]/2.0: coords[count_i,1] -= L[1]; placed += [count_i,j]; flag = 1
#                 if delta[1] < -L[1]/2.0: coords[count_i,1] += L[1]; placed += [count_i,j]; flag = 1
#                 if delta[2] >  L[2]/2.0: coords[count_i,2] -= L[2]; placed += [count_i,j]; flag = 1
#                 if delta[2] < -L[2]/2.0: coords[count_i,2] += L[2]; placed += [count_i,j]; flag = 1

#                 # If the count_i atom has been moved, then break to avoid multiple translations
#                 if flag == 1:
#                     break                           

#         if placed == []:
#             break

#     # Find centroid of each molecule:    
#     centroids = {}
#     for i in sorted(set(mol)):
#         ind = [ count_j for count_j,j in enumerate(mol) if j == i ]
#         centroids[i] = mean(coords[ind,:],axis=0)

#     # Cycle through molecules and translate so that the centroid is within the simulation box
#     for i in set(mol):
#         ind = [ count_j for count_j,j in enumerate(mol) if j == i ]
#         if centroids[i][0] > box[1]: coords[ind,0] = coords[ind,0] - ceil((centroids[i][0]-box[1])/L[0])*L[0]
#         if centroids[i][0] < box[0]: coords[ind,0] = coords[ind,0] + ceil((box[0]-centroids[i][0])/L[0])*L[0]
#         if centroids[i][1] > box[3]: coords[ind,1] = coords[ind,1] - ceil((centroids[i][1]-box[3])/L[1])*L[1]
#         if centroids[i][1] < box[2]: coords[ind,1] = coords[ind,1] + ceil((box[2]-centroids[i][1])/L[1])*L[1]
#         if centroids[i][2] > box[5]: coords[ind,2] = coords[ind,2] - ceil((centroids[i][2]-box[5])/L[2])*L[2]
#         if centroids[i][2] < box[4]: coords[ind,2] = coords[ind,2] + ceil((box[4]-centroids[i][2])/L[2])*L[2]

#     return coords

# Unwrap the pbc geometry of a frame (UPDATED function)
def unwrap(coords,adj_mat,mol,box):

    # Calculate box lengths
    L = array([box[1]-box[0],box[3]-box[2],box[5]-box[4]])

    # Cycle through atoms and unwrap long bonds. Only the count_i atom is moved in any step
    # and once an atom has been moved based on a bond, both bonding atoms are left immobile. 
    # The cycle is continued until the geometry is looped over without any translations being
    # executed.

    # Initialize the set of placed atoms and the list of queued atoms
    placed = set([])

    # Loop over all atoms to ensure everything gets placed
    for count_i,i in enumerate(coords):
        
        # Skip atoms that have already been placed and if the queue is empty seed it with the current value
        if count_i in placed: continue
        queue = [count_i]

        # Loop over queued atoms
        for j in queue:                

            # Find bonds for the current atom (j) and add new bonded atoms to the end of the queue
            bond_ind = [ count_k for count_k,k in enumerate(adj_mat[j]) if k == 1 and count_k not in placed ]
            queue += [ k for k in bond_ind if k not in queue ]
            placed.update(bond_ind)

            # Loop over the bonded atoms and check if any need displacements
            for k in bond_ind:

                # Calculate bond displacment
                delta = coords[k]-coords[j]

                # Apply periodic boundary conditions
                if delta[0] >  L[0]/2.0: coords[k,0] -= L[0]
                if delta[0] < -L[0]/2.0: coords[k,0] += L[0]
                if delta[1] >  L[1]/2.0: coords[k,1] -= L[1]
                if delta[1] < -L[1]/2.0: coords[k,1] += L[1]
                if delta[2] >  L[2]/2.0: coords[k,2] -= L[2]
                if delta[2] < -L[2]/2.0: coords[k,2] += L[2]

    # Find centroid of each molecule:    
    centroids = {}
    for i in sorted(set(mol)):
        ind = [ count_j for count_j,j in enumerate(mol) if j == i ]
        centroids[i] = mean(coords[ind,:],axis=0)

    # Cycle through molecules and translate so that the centroid is within the simulation box
    for i in set(mol):
        ind = [ count_j for count_j,j in enumerate(mol) if j == i ]
        if centroids[i][0] > box[1]: coords[ind,0] = coords[ind,0] - ceil((centroids[i][0]-box[1])/L[0])*L[0]
        if centroids[i][0] < box[0]: coords[ind,0] = coords[ind,0] + ceil((box[0]-centroids[i][0])/L[0])*L[0]
        if centroids[i][1] > box[3]: coords[ind,1] = coords[ind,1] - ceil((centroids[i][1]-box[3])/L[1])*L[1]
        if centroids[i][1] < box[2]: coords[ind,1] = coords[ind,1] + ceil((box[2]-centroids[i][1])/L[1])*L[1]
        if centroids[i][2] > box[5]: coords[ind,2] = coords[ind,2] - ceil((centroids[i][2]-box[5])/L[2])*L[2]
        if centroids[i][2] < box[4]: coords[ind,2] = coords[ind,2] + ceil((box[4]-centroids[i][2])/L[2])*L[2]

    return coords

def gen_orca(Save_Folder,Name,A_Elements,A_Types,A_Geometry,B_Elements,B_Types,B_Geometry,functional,D3_option,charge_A=0,multiplicity_A=1,charge_B=1,multiplicity_B=1,procs=8,QC_types=['dft','mp2']):

    # Print diagnostic message
    Save_Name = Save_Folder+str(Name)+'.in'

    # Check for Li Be He and H, so that the frozen core approximation can be avoided
    No_Frozen_List = ["Li","Be","He","H"]    
    No_Frozen_Flag = 0
    for i in No_Frozen_List:
        if ( len(A_Elements) == 1 and i in A_Elements ) or ( len(B_Elements) == 1 and i in B_Elements ):
            No_Frozen_Flag = 1

    # Make subdirectory to hold the current vdw files
    if os.path.isdir(Save_Folder) == False:
        os.makedirs(Save_Folder)

    # Open input file
    with open(Save_Name,'w') as f:

        ############################################################
        # Run DFT and MP2 single-points of the dimer configuration #
        ############################################################

        # Write DFT Header
        if 'dft' in QC_types:

            # Avoid dispersion correction
            if No_Frozen_Flag == 1:
                if procs == 1:
                    if charge_A < 0 or charge_B < 0:
                        f.write('#Run DFT single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF CHELPG PMODEL\n\n%base "{}_DFT_AB"'.format(functional,Name))
                    else:
                        f.write('#Run DFT single point of the dimer\n! {} def2-TZVP TIGHTSCF CHELPG PMODEL\n\n%base "{}_DFT_AB"'.format(functional,Name))
                else:
                    if charge_A < 0 or charge_B < 0:
                        f.write('#Run DFT single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF CHELPG PMODEL PAL{}\n\n%base "{}_DFT_AB"'.format(functional,procs,Name))
                    else:
                        f.write('#Run DFT single point of the dimer\n! {} def2-TZVP TIGHTSCF CHELPG PMODEL PAL{}\n\n%base "{}_DFT_AB"'.format(functional,procs,Name))

            # Use dispersion correction
            # If --noD3 flag is on them force no dispersion correction (this is for wB97X-D3)
            else:
                if procs == 1:
                    if charge_A < 0 or charge_B < 0:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF {}CHELPG PMODEL\n\n%base "{}_DFT_AB"'.format(functional,D3_option,Name))
                    else:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} def2-TZVP TIGHTSCF {}CHELPG PMODEL\n\n%base "{}_DFT_AB"'.format(functional,D3_option,Name))
                else:
                    if charge_A < 0 or charge_B < 0:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF {}CHELPG PMODEL PAL{}\n\n%base "{}_DFT_AB"'.format(functional,D3_option,procs,Name))
                    else:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} def2-TZVP TIGHTSCF {}CHELPG PMODEL PAL{}\n\n%base "{}_DFT_AB"'.format(functional,D3_option,procs,Name))

            # Increase default max scf iterations
            f.write('\n\n%scf\nMaxIter 1000\nend\n')

            # Write Coordinate header
            f.write('\n* xyz {} {}\n'.format(charge_A+charge_B,max(multiplicity_A,multiplicity_B)))

            # Write A_Geometry 
            for count_j,j in enumerate(A_Geometry):
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))

            # Write B_Geometry 
            for count_j,j in enumerate(B_Geometry):
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

        # Write RI-MP2 single point calculation on the dimer
        if 'mp2' in QC_types:
            if procs == 1:
                if No_Frozen_Flag == 1:
                    if charge_A < 0 or charge_B < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(Name))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(Name))
                else:
                    if charge_A < 0 or charge_B < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG\n'.format(Name))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG\n'.format(Name))
            else:
                if No_Frozen_Flag == 1:
                    if charge_A < 0 or charge_B < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore PAL{}\n'.format(Name,procs))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore PAL{}\n'.format(Name,procs))
                else:
                    if charge_A < 0 or charge_B < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG PAL{}\n'.format(Name,procs))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG PAL{}\n'.format(Name,procs))

            # Increase default max scf iterations
            f.write('\n%scf\nMaxIter 1000\nend\n')

            # Write coordinate header
            f.write('\n* xyz {} {}\n'.format(charge_A+charge_B,max(multiplicity_A,multiplicity_B)))

            # Write A_Geometry 
            for count_j,j in enumerate(A_Geometry):
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))

            # Write B_Geometry 
            for count_j,j in enumerate(B_Geometry):
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

        ############################################################################
        # Run DFT and MP2 single-points of the molecule A with the dimer basis set #
        ############################################################################

        # Write DFT Header
        if 'dft' in QC_types:

            # Avoid dispersion correction
            if No_Frozen_Flag == 1:

                f.write('\n#Run DFT single point of the molecule with dimer basis set\n')
                if charge_A < 0 or charge_B < 0:
                    f.write('$new_job\n%base "{}_DFT_A"\n! {} ma-def2-TZVP TIGHTSCF CHELPG PMODEL\n'.format(Name,functional))
                else:
                    f.write('$new_job\n%base "{}_DFT_A"\n! {} def2-TZVP TIGHTSCF CHELPG PMODEL\n'.format(Name,functional))

            # Use dispersion correction
            else:

                f.write('\n#Run DFT-D3 single point of the molecule with dimer basis set\n')
                if charge_A < 0 or charge_B < 0:
                    f.write('$new_job\n%base "{}_DFT_A"\n! {} ma-def2-TZVP TIGHTSCF {}CHELPG PMODEL\n'.format(Name,functional,D3_option))
                else:
                    f.write('$new_job\n%base "{}_DFT_A"\n! {} def2-TZVP TIGHTSCF {}CHELPG PMODEL\n'.format(Name,functional,D3_option))

            # Increase default max scf iterations
            f.write('\n%scf\nMaxIter 1000\nend\n')

            # Write molecule A geometry, and molecule B geometry with ghost basis functions (: notation)
            f.write('\n* xyz {} {}\n'.format(charge_A,multiplicity_A))
            for count_j,j in enumerate(A_Geometry):                    
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))
            for count_j,j in enumerate(B_Geometry):                    
                f.write('  {:<20s} : {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

        # Write RI-MP2 single point calculation on the dimer
        if 'mp2' in QC_types:
            if No_Frozen_Flag == 1:
                if charge_A < 0 or charge_B < 0:
                    f.write('\n#Run RI-MP2 single point of the molecule with dimer basis set\n$new_job\n%base "{}_RI-MP2_A"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(Name))
                else:
                    f.write('\n#Run RI-MP2 single point of the molecule with dimer basis set\n$new_job\n%base "{}_RI-MP2_A"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(Name))
            else:
                if charge_A < 0 or charge_B < 0:
                    f.write('\n#Run RI-MP2 single point of the molecule with dimer basis set\n$new_job\n%base "{}_RI-MP2_A"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG\n'.format(Name))
                else:
                    f.write('\n#Run RI-MP2 single point of the molecule with dimer basis set\n$new_job\n%base "{}_RI-MP2_A"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG\n'.format(Name))

            # Increase default max scf iterations
            f.write('\n%scf\nMaxIter 1000\nend\n')

            # Write molecule A geometry, and molecule B geometry with ghost basis functions (: notation)
            f.write('\n* xyz {} {}\n'.format(charge_A,multiplicity_A))
            for count_j,j in enumerate(A_Geometry):                    
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))
            for count_j,j in enumerate(B_Geometry):                    
                f.write('  {:<20s} : {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

        ##########################################################################
        # Run DFT and MP2 single-points of the particle with the dimer basis set #
        ##########################################################################

        # Write DFT Header
        if 'dft' in QC_types:

            # Avoid dispersion correction
            if No_Frozen_Flag == 1:

                f.write('\n#Run DFT single point of the particle with dimer basis set\n')
                if charge_A < 0 or charge_B < 0:
                    f.write('$new_job\n%base "{}_DFT_B"\n! {} ma-def2-TZVP TIGHTSCF CHELPG PMODEL\n'.format(Name,functional))
                else:
                    f.write('$new_job\n%base "{}_DFT_B"\n! {} def2-TZVP TIGHTSCF CHELPG PMODEL\n'.format(Name,functional))

            # Use dispersion correction
            else:

                f.write('\n#Run DFT-D3 single point of the particle with dimer basis set\n')
                if charge_A < 0 or charge_B < 0:
                    f.write('$new_job\n%base "{}_DFT_B"\n! {} ma-def2-TZVP TIGHTSCF {}CHELPG PMODEL\n'.format(Name,functional,D3_option))
                else:
                    f.write('$new_job\n%base "{}_DFT_B"\n! {} def2-TZVP TIGHTSCF {}CHELPG PMODEL\n'.format(Name,functional,D3_option))

            # Increase default max scf iterations
            f.write('\n%scf\nMaxIter 1000\nend\n')

            # Write molecule A geometry with ghost basis functions (: notation), and molecule B geometry as the active molecule.
            f.write('\n* xyz {} {}\n'.format(charge_B,multiplicity_B))
            for count_j,j in enumerate(A_Geometry):                    
                    f.write('  {:<20s} : {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))
            for count_j,j in enumerate(B_Geometry):                    
                    f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

        # Write RI-MP2 single point calculation on the dimer
        if 'mp2' in QC_types:
            if No_Frozen_Flag == 1:
                if charge_A < 0 or charge_B < 0:
                    f.write('\n#Run RI-MP2 single point of the particle with dimer basis set\n$new_job\n%base "{}_RI-MP2_B"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(Name))
                else:
                    f.write('\n#Run RI-MP2 single point of the particle with dimer basis set\n$new_job\n%base "{}_RI-MP2_B"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(Name))
            else:
                if charge_A < 0 or charge_B < 0:
                    f.write('\n#Run RI-MP2 single point of the particle with dimer basis set\n$new_job\n%base "{}_RI-MP2_B"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG\n'.format(Name))
                else:
                    f.write('\n#Run RI-MP2 single point of the particle with dimer basis set\n$new_job\n%base "{}_RI-MP2_B"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG\n'.format(Name))

            # Increase default max scf iterations
            f.write('\n%scf\nMaxIter 1000\nend\n')

            # Write molecule A geometry with ghost basis functions (: notation), and molecule B geometry as the active molecule.
            f.write('\n* xyz {} {}\n'.format(charge_B,multiplicity_B))
            for count_j,j in enumerate(A_Geometry):                    
                f.write('  {:<20s} : {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))
            for count_j,j in enumerate(B_Geometry):                    
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

    return

# Generates an input file to run a SAPT job in PSI4.
# Stephen Shiring, October 25, 2017
#SAPT_order = sapt0, ssapt0, sapt2
def gen_SAPT(Save_Folder,Name,A_Elements,A_Types,A_Geometry,B_Elements,B_Types,B_Geometry,charge_A=0,multiplicity_A=1,charge_B=1,multiplicity_B=1,SAPT_order='ssapt0',SAPT_complex_name='complex'):
    # Print diagnostic message
    Save_Name = Save_Folder+str(Name)+'.in'

    # Make subdirectory to hold the current vdw files
    if os.path.isdir(Save_Folder) == False:
        os.makedirs(Save_Folder)

    # Open input file
    with open(Save_Name,'w') as f:
        f.write('molecule ' + str(SAPT_complex_name) + ' {\n')
        f.write('     ' + str(charge_A) + ' ' + str(multiplicity_A) + '\n')
        
        for count_j,j in enumerate(A_Geometry):
            f.write('     {:5s} {:< 10.6f} {:< 10.6f} {:< 10.6f}\n'.format(A_Elements[count_j],j[0],j[1],j[2]))
        
        f.write('     --\n')
        f.write('     ' + str(charge_B) + ' ' + str(multiplicity_B) + '\n')
        
        for count_j,j in enumerate(B_Geometry):    
            f.write('     {:5s} {:< 10.6f} {:< 10.6f} {:< 10.6f}\n'.format(B_Elements[count_j],j[0],j[1],j[2]))
        
        f.write('     units angstrom\n}\n\n')
        
        if str(SAPT_order) == 'sapt0':
            f.write('set globals {\n    basis          jun-cc-pvdz\n    df_basis_scf   jun-cc-pvdz-jkfit\n    df_basis_sapt  jun-cc-pvdz-ri\n    guess          sad\n    scf_type       df\n}\n\n')
            f.write('set sapt {\n    print          1\n}\n\n')
            f.write('energy(\'sapt0\')')
        
        elif str(SAPT_order) == 'ssapt0':
            f.write('set globals {\n    basis          jun-cc-pvdz\n    df_basis_scf   jun-cc-pvdz-jkfit\n    df_basis_sapt  jun-cc-pvdz-ri\n    guess          sad\n    scf_type       df\n}\n\n')
            f.write('set sapt {\n    print          1\n}\n\n')
            f.write('energy(\'ssapt0\')')

        elif str(SAPT_order) == 'sapt2':
            f.write('set globals {\n    basis         aug-cc-pvdz\n    df_basis_scf  aug-cc-pvdz-jkfit\n    df_basis_sapt aug-cc-pvdz-ri\n    guess         sad\n    scf_type      df\n}\n\n')
            f.write('set sapt {\n    print         1\n    freeze_core   true\n}\n\n')
            f.write('energy(\'sapt2\')')
        
    return


# Parses the box,id#s,atom_types,mol#s, and geometry from the requested frame in the supplied trajectory
def get_frame(frame,traj,num_atoms):

    ids = [0]*num_atoms
    types = ["X"]*num_atoms
    geo = zeros([num_atoms,3])
    mol = ["X"]*num_atoms
    box = array([0.0,0.0,0.0,0.0,0.0,0.0])

    # Check for atom dump type
    with open(traj,'r') as f:
        timestep_flag = -1        # zero indexing
        parse_flag    = 0         # toggles the geometry parse
        box_flag      = 0
        atom_counter  = 0         # keeps track of the number of placed atoms
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) > 0 and fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                timestep_flag += 1
            if len(fields) > 0 and fields[0] == "ITEM:" and fields[1] == "TIMESTEP" and parse_flag == 1:
                break

            # Geo parse statements
            if timestep_flag == frame and len(fields) >= 8 and fields[0] == "ITEM:" and fields[1] == "ATOMS" and "mol" in fields:
                mol_ind = fields.index("mol") - 2
                parse_flag = 1
                continue
            if parse_flag == 1:
                ids[atom_counter]   = int(fields[0])
                types[atom_counter] = fields[1]
                geo[atom_counter,:] = array([float(fields[2]),float(fields[3]),float(fields[4])])
                mol[atom_counter]   = fields[mol_ind]
                atom_counter += 1

            # Box parse statements
            if timestep_flag == frame and len(fields) >= 2 and fields[0] == "ITEM:" and fields[1] == "BOX":
                box_flag = 1
                continue
            if box_flag > 0 and box_flag <= 5:
                box[box_flag-1] = float(fields[0])
                box[box_flag]   = float(fields[1])
                box_flag += 2


    return ids,types,geo,mol,box

# Check dump formatting and scrap lammps atom types and mol ids
def check_dump_format(Name):

    # Initialize a list of statuses for each dump file
    status    = 0
    num_atoms = 0
    types     = []
    mols      = []

    # Check for atom dump type
    with open(Name,'r') as f:
        timestep_flag = 0
        atoms_flag     = 0
        for lc,lines in enumerate(f):
            if lc >= 1000:
                status = 1
                break
            fields = lines.split()
            if len(fields) > 0 and fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                timestep_flag = 1
            if len(fields) >= 6 and fields[0] == "ITEM:" and fields[1] == "ATOMS" and "mol" in fields:
                atoms_flag = 1
                mol_ind = fields.index("mol") - 2
            if timestep_flag == 1 and atoms_flag == 1:
                break

    # Continue if the current dump file isn't formatted as "atom" type
    if status == 1:
        return status,types,num_atoms

    # Parse atom types and mol ids
    with open(Name,'r') as f:
        time_count=0
        write_flag=0
        atom_count=0
        for lines in f:
            fields = lines.split()

            # ID start of frame
            if fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                time_count += 1
                if time_count == 2:
                    break

            # ID start of the geometry
            if fields[0] == "ITEM:" and fields[1] == "ATOMS":
                write_flag = 1
                continue

            # Save atom types
            if write_flag == 1:                    
                atom_count += 1
                if fields[1] not in types:
                    types += [fields[1]]
                if fields[mol_ind] not in mols:
                    mols  += [fields[mol_ind]]

    return status,types,mols,atom_count

def Read_Adj_Mat(map_file):
    # Read number of atoms
    with open(map_file,'r') as f:
        for count_i,i in enumerate(f):
            fields = i.split()
            if count_i == 0:
                adj_mat = zeros([int(fields[0]),int(fields[0])])
            if count_i > 1:
                if len(fields) > 5:
                    ind = [ int(j) for j in fields[5:] ]
                    adj_mat[count_i-2,ind] = 1
                    adj_mat[ind,count_i-2] = 1
    return adj_mat

# This function returns 2 dictionaries that are used to map from the 
# database atomtypes to the datafile atomtypes. The first dictionary
# returned, map2data_dict, has the database atomtypes for keys that
# return the datafile atomtype as entries. The second dictionary returned,
# data2map_dict, has the datafile atomtypes for keys that return the
# database atomtypes as entries.
def Generate_Map(traj_file,map_file):

    # Parse the mapfile for the list of atom_types, elements, and charges
    with open(map_file,'r') as f:
        for count_i,i in enumerate(f):

            fields = i.split()

            if count_i == 0:
                types = ["X"]*int(fields[0])
                number_of_chains = int(fields[1])
                elements = ["X"]*int(fields[0])
                charges  = zeros(int(fields[0]))

            elif count_i > 1 and len(fields) >= 5: 
                types[count_i-2] = fields[0]
                elements[count_i-2] = fields[1]
                charges[count_i-2] = float(fields[4])

    # Find the indices in types for the first occurence 
    # of each type. When parsing the datafile only these
    # atoms are used for mapping between the .data and .map
    # atomtypes.
    ind = []
    type_set = []
    element_set = []
    charge_set = []
    for i in set(types):
        for count_j,j in enumerate(types):
            if j == i:
                ind += [count_j+1]
                type_set += [i]
                element_set += [elements[count_j]]
                charge_set += [charges[count_j]]
                break

    # Iterate over the trajectory file until all atomtypes have been mapped
    data2map_dict    = {} # keys are the atom types in the datafile, entries are the atom types in the mapfile
    map2data_dict    = {} # keys are the atom types in the mapfile, entries are the atom types in the datafile
    data2elem_dict   = {} # keys are the atom types in the datafile, entries are the elements mapfile
    data2charge_dict = {} # keys are the atom types in the datafile, entries are the charges mapfile
    num_assign=0    
    num_types=len(type_set)
    parse_flag=0
    time_count = 0
    with open(traj_file,'r') as f:
        for count_i,i in enumerate(f):
            fields = i.split()

            # ID start of frame
            if fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                time_count += 1
                if time_count == 2:
                    break

            # ID start of the geometry
            if fields[0] == "ITEM:" and fields[1] == "ATOMS":
                parse_flag = 1
                continue

            if parse_flag == 1:

                # Parse the atom info if its index is in ind
                if int(fields[0]) in ind:
                    for count_j,j in enumerate(ind):
                        if j == int(fields[0]):
                            map2data_dict[str(type_set[count_j])]=int(fields[1])
                            data2map_dict[fields[1]]=type_set[count_j]
                            data2elem_dict[fields[1]]=element_set[count_j]
                            data2charge_dict[fields[1]]=charge_set[count_j]

    return map2data_dict,data2map_dict,data2elem_dict,data2charge_dict

# Returns the total number of frames in the trajectory
def get_num_frames(Name):

    # Parse number of frames    
    with open(Name,'r') as f:
        frame_count = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) > 0 and fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                frame_count += 1

    return frame_count

def gen_xyz(Save_Folder,Name,A_Elements,B_Elements,A_Geometry,B_Geometry,A_Charges,B_Charges,A_Types,B_Types):

    # Print diagnostic message
    Save_Name = Save_Folder+str(Name)+'.in'

    # Make subdirectory to hold the current vdw files
    if os.path.isdir(Save_Folder) == False:
        os.makedirs(Save_Folder)

    # Save an xyz for viewing the configuration in VMD
    with open(Save_Folder+str(Name)+'.xyz','w') as f:
        f.write('{}\n\n'.format(len(A_Geometry)+len(B_Geometry)))
        for count_j,j in enumerate(A_Geometry):
            f.write('  {:<10s} {:< 20.6f} {:< 20.6f} {:< 20.6f} {:<12d} {:< 12.6f} {:<60s}\n'.format(A_Elements[count_j],j[0],j[1],j[2],0,A_Charges[count_j],A_Types[count_j]))
        for count_j,j in enumerate(B_Geometry):
            f.write('  {:<10s} {:< 20.6f} {:< 20.6f} {:< 20.6f} {:<12d} {:< 12.6f} {:<60s}\n'.format(B_Elements[count_j],j[0],j[1],j[2],1,B_Charges[count_j],B_Types[count_j]))

def Generate_Pair_Info(traj_file,A_list,B_list,data2map_dict,initial_params,r_min_scale=1.5):

    # Iterate over the trajectory file until all unique pairs between molecules in A_list and 
    # molecules in B_list have been parsed.
    pair_dict        = {} # keys are the atom types in the datafile, entries are the atom types in the mapfile
    types_A          = []
    types_B          = []
    parse_flag=0
    time_count = 0
    with open(traj_file,'r') as f:
        for count_i,i in enumerate(f):
            fields = i.split()

            # ID start of frame
            if fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                time_count += 1
                if time_count == 2:
                    break

            # ID start of the geometry
            if fields[0] == "ITEM:" and fields[1] == "ATOMS":
                parse_flag = 1
                continue

            if parse_flag == 1:

                # Parse the atom info if its index is in ind
                if fields[5] in A_list:
                    types_A += [data2map_dict[fields[1]]]
                if fields[5] in B_list:
                    types_B += [data2map_dict[fields[1]]]

    # Convert the types to unique lists
    types_A = sorted(set(types_A))
    types_B = sorted(set(types_B))

    # Generate a unique list of pairs
    pairs = []
    for i in types_A:
        for j in types_B:
            if i >= j and (i,j) not in pairs:
                pairs += [(i,j)]
            elif j > i and (j,i) not in pairs:
                pairs += [(j,i)]

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # Initialize pair_dict. The r_thresh and configuration count for each pair are held in subdictionaries keyed to the pair type
    pair_dict = {}
    initial_guess_dict = {}   # Used for writing the initial guess pair database file for the job. 
    for i in pairs:
        pair_dict[i] = {}

        # If the parameters were found in the *in.settings file then use those
        if i in list(initial_params.keys()):
            pair_dict[i]["r_thresh"] = initial_params[i][3]*2.**(1./6.) * r_min_scale
            initial_guess_dict[i] = initial_params[i][1:]

        # Default to UFF 
        else:
            type_1 = int(i[0].split('[')[1].split(']')[0])
            type_2 = int(i[1].split('[')[1].split(']')[0])
            r_min = ((UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0)*2.**(1./6.)    # convert sigma into r_min of the potential using the 2^(1/6) factor
            pair_dict[i]["r_thresh"]       = r_min*r_min_scale                     # Scale r_min for configurational sampling and key it to r_thresh
            initial_guess_dict[i] = [ 'lj', ( UFF_dict[type_1][0]*UFF_dict[type_2][0] )**(0.5), (UFF_dict[type_1][1] + UFF_dict[type_2][1] ) / 2.0 ]

        # Create an entry to count configurations
        pair_dict[i]["configurations"] = 0                                     

    return pair_dict,initial_guess_dict

# Description: Reads in the FF parameters
def get_FF_data(FF_files,keep_types=[ "atom", "vdw", "bond", "angle", "torsion", "dihedral", "charge" ]):

    Data = {"atoms":{},"bonds":{},"angles":{},"dihedrals":{},"vdws":{},"charges":{}}
    
    # Loop over the files and read in FF parameters
    for i in FF_files:
        with open(i,'r') as f:
            for lines in f:
                fields = lines.split()

                if len(fields) == 0: continue
                if fields[0] == "atom" and "atom" in keep_types:
                    Data["atoms"][fields[1]] = fields[1:]
                if fields[0] == "vdw" and "vdw" in keep_types:
                    Data["vdws"][(fields[1],fields[2])] = fields[1:]
                if fields[0] == "bond" and "bond" in keep_types:
                    Data["bonds"][(fields[1],fields[2],fields[3])] = fields[1:]
                if fields[0] == "angle" and "angle" in keep_types:
                    Data["angles"][(fields[1],fields[2],fields[3],fields[4])] = fields[1:]
                if fields[0] == "diehdral" or fields[0] == "torsion" and ( "diehdral" in keep_types or "torsion" in keep_types ):
                    Data["dihedrals"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = fields[1:]
                if fields[0] == "charge" and "charge" in keep_types:
                    Data["charges"][fields[1]] = fields[1:]

    return Data

# this function converts the input string into a list of molecule indices
# the bulk of the code is dedicated to the ":" mini-language that allows the user
# to expand the molecule indices using ranges (e.g., 1:6 7 10 means 1,2,3,4,5,6,7,10
# and 2:4:2 3:3:9 means 2,4,3,6,9)
def parse_mol_list(mol_list):

    entries = mol_list.split()
    for count_i,i in enumerate(entries):
        if ":" in i:
            entries[count_i] = [ int(j) for j in i.split(":") ]
        else:
            entries[count_i] = int(i)

    expanded_list = []
    for i in entries:
        if type(i) == int:
            expanded_list += [i]
        elif len(i) == 2:
            expanded_list += list(range(i[0],i[1]+1))
        elif len(i) == 3:
            expanded_list += list(range(i[0],i[2]+1,i[1]))

    return [ i for i in expanded_list ]

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/gen_jobs_for_vdw.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass




# Description: This function searches through the configurations being parsed for unique
#              molecules. All FF_parameters are linked to examples of the molecule(s) used
#              for their parametrization. 
# Note:        identifying unique molecules is generally a graph isomorphism problem
#              (P difficult). The algorithm used here assumes that replica molecules have
#              same atomic ordering (which the vdw_gen.py parser should guarrantee) which 
#              makes identifying the molecule types much easier.
#             
# Inputs:      folder: folder for saving the molecule.db file
#              Data:   the dictionary holding the configuration data 
def parse_molecules(mapfile,verbose=0):

    # Print banner
    if verbose != 0:
        print("\n{}".format("*"*104))
        print("* {:^100s} *".format("Parsing unique molecules from {}".format(mapfile)))
        print("{}".format("*"*104))

    # Build adjacency matrix from map file
    with open(mapfile,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            if lc == 0: 
                adj_mat = zeros([int(fields[0]),int(fields[0])])
                atomtypes = ["X"]*int(fields[0])
                charges = [0.0]*int(fields[0])
                continue
            if lc == 1: continue
            if len(fields) == 0: continue
            elif len(fields) < 6: atomtypes[lc-2]=fields[0]; charges[lc-2]=float(fields[4])
            else: 
                atomtypes[lc-2]=fields[0]
                charges[lc-2]=float(fields[4])
                ind = [ int(i) for i in fields[5:] ]
                adj_mat[lc-2,ind] = 1
                adj_mat[ind,lc-2] = 1

    # Find subgraphs
    subgraph_list = []
    for count_i,i in enumerate(adj_mat):

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

        # Add new subgraphs to the master list
        if subgraph not in subgraph_list:
            subgraph_list += [subgraph]

    # Create molecule tuples holding the list of atomtypes as the first element and the adj mat as the second element
    mol_id2type = {}
    adj_mat_list  = []
    atomtype_list = []
    charges_list  = []
    for count_i,i in enumerate(subgraph_list):
        mol_adj_mat = zeros([len(i),len(i)])
        for count_j,j in enumerate(i):
            mol_adj_mat[count_j,:] = adj_mat[j,i] # NOTE: i is a list whereas j is an int
            mol_adj_mat[:,count_j] = adj_mat[i,j] # NOTE: i is a list whereas j is an int
        mol_atomtypes = [ atomtypes[j] for j in i ]
        mol_charges   = [ charges[j] for j in i ]

        # Check if this is a new molecule based on comparisons against adj_mats, atom types, and charges. 
        match = 0
        for count_j,j in enumerate(adj_mat_list):
            if array_equal(mol_adj_mat,j) and atomtype_list[count_j] == mol_atomtypes and charges_list[count_j] == mol_charges:
                mol_id2type[count_i] = count_j
                match = 1
                break

        # If no match was found then add the new molecule info to the appropriate lists
        if match == 0:
            adj_mat_list += [mol_adj_mat]
            atomtype_list += [mol_atomtypes]
            charges_list += [mol_charges]
            mol_id2type[count_i] = len(adj_mat_list)-1

    # Reassign moltypes based on size (smallest first) then charge 
    # NOTE: scaling the size by 1000 opens up space to rank by charge as a secondary measure
    # NOTE: -sum(charges) leads to cations being ranked first.   
    scores = [ len(i)*1000.0-sum(charges_list[count_i]) for count_i,i in enumerate(atomtype_list) ]
    ind_map = sorted(list(range(len(scores))), key=lambda k: scores[k])
    for i in mol_id2type: mol_id2type[i] = ind_map.index(mol_id2type[i])
    scores,adj_mat_list,atomtype_list,charges_list = list(zip(*sorted(zip(scores,adj_mat_list,atomtype_list,charges_list))))

    molecules = []
    for count_i,i in enumerate(atomtype_list):
        molecules += [(i,charges_list[count_i],adj_mat_list[count_i])]

    # Print diagnostic and return
    if verbose != 0:
        print("Number of unique molecules in the fit configurations: {}".format(len(molecules)))

        for count_i,i in enumerate(molecules):
            print("\nMolecule: {}".format(count_i))
            print("\n\t{:<40s}  {:<40s}".format("Atom_types","Charges"))
            for count_j,j in enumerate(i[0]):
                print("\t{:<40s} {:< 12.6f}".format(j,i[1][count_j]))
            list_to_print = [ str(j) for j in list(mol_id2type.keys()) if mol_id2type[j]==count_i ]
            print("\n\t{:<40s} {}".format("Correspondence with trajectory mol_id:",', '.join(list_to_print[:10])))
            if len(list_to_print) > 10:
                for j in range(int(ceil((len(list_to_print[10:])/10.0)))):
                    print("\t{} {}".format(" "*40,', '.join(list_to_print[10*(j+1):10*(j+2)])))
            print("\tTotal charge: {:< d}".format(int(round(sum(i[1])))))

    return molecules,mol_id2type

# Parse pair parameters from *.in.settings file if it is present
# basic protocols are in place to detect hybrid types / non-lj pair types. 
def parse_pair_params(map_file,data2map_dict):

    settings_file = ".".join([ i for i in map_file.split(".")[:-1] ])+'.in.settings'
    pair_file = ".".join([ i for i in map_file.split(".")[:-1] ])+'.pairs'
    initial_params = {}

    # Parse parameters from the pair file if it is present
    if os.path.isfile(pair_file):
        check_flag = 1
        with open(pair_file,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "pair_coeff":

                    # Avoid redundant addition
                    if (data2map_dict[fields[1]],data2map_dict[fields[2]]) in list(initial_params.keys()):
                        continue

                    # Check for hybrid formatting
                    if check_flag == 1:
                        try:
                            float(fields[3])
                            hybrid_flag = 0
                        except ValueError:
                            hybrid_flag = 1
                        check_flag = 0 

                    # Parse parameters
                    if hybrid_flag == 0:                        
                        initial_params[ (data2map_dict[fields[1]],data2map_dict[fields[2]]) ] = [ 'pair_file','lj', float(fields[3]), float(fields[4]) ]
                        initial_params[ (data2map_dict[fields[2]],data2map_dict[fields[1]]) ] = [ 'pair_file','lj', float(fields[3]), float(fields[4]) ]
                    elif hybrid_flag == 1:
                        if 'lj' not in fields[3]:
                            print("WARNING: the FF parsing algorithms only have compatibility with lj type potentials.")
                        initial_params[ (data2map_dict[fields[1]],data2map_dict[fields[2]]) ] = [ 'pair_file','lj', float(fields[4]), float(fields[5]) ]
                        initial_params[ (data2map_dict[fields[2]],data2map_dict[fields[1]]) ] = [ 'pair_file','lj', float(fields[4]), float(fields[5]) ]

    # Parse parameters from the settings file if it is present
    if os.path.isfile(settings_file):
        check_flag = 1
        with open(settings_file,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "pair_coeff":

                    # Avoid redundant addition
                    if (data2map_dict[fields[1]],data2map_dict[fields[2]]) in list(initial_params.keys()):
                        continue

                    # Check for hybrid formatting
                    if check_flag == 1:
                        try:
                            float(fields[3])
                            hybrid_flag = 0
                        except ValueError:
                            hybrid_flag = 1
                        check_flag = 0 

                    # Parse parameters
                    if hybrid_flag == 0:
                        initial_params[ (data2map_dict[fields[1]],data2map_dict[fields[2]]) ] = [ 'settings','lj', float(fields[3]), float(fields[4]) ]
                        initial_params[ (data2map_dict[fields[2]],data2map_dict[fields[1]]) ] = [ 'settings','lj', float(fields[3]), float(fields[4]) ]
                    elif hybrid_flag == 1:
                        if 'lj' not in fields[3]:
                            print("WARNING: the FF parsing algorithms only have compatibility with lj type potentials.")
                        initial_params[ (data2map_dict[fields[1]],data2map_dict[fields[2]]) ] = [ 'settings','lj', float(fields[4]), float(fields[5]) ]
                        initial_params[ (data2map_dict[fields[2]],data2map_dict[fields[1]]) ] = [ 'settings','lj', float(fields[4]), float(fields[5]) ]

    return initial_params

if __name__ == "__main__":
   main(sys.argv[1:])
