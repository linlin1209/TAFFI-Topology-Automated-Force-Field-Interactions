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
    
    parser = argparse.ArgumentParser(description='Molecular configurations from a LAMMPS trajectory and generates an input for an orca single-point energy calculations.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('traj_file', help = 'A lammps trajectory file (usually outputted from a vdw_self_gen.py job).')

    parser.add_argument('map_file', help = 'The map file holding the correspondence between the atoms in the data file and the atom types in the VDW file (usually outputted from vdw_self_gen.py).')

    #optional arguments
    parser.add_argument('-start', dest='frame_start', default=0,
                        help = 'Frame to start parsing from in each trajectory (first=0; default: 0)')

    parser.add_argument('-end', dest='frame_end', default=0,
                        help = 'Frame to stop parsing from the trajectory (zero-indexing, like the start option; 0 defaults to last frame;  default: last)')

    parser.add_argument('-every', dest='every', default=1,
                        help = 'Every N frames will be parsed, where N is this variable (i.e. if mod(frame_num,every) == 0 then it is parsed) (parse all = 1; default: 1)')

    parser.add_argument('-N', dest='N_configs', default=100,
                        help = 'Holds the number of configurations to be generated (for each pair type) before exiting (default: 100)')

    parser.add_argument('-q', dest='q', default = 0,
                        help = 'Holds the charge of molecules being sampled. (elementary charge; default: 0)')

    parser.add_argument('-m', dest='m', default = 1,
                        help = 'Holds the multiplicity of the molecules being sampled. (default: 1)')

    parser.add_argument('-mol_list', dest='mol_list', default=[],
                        help = 'Holds a list of mol types, based on the trajectory file, included for single-point configurational averaging. '+\
                               'IDs are based on the output of the find_molecules subroutine, its output can be read using the --print_mol option. (by default all molecules are included ')

    parser.add_argument('-p', dest='procs', default=8,
                        help = 'holds the number of processors to run each orca job on (default: 8)')

    parser.add_argument('-o', dest='filename', default="CHELPG_jobs",
                        help = 'holds the name of the folder where the logfile and output data are stored. (default: vdw_orca_jobs)')
    
    parser.add_argument('-f', dest='functional', default='wB97X-D3',
                        help = 'Sets the functional for the calculation. (default: wB97X-D3; other typical options are B3LYP, and M062X)')
    
    parser.add_argument('--no_D3', dest='D3_option', default="D3BJ ", action='store_const', const="",
                        help = 'When this flag is present, the D3-dispersion correction is disabled. By default the dispersion correction is used.')

    parser.add_argument('--write_unwrapped', dest='write_unwrapped', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the unwrapped coordinates of each parsed frame are saved. Useful for debugging. (default: off)')

    parser.add_argument('-QC_types', dest='QC_types', default='dft mp2',
                        help = 'This variable controls whether dft, mp2, or both types of jobs are initialized. The variables expects a case-insensitive '+\
                               'space-delimited string (default: "dft mp2")')

    parser.add_argument('--print_mol', dest='print_mol', default=0, action='store_const', const=1,
                        help = 'When this flag is present, the unique molecule info is parsed from the simulation data and printed to screen. The molecule ids printed here can then be ' +\
                               'used for setting the A and B list arguments. (default: off)')
    

    # Assign stdin arguments to variables
    args=parser.parse_args(argv)
    args.frame_start = int(args.frame_start)
    args.frame_end   = int(args.frame_end)
    args.every       = int(args.every)
    args.filename    = str(args.filename)
    args.N_configs   = int(args.N_configs)
    args.q           = int(float(args.q))
    args.m           = int(float(args.m))
    args.procs       = int(args.procs)
    args.QC_types    = [ str(i).lower() for i in args.QC_types.split() ]
    if len(args.QC_types) == 0: print("ERROR: -QC_types requires at least one argument ('dft' 'mp2' or 'dft mp2'). Exiting..."); quit()
    if len([ i for i in args.QC_types if i not in ['dft','mp2']]) > 0: print("ERROR: -QC_types only accepts 'mp2' or 'dft' as arguments. Exiting..."); quit()
    if args.mol_list != []:
        args.mol_list = parse_mol_list(args.mol_list)

    # Check to ensure that data won't be overwritten
    if os.path.exists(args.filename):
        print('ERROR: A folder named {} already exists in this directory. Exiting to avoid overwriting data...'.format(args.filename))
        quit()
    elif args.print_mol == 0:
        os.makedirs(args.filename)
        os.makedirs(args.filename+'/configs/')     
        sys.stdout = Logger(args.filename)
        print("PROGRAM CALL: python gen_jobs_for_charges.py {}\n".format(' '.join([ i for i in argv])))

    # Check that the input trajectory exists
    if os.path.isfile(args.traj_file) == False:
        print("ERROR: Couldn't find {}. Double-check the name and location of the trajectory file. Exiting...".format(args.traj_file))
        quit()

    # Check that the files are properly formatted
    status,types,mols,num_atoms = check_dump_format(args.traj_file)

    if status==1:
        problems = [ Names[count_i] for count_i,i in enumerate(status) if i == 1 ]
        print("{}".format("*"*80))
        print("* {:76s} *".format("ERROR: the parser requires *sorted* *atom* dump types with the mol id in the 8th"))
        print("* {:76s} *".format("       column. Scripts exist to sort the dump files during post-processing."))
        print("{}".format("*"*80))        
        print("\nThe following file(s) had formatting errors:")
        for i in problems:
            print("\t{}".format(i))
        print("")
        quit()

    # Parse the unique molecules from the simulation data
    molecules,mol_id2type = parse_molecules(args.map_file)

    # If the user just wanted to see the molecule info then the program is exited here.
    if args.print_mol == 1:
        quit()

    # Auto-populate mol_list 
    if args.mol_list == []:
        args.mol_list = [ str(i) for i in list(mol_id2type.keys()) ]
    else:
        tmp = []
        for i in args.mol_list:            
            if i > len(molecules):
                print("ERROR: mol number {} was in mol_list, but not found in the list of unique molecules. Exiting...".format(i))
                quit()
            else:
                tmp += [ str(j) for j in list(mol_id2type.keys()) if mol_id2type[j] == i ]
        args.mol_list = tmp

    # Look up the mapping between the datafile and the mapfile
    map2data_dict,data2map_dict,data2elem_dict,data2charge_dict = Generate_Map(args.traj_file,args.map_file)

    # Generate adjacency matrix
    adj_mat =  Read_Adj_Mat(args.map_file)
    
    ###################################################
    # Parse the configurations from LAMMPS trajectory #
    ###################################################
    print("\n{}".format("*"*104))
    print("* {:^100s} *".format("Parsing configurations from LAMMPS trajectory "))
    print("{}\n".format("*"*104))

    # Get the number of frames in each trajectory
    num_frames = get_num_frames(args.traj_file)

    # Appropriately set the termination frame for each trajectory
    if args.frame_end == 0 or num_frames-1 < args.frame_end: # option for parsing all frames in all trajectories, or if less frames than the requested are present parse what's there. 
            args.frame_end = num_frames-1                    # set last frame to the total number of frames (-1 for zero indexing)

    # Parse configurations
    config_count = 0
    complete_flag = 0
    print("{:45} Configurations_Drawn".format(''))
    for frame in range(args.frame_start,args.frame_end+1,args.every):

        # Print diagnostic
        print("parsing configurations from frame {:<10}: {}".format(frame,config_count))

        # Grab frame data (num_atoms allows more efficient initialization of the parser arrays)
        ids,types,geo,mol,box = get_frame(frame,args.traj_file,num_atoms)

        # Unwrap molecules
        geo = unwrap(geo,adj_mat,mol,box)
        
        # Write unwrapped coords
        if args.write_unwrapped == 1:
            with open(args.filename+'/unwrapped_frames.xyz'.format(frame),'a') as f:
                f.write("{}\n\n".format(len(geo)))
                for count_i,i in enumerate(geo):
                    f.write(" {:20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n".format(data2elem_dict[types[count_i]],i[0],i[1],i[2]))
        
        # Grab relevant molecules
        red_list = []
        tmp_count = 0
        for count_i,i in enumerate(args.mol_list):

            # Intialize subgeometries for each molecule
            mol_inds     = [ count_k for count_k,k in enumerate(mol) if k == i ]
            mol_geo      = geo[mol_inds,:]

            # Generate element lists
            mol_elem = [ data2elem_dict[types[m]] for m in mol_inds ]

            # Generate type lists
            mol_types = [ data2map_dict[types[m]] for m in mol_inds ]

            # Generate charge lists
            mol_charges = [ data2charge_dict[types[m]] for m in mol_inds ] 

            # Generate xyz with geometry, types, mol_ids, and charges
            gen_xyz(args.filename+'/configs/'+str(config_count)+'/',config_count,mol_elem,mol_geo,mol_charges,mol_types)

            # Generate orca input file
            gen_orca(args.filename+'/configs/'+str(config_count)+'/',config_count,mol_elem,mol_types,args.functional,args.D3_option,mol_geo,args.q,args.m,args.procs,args.QC_types)

            # Iterate configuration counter
            config_count += 1

            # Check for break condition (inner loop over molecules)
            if config_count >= args.N_configs: 
                complete_flag = 1;
                break

        # Check for break condition (outer loop over frames)
        if complete_flag == 1:
            break

    print("\nTotal configurations parsed: {}".format(config_count))
                
    if complete_flag == 1:
        print("\n{}".format("*"*104))
        print("* {:^100s} *".format("Parse Complete! The configurations are waiting in {}".format(args.filename)))
        print("{}\n".format("*"*104))
    else:
        print("\n{}".format("*"*104))
        print("* {:^100s} *".format("Ran out of trajectory data before finding enough useful"))
        print("* {:^100s} *".format("configurations. Try running more MD or changing the -every setting."))
        print("* {:^100s} *".format("The configurations that were found are in {}".format(args.filename)))                                    
        print("{}\n".format("*"*104))


    return
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

# A function for writing orca input files for the charge calculation
def gen_orca(Save_Folder,config_count,Elements,Types,functional,D3_option,Geometry,charge=0,multiplicity=1,procs=8,QC_types=['dft','mp2']):

    # Print diagnostic message
    Save_Name = Save_Folder+str(config_count)+'.in'

    # Check for Li Be He and H, so that the frozen core approximation can be avoided
    No_Frozen_List = ["Li","Be","He","H"]    
    No_Frozen_Flag = 0
    for i in No_Frozen_List:
        if ( len(Elements) == 1 and i in Elements ):
            No_Frozen_Flag = 1

    # Make subdirectory to hold the current vdw files
    if os.path.isdir(Save_Folder) == False:
        os.makedirs(Save_Folder)

    # Open input file
    with open(Save_Name,'w') as f:

        #################################
        # Run DFT and MP2 single-points #
        #################################

        # Write DFT Header
        if 'dft' in QC_types:

            # Avoid dispersion correction
            if No_Frozen_Flag == 1:
                if procs == 1:
                    if charge < 0:
                        f.write('#Run DFT single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF CHELPG\n\n%base "{}_charges"'.format(functional,config_count))
                    else:
                        f.write('#Run DFT single point of the dimer\n! {} def2-TZVP TIGHTSCF CHELPG\n\n%base "{}_charges"'.format(functional,config_count))
                else:
                    if charge < 0:
                        f.write('#Run DFT single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF CHELPG PAL{}\n\n%base "{}_charges"'.format(functional,procs,config_count))
                    else:
                        f.write('#Run DFT single point of the dimer\n! {} def2-TZVP TIGHTSCF CHELPG PAL{}\n\n%base "{}_charges"'.format(functional,procs,config_count))
                
            # Use dispersion correction
            # If --noD3 flag is on then force no dispersion correction (this is for wB97X-D3)
            else:
                if procs == 1:
                    if charge < 0:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF {} CHELPG\n\n%base "{}_charges"'.format(functional,D3_option,config_count))
                    else:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} def2-TZVP TIGHTSCF {} CHELPG\n\n%base "{}_charges"'.format(functional,D3_option,config_count))
                else:
                    if charge < 0:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} ma-def2-TZVP TIGHTSCF {} CHELPG PAL{}\n\n%base "{}_charges"'.format(functional,D3_option,procs,config_count))
                    else:
                        f.write('#Run DFT-D3 single point of the dimer\n! {} def2-TZVP TIGHTSCF {} CHELPG PAL{}\n\n%base "{}_charges"'.format(functional,D3_option,procs,config_count))

            # Increase default max scf iterations
            f.write('\n\n%scf\nMaxIter 1000\nend\n')

            # Add print statements for different charge definitions
            f.write('%output\n')
            f.write('  {}\n  {}\n  {}\n  {}\n  {}\nend\n'.format('Print[ P_Mayer ] 1','Print[ P_NatPop ] 1','Print[ P_Hirshfeld ] 1', 'Print[ P_Mulliken ] 1','Print [ P_Loewdin ] 1'))

            # Write Coordinate header
            f.write('\n* xyz {} {}\n'.format(charge,multiplicity))

            # Write Geometry 
            for count_j,j in enumerate(Geometry):
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

        # Write RI-MP2 single point calculation on the dimer
        if 'mp2' in QC_types:
            if procs == 1:
                if No_Frozen_Flag == 1:
                    if charge < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(i))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore\n'.format(i))
                else:
                    if charge < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG\n'.format(i))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG\n'.format(i))
            else:
                if No_Frozen_Flag == 1:
                    if charge < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore PAL{}\n'.format(i,procs))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG NoFrozenCore PAL{}\n'.format(i,procs))
                else:
                    if charge < 0:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 aug-cc-pVTZ aug-cc-pVTZ/C TIGHTSCF CHELPG PAL{}\n'.format(i,procs))
                    else:
                        f.write('\n#Run RI-MP2 single point of the dimer\n$new_job\n%base "{}_RI-MP2_AB"\n! RI-MP2 cc-pVTZ cc-pVTZ/C TIGHTSCF CHELPG PAL{}\n'.format(i,procs))

            # Increase default max scf iterations
            f.write('\n%scf\nMaxIter 1000\nend\n')

            # Add print statements for different charge definitions
            f.write('%output\n')
            f.write('  {}\n  {}\n  {}\n  {}\n  {}\nend\n'.format('Print[ P_Mayer ] 1','Print[ P_NatPop ] 1','Print[ P_Hirshfeld ] 1', 'Print[ P_Mulliken ] 1','Print [ P_Loewdin ] 1'))

            # Write coordinate header
            f.write('\n* xyz {} {}\n'.format(charge,multiplicity))

            # Write Geometry 
            for count_j,j in enumerate(Geometry):
                f.write('  {:<20s}   {: <20.6f} {: <20.6f} {: <20.6f}\n'.format(Elements[count_j],j[0],j[1],j[2]))
            f.write('*\n')

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

# Function checks the formatting of the LAMMPS dump file
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
        print("ERROR: the program expects an 'atom' formatted lammps trajectory with the atom id, atom type,")
        print("       molecule id, and cartesian coordinates for each atom. Exiting...")
        quit()        

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

# Generates the adjacency matrix from a map file
# XXX consider turning one of the incarnations of this function into a Lib function
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

# XXX Replace this function with the corresponding Lib function 
def gen_xyz(Save_Folder,config_count,Elements,Geometry,Charges,Types):

    # Print diagnostic message
    Save_Name = Save_Folder+str(config_count)+'.in'

    # Make subdirectory to hold the current vdw files
    if os.path.isdir(Save_Folder) == False:
        os.makedirs(Save_Folder)

    # Save an xyz for viewing the configuration in VMD
    with open(Save_Folder+str(config_count)+'.xyz','w') as f:
        f.write('{}\n\n'.format(len(Geometry)))
        for count_j,j in enumerate(Geometry):
            f.write('  {:<10s} {:< 20.6f} {:< 20.6f} {:< 20.6f} {:<12d} {:< 12.6f} {:<60s}\n'.format(Elements[count_j],j[0],j[1],j[2],0,Charges[count_j],Types[count_j]))

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
        self.log = open(folder+"/gen_jobs_for_charges.log", "a")

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
def parse_molecules(mapfile):

    # Print banner
    print("\n{}".format("*"*104))
    print("* {:^100s} *".format("Parsing unique molecules from the simulation data"))
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

    # 


    molecules = []
    for count_i,i in enumerate(atomtype_list):
        molecules += [(i,charges_list[count_i],adj_mat_list[count_i])]

    # Print diagnostic and return
    print("Number of unique molecules in the fit configurations: {}".format(len(molecules)))

    for count_i,i in enumerate(molecules):
        print("\nMolecule: {}".format(count_i))
        print("\n\t{:<40s}  {:<40s}".format("Atom_types","Charges"))
        for count_j,j in enumerate(i[0]):
            print("\t{:<40s} {:< 12.6f}".format(j,i[1][count_j]))
        print("\n\tCorrespondence with trajectory mol_id: {}".format(', '.join([ str(j) for j in list(mol_id2type.keys()) if mol_id2type[j]==count_i ])))
        print("\tTotal charge: {:< d}".format(int(round(sum(i[1])))))

    return molecules,mol_id2type

if __name__ == "__main__":
   main(sys.argv[1:])
