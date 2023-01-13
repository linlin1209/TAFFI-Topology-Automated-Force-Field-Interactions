#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *
from gen_md_for_sampling import initialize_VDW, write_map, write_molecule
#from gen_md_for_drude import get_data, Write_data, Pack_Box
from gen_md_for_drude import get_data, Pack_Box
from gen_md_for_drude_DCE import Write_data

# NOTE: the order of the imports matters until I resolve the namespace issues (XXX DON'T USE import *)
import random
from numpy import *
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file and writes inputs for a LAMMPS job for a single molecule shrink-wrapped simulation.'
                                                 'These simulations are meant to be used in conjunction with free-energy trajectories for subtracting out the intramolecular '+\
                                                 'electrostatics from the final results.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'A quoted list of the *.xyz files to be included in the sim. '+\
                                              ' (currently must be an xyz with the atom types in the fourth column.')

    parser.add_argument('FF_files', help = 'A quoted list of the force-field files. '+\
                                           'Formatting of the force-field files is as produced by the FF_gen.py/FF_extract.py programs)')

    #optional arguments    
    parser.add_argument('-N', dest='N_mol', default="1",
                        help = 'Controls the number of molecules put in the simulation cell. The program expects a list of integers. If less integers '+\
                               'are supplied than the number of coord files, then the list is automatically expanded using the last supplied (or default) '+\
                               'element. A single integer rather than a list is also accepted. (default: 1, i.e. 1 of each molecule supplied are placed in the simulation cell)')

    parser.add_argument('-T', dest='T_equil', default=400,
                        help = 'Controls the temperature of the MD simulation. (in Kelvin; default: 400)')

    parser.add_argument('-t', dest='t_equil', default=1E6,
                        help = 'Controls the length of the MD equilibration segment. (default: 1E6)')

    parser.add_argument('-T_A', dest='T_anneal', default=400,
                        help = 'Controls the temperature of the MD simulation. (in Kelvin; default: 400)')

    parser.add_argument('-t_A', dest='t_anneal', default=1E6,
                        help = 'Controls the length of the MD equilibration segment. (default: 1E6)')

    parser.add_argument('-d', dest='density', default=0,
                        help = 'Sets the density for the MD simulation. '+\
                               'This option is incompatible with the -d_N option, only one can be supplied at a time (default: 0, unfixed, units of g/cubic cm)')

    parser.add_argument('-d_N', dest='N_density', default=0,
                        help = 'Sets the number density for the MD simulation. '+\
                               'This option is incompatible with the -d option, only one can be supplied at a time (default: 0, unfixed; units of atoms/cubic Angstrom)')

    parser.add_argument('-f', dest='frequency', default=1000,
                        help = 'Controls the sampling frequency during the MD simulation. (default: 1000)')

    parser.add_argument('-q', dest='q_list', default="none",
                        help = 'Controls the total charge on the molecule. By default, the rounded integer charge is used for each molecule (round). The program expects a list of integers (or none, or round for the default). '+\
                               'If less integers are supplied than the number of coord files, then the list is automatically expanded using the last supplied (or default) element. "none" makes no adjustment to the charges '+\
                               '(default: None)') 

    parser.add_argument('-o', dest='outputname', default='Shrink_traj',
                        help = 'Sets the output filename prefix (default: Shrink_traj)')

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'When atomtypes are being automatically determined (i.e. when the supplied *.xyz files do not have atomtype information) this variable controls the bond depth used to define unique atomtypes. (default: 2)')

    parser.add_argument('-onefourscale', dest='onefourscale', default=0.0,
                        help = 'Sets scaling for 1-4 VDW and Electrostatic interactions (default: 0.0)')

    parser.add_argument('-eps_scale', dest='eps_scale', default=1.0,
                        help = 'Sets scaling for the default UFF eps parameters. To increase configurational sampling, it can be advantageous '+\
                               'to reduce the stiffness of the VDW potentials. (default: 1.0; full value)')

    parser.add_argument('-sigma_scale', dest='sigma_scale', default=1.0,
                        help = 'Sets scaling for the default UFF sigma parameters. To increase configurational sampling, it can be advantageous '+\
                               'to reduce the LJ minima of the VDW potentials. (default: 1.0; full value)')

    parser.add_argument('-charge_scale', dest='charge_scale', default=1.0,
                        help = 'Sets scaling for the simulation charges. To increase configurational sampling, it can be advantageous '+\
                               'to reduce the charges (especially for ions). (default: 1.0; full value)')

    parser.add_argument('-mixing_rule', dest='mixing_rule', default="none",
                        help = 'Defines the mixing rule to be used for missing LJ parameters. Waldman-Hagler (wh) and Lorentz-Berthelot (lb) and "none", are valid options. When set to "none", the program will print a message and '+\
                               'exit if there are missing parameters. (default: "none")')

    parser.add_argument('-pair_styles', dest='pair_styles', default='lj/cut/coul/cut 15.0 15.0',
                        help = 'Supplies the information for the pair style setting(s). Supplied as a string and should be formatted for comply '+\
                               'with the LAMMPS "pair_style hybrid" command and should match the cutoffs for the soft potential used in the free energy trajectory '+\
                               '(default: "lj/cut/coul/cut 14.0 14.0")')

    parser.add_argument('--impropers', dest='improper_flag', default=False, action='store_const', const=True,
                        help = 'When present, the simulation will use improper dihedral terms. (default: False)')

    parser.add_argument('--UFF', dest='force_UFF', default=0, action='store_const', const=1,
                        help = 'Forces the use of UFF parameters (regardless of the presence of corresponding parameters in the read FF files). (default: off)')

    # Make parse inputs
    args=parser.parse_args(argv)

    # Converting each list input argument. 
    args.coord_files = args.coord_files.split()
    args.FF_files = args.FF_files.split()
    args.N_mol = [ int(i) for i in args.N_mol.split() ]
    args.q_list = [ str(i) for i in args.q_list.split() ]

    # Extend N_mol and q_list to match the length of coord_file
    while len(args.N_mol) < len(args.coord_files): args.N_mol = args.N_mol + [args.N_mol[-1]] 
    while len(args.q_list) < len(args.coord_files): args.q_list = args.q_list + [args.q_list[-1]] 

    # Parse input and make N values global 
    if type(args.N_mol) != type([]): args.N_mol = ast.literal_eval(args.N_mol)

    args.gens = int(args.gens)
    args.T_equil = float(args.T_equil)
    args.T_anneal = float(args.T_anneal)
    args.t_equil = int(float(args.t_equil))
    args.t_anneal = int(float(args.t_anneal))
    args.frequency = int(float(args.frequency))
    args.eps_scale = float(args.eps_scale)
    args.sigma_scale = float(args.sigma_scale)
    args.charge_scale = float(args.charge_scale)
    args.density = float(args.density)
    args.N_density = float(args.N_density)
    args.mixing_rule = str(args.mixing_rule).lower()
    if args.density != 0 and args.N_density != 0: print("ERROR: Setting both the -d and -d_N options is inconsistent (both specify density). Please revise. Exiting..."); quit()
    if args.mixing_rule not in ["wh","lb","none"]:
        print("ERROR: '{}' is not a valid argument for -mixing_rule. Valid options are 'wh' 'lb' and 'none'. Exiting...".format(args.mixing_rule))
        quit()

    # Check that the input is an .xyz file. 
    for i in args.coord_files:
        if i.split('.')[-1] != 'xyz':
            print("ERROR: Check to ensure that the input coordinate file(s) are in .xyz format.")
            quit()

    # Generate filename and directory for output
    # Make directory to hold the output
    if args.outputname:
        Filename = args.outputname
    else:
        Filename = args.coord_file[0].split('.')[0]

    # Make the output directory if it doesn't already exist.
    if os.path.exists(Filename):
        print('ERROR: The desired output folder ({}) already exists. Exiting to avoid overwriting data...'.format(Filename))
        return
    else:
        os.makedirs(Filename)
        sys.stdout = Logger(Filename)
        print("PROGRAM CALL: python gen_mol_traj.py {}\n".format(' '.join([ i for i in argv])))
        
    # Check that the supplied FF files exist
    for i in args.FF_files:
        if os.path.isfile(i) != True:
            print("ERROR: Specified FF_file ({}) doesn't exist. Please check the path. Exiting...".format(i))
            quit()

    # Catenate the FF files and copy them to the output folder
    FF_all = Filename+"/"+Filename.split('/')[-1]+'.db'
    with open(FF_all,'w') as f:
        for i in args.FF_files:
            with open(i,'r') as ff_file:
                for j in ff_file:
                    f.write(j)
            f.write("\n")

    # Grab data on each molecule being added to the MD simulation.
    Data = get_data(FF_all,args.coord_files,args.N_mol,args.q_list,args.gens,Improper_flag=args.improper_flag)
    alpha_dict = {}
    for i in Data:
         print(Data[i]['Charges'])
         for count_j,j in enumerate(Data[i]['Atom_types']):
            alpha_dict[j] = Data[i]['alphas'][count_j]

    # Pack the simulation cell
    Elements,Atom_types,Geometry,Bonds,Bond_types,Angles,Angle_types,Dihedrals,Dihedral_types,Impropers,Improper_types,Charges,Molecule,Molecule_names,Adj_mat,Sim_Box,b_min,alphas = \
    Pack_Box(Data,Box_offset=0,Density=args.density,N_Density=args.N_density)

    # Generate VDW parameters    
    VDW_params = initialize_VDW(sorted(set(Atom_types)),sigma_scale=args.sigma_scale,eps_scale=args.eps_scale,VDW_FF=Data[args.coord_files[0]]["VDW_params"],\
                                Force_UFF=args.force_UFF,mixing_rule=args.mixing_rule)    

    # Update pair style if supplied by the user
    if args.pair_styles is not None:
        style = args.pair_styles.split()[0]
        for i in list(VDW_params.keys()):
            VDW_params[i] = tuple([style]+VDW_params[i][1:])

    # Generate Simulation Dictionaries
    # The bond, angle, and diehdral parameters for each molecule are combined into one dictionary
    Bond_params = {}; Angle_params = {}; Dihedral_params = {}; Improper_params = {}; Masses = {}
    for i in list(Data.keys()):
        for j in list(Data[i]["Bond_params"].keys()): Bond_params[j] = Data[i]["Bond_params"][j]
        for j in list(Data[i]["Angle_params"].keys()): Angle_params[j] = Data[i]["Angle_params"][j]
        for j in list(Data[i]["Dihedral_params"].keys()): Dihedral_params[j] = Data[i]["Dihedral_params"][j]
        for j in list(Data[i]["Improper_params"].keys()): Improper_params[j] = Data[i]["Improper_params"][j]
        for j in list(Data[i]["Masses"].keys()): Masses[j] = Data[i]["Masses"][j]
        
    # Write the lammps datafile (the function returns a dictionary, Atom_type_dict, holding the mapping between atom_types and the lammps id numbers; this mapping is needed for setting fixes)
    print("Writing LAMMPS datafile ({})...".format(Filename+'.data'))
    print("Writing LAMMPS settings file ({})...".format(Filename+'.in.settings'))
    # ?????
    Atom_type_dict,fixed_modes,drude_tag = Write_data(Filename,Atom_types,Sim_Box,Elements,Geometry,Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,\
                                            Dihedrals,Dihedral_types,Dihedral_params,Impropers,Improper_types,Improper_params,Charges,VDW_params,Masses,Molecule,args.N_mol[0],alphas,alpha_dict,Improper_flag = args.improper_flag)
    
    # Generate fix statement for scaled charges
    fixes = scale_charges(args.charge_scale,Atom_type_dict,Atom_types,Charges)

    # Gather the different dihedral_styles
    if len(Dihedrals) > 0:
        Dihedral_styles = set([ Dihedral_params[i][0] for i in list(Dihedral_params.keys()) ])
    else:
        Dihedral_styles = []

    # Write the lammps input files
    print("Writing LAMMPS input file ({})...".format(Filename+'.in.init'))
    Write_input(Filename,args.T_equil,args.t_equil,args.T_anneal,args.t_anneal,args.frequency,args.onefourscale,args.pair_styles,Dihedral_styles,args.density,args.N_density,fixes,drude_tag)

    # Write map file to more easily post-process the trajectory
    print("Writing mapfile ({})...".format(Filename+'.map'))
    write_map(Filename,Elements,Atom_types,Charges,Masses,Adj_mat,zeros([len(Elements)]),args.N_mol)

    # Writing molecule map file
    print("Writing molfile ({})...".format(Filename+'.mol.txt'))
    write_molecule(Filename,Molecule_names)
    
    # Remove concatenated FF file from run folder
    os.remove(FF_all)

    # Print banner
    print("\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Success! Have a Nice Day!","*"*167))

# Description: A wrapper for the write commands for generating the lammps input file
def Write_input(Filename,equil_temp,equil_time,anneal_temp,anneal_time,frequency,Onefourscale,Pair_styles,Dihedral_styles,Density,N_Density,fixes,drude_tag):

    with open(Filename+'/'+Filename.split('/')[-1]+'.in.init','w') as f:
        f.write("# lammps input file for polymer simulation with dilute ions\n\n"+\
                "# VARIABLES\n"+\
                "variable        data_name       index   {}\n".format(Filename.split('/')[-1]+'.data')+\
                "variable        settings_name   index   {}\n".format(Filename.split('/')[-1]+'.in.settings')+\
                "variable        log_name        index   {}\n".format(Filename.split('/')[-1]+'.log')+\
                "variable        nSteps_ramp     index   {} # number of data steps for the ramped anneal\n".format(int(anneal_time/2))+\
                "variable        nSteps_equil    index   {} # number of data steps for the ramped anneal\n".format(equil_time)+\
                "variable        avg_freq        index   1000\n".format(frequency)+\
                "variable        min_freq        index   1\n"+\
                "variable        coords_freq     index   1000\n".format(frequency)+\
                "variable        thermo_freq     index   1000\n".format(frequency)+\
                "variable        dump4avg        index   100\n"+\
                "variable        vseed           index   {: <6d}\n".format(int(random.random()*100000))+\
                "variable        ANNEAL_TEMP     index   {} # Temperature during the initial anneal\n".format(anneal_temp)+\
                "variable        FINAL_TEMP      index   {} # Temperature ramped to during the final anneal\n".format(equil_temp)+\

                "# Change the name of the log output #\n"+\
                "log ${log_name}\n\n"+\
                
                "#===========================================================\n"+\
                "# GENERAL PROCEDURES\n"+\
                "#===========================================================\n"+\
                "units		real	# g/mol, angstroms, fs, kcal/mol, K, atm, charge*angstrom\n"+\
                "dimension	3	# 3 dimensional simulation\n"+\
                "newton		off	# use Newton's 3rd law\n"+\
                "boundary	s s s	# periodic boundary conditions \n"+\
                "atom_style	full    # molecular + charge\n\n"+\

                "#===========================================================\n"+\
                "# FORCE FIELD DEFINITION\n"+\
                "#===========================================================\n"+\
                'special_bonds  lj   0.0 0.0 0.0     # NO     1-4 LJ interactions\n'+\
                'special_bonds  coul 0.0 0.0 {}     # REDUCE 1-4 electrostatics\n'.format(Onefourscale)+\
                'pair_style     hybrid/overlay {} thole 2.6 14.0 # outer_LJ outer_Coul (cutoff values, see LAMMPS Doc)\n'.format(Pair_styles)+\
                'bond_style     harmonic             # parameters needed: k_bond, r0\n'+\
                'angle_style    harmonic             # parameters needed: k_theta, theta0\n')
        if len(Dihedral_styles) >= 1:
            f.write('dihedral_style hybrid {}            # parameters needed: V1, V2, V3, V4\n'.format(' '.join(Dihedral_styles)))            
        f.write('pair_modify    shift yes mix arithmetic       # using Lorenz-Berthelot mixing rules\n\n'+\

                "#===========================================================\n"+\
                "# SETUP SIMULATIONS\n"+\
                "#===========================================================\n\n"+\
                "# READ IN COEFFICIENTS/COORDINATES/TOPOLOGY\n"+\
                'read_data ${data_name} extra/special/per/atom 1\n'+\
                'include ${settings_name}\n\n'+\

                '# INVOKE Drude\n'+\
                'fix DRUDE all drude {}\n\n'.format(" ".join(drude_tag))+\

                "# SET RUN PARAMETERS\n"+\
                "timestep	1.0		# fs\n"+\
                "run_style	verlet 		# Velocity-Verlet integrator\n"+\
                "neighbor 2.0 nsq\n"+\
                "neigh_modify every 1 delay 0 check no # More relaxed rebuild criteria can be used\n\n"+\
                "# atom groups convenient for thermostats (see package documentation), etc.\n"+\
                "group ATOMS type {}\n".format(" ".join([str(count_i+1) for count_i,i in enumerate(drude_tag) if i != 'D']))+\
                "group CORES type {}\n".format(" ".join([str(count_i+1) for count_i,i in enumerate(drude_tag) if i == 'C']))+\
                "group DRUDES type {}\n".format(" ".join([str(count_i+1) for count_i,i in enumerate(drude_tag) if i == 'D']))+\
                "group NotDrude type {}\n\n".format(" ".join([str(count_i+1) for count_i,i in enumerate(drude_tag) if i != 'D'])))
        
        if fixes != []:
            f.write("{}".format(fixes))

        f.write("#===========================================================\n"+\
                "# RUN CONSTRAINED RELAXATION\n"+\
                "#===========================================================\n\n"+

                "# SET OUTPUTS\n"+\
                "compute TDRUDE all temp/drude\n"+\
                "thermo_style    custom step temp c_TDRUDE[1] c_TDRUDE[2] vol density etotal pe ebond eangle edihed ecoul elong evdwl enthalpy press\n"+\
                "thermo_modify   format float %14.6f\n"+\
                "thermo ${thermo_freq}\n\n"+\

                "# DECLARE RELEVANT OUTPUT VARIABLES\n"+\
                "variable        my_step equal   step\n"+\
                "variable        my_temp equal   temp\n"+\
                "variable        my_temp_true equal    c_TDRUDE[1]\n"+\
                "variable        my_temp_drude equal   c_TDRUDE[2]\n"+\
                "variable        my_rho  equal   density\n"+\
                "variable        my_pe   equal   pe\n"+\
                "variable        my_ebon equal   ebond\n"+\
                "variable        my_eang equal   eangle\n"+\
                "variable        my_edih equal   edihed\n"+\
                "variable        my_evdw equal   evdwl\n"+\
                "variable        my_eel  equal   (ecoul+elong)\n"+\
                "variable        my_ent  equal   enthalpy\n"+\
                "variable        my_P    equal   press\n"+\
                "variable        my_vol  equal   vol\n\n"+\

                "fix  averages all ave/time ${dump4avg} $(v_avg_freq/v_dump4avg) ${avg_freq} v_my_temp v_my_temp_true v_my_temp_drude v_my_rho v_my_vol v_my_pe v_my_edih v_my_evdw v_my_eel v_my_ent v_my_P file thermo.anneal.avg\n\n"+\

                "#===========================================================\n"+\
                "# RUN DRUDE RELAXATION\n"+\
                "#===========================================================\n\n"+
                "fix LANG all langevin/drude ${FINAL_TEMP} 100 12435 1. 20 13977\n"+\
                "comm_modify vel yes\n"+\
                "velocity        ATOMS create ${FINAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n"+\
                "velocity        DRUDES create 1. ${vseed} mom yes rot yes     # DRAW VELOCITIES\n"+\
                "fix mom all momentum 1000 linear 1 1 1 angular # Zero out system linear and angular momentum every ps\n"+\
                "fix MIN all nve\n"+\
                "#fix FREEZE NotDrude setforce 0.0 0.0 0.0\n"+\
                "dump MIN all custom ${min_freq} min.lammpstrj id type x y z mol\n"+\
                "dump_modify MIN sort id format float %20.10g\n"+\
                "min_style sd\n"+\
                "minimize 1e-6 1e-9 1000 100000\n"+\
                "undump MIN\n"+\
                "unfix MIN\n"+\
                "#unfix FREEZE\n\n"+\
   
                #"# INITIALIZE VELOCITIES AND CREATE THE CONSTRATINED RELAXATION FIX\n"+\
                #"velocity        all create ${ANNEAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n\n"+\

                #"fix relax all nve/limit 0.1\n"+\
                #"run             1000\n"+\
                #"unfix relax\n\n"+\

                "#===========================================================\n"+\
                "# RUN RAMPED ANNEAL\n"+\
                "#===========================================================\n\n"+\

                "# REINITIALIZE THE VELOCITIES AND CREATE THE ANNEALING FIX\n"+\
                "velocity        ATOMS create ${FINAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n"+\
                "velocity        DRUDES create 1. ${vseed} mom yes rot yes\n")
        f.write("comm_modify vel yes\n")
        f.write("fix anneal all nve # NVT, langvin\n\n")

        f.write("# CREATE COORDINATE DUMPS FOR ANNEAL\n"+\
                "dump anneal all custom ${coords_freq} anneal.lammpstrj id type x y z mol\n"+\
                "dump_modify anneal sort  id\n\n"+\

                "# RUN RAMPED ANNEAL\n"+\
                "run		${nSteps_ramp}\n"+\
                "unfix anneal\n\n"+\

                "# RUN EQUILIBRATION PHASE AT FINAL TEMP\n")
        f.write("fix anneal all nve # NVT, langevin\n\n")
        f.write("run ${nSteps_ramp}\n"+\
                "unfix anneal\n"+\
                "unfix averages\n"+\
                "undump anneal\n\n"+\

                "#===========================================================\n"+\
                "# RUN EQUILIBRIUM SIM\n"+\
                "#===========================================================\n\n"+

                "# UPDATE RUN PARAMETERS AND CREATE FIX\n")
        f.write("comm_modify vel yes\n")
        f.write("fix equil all nve # NVT, langevin\n\n")

        f.write("# CREATE COORDINATE DUMPS FOR EQUILIBRIUM\n"+\
                "fix  averages all ave/time ${dump4avg} $(v_avg_freq/v_dump4avg) ${avg_freq} v_my_temp v_my_temp_true v_my_temp_drude v_my_rho v_my_vol v_my_pe v_my_edih v_my_evdw v_my_eel v_my_ent v_my_P file thermo.equil.avg\n\n"+\
                "dump equil all custom ${coords_freq} equil.lammpstrj id type x y z mol\n"+\
                "dump_modify equil sort  id\n\n"+\

                "# RUN EQUIL\n"+\
                "run		${nSteps_equil}\n"+\
                "unfix equil\n"+\
                "undump equil\n\n"+\
                    
                "# WRITE RESTART FILES, CLEANUP, AND EXIT\n"+\
                "write_restart   {}\n".format(Filename.split('/')[-1]+'.end.restart')+\
                "write_data      {} pair ii\n".format(Filename.split('/')[-1]+'.end.data')+\
                "unfix		averages\n")
    f.close()

    return

# Description: finds the number of disconnected subnetworks in the 
#              adjacency matrix, which corresponds to the number of 
#              separate molecules.
#
# Inputs:      adj_mat: numpy array holding a 1 in the indices of bonded
#                        atom types. 
#
# Returns:     mol_count: scalar, the number of molecules in the adj_mat
def mol_count(adj_mat):
    
    # Initialize list of atoms assigned to molecules and a counter for molecules
    placed_idx = []    
    mol_count = 0

    # Continue the search until all the atoms have been assigned to molecules
    while len(placed_idx)<len(adj_mat):

        # Use sequential elements of the adj_mat as seeds for the spanning network search
        for count_i,i in enumerate(adj_mat):

            # Only proceed with search if the current atom hasn't been placed in a molecule
            if count_i not in placed_idx:

                # Increment mol_count for every new seed and add the seed to the list of placed atoms
                mol_count += 1               
                placed_idx += [count_i]
                
                # Find connections
                idx = [ count_j for count_j,j in enumerate(i) if j==1 and count_j not in placed_idx ]
                
                # Continue until no new atoms are found
                while len(idx) > 0:
                    current = idx.pop(0)
                    if current not in placed_idx:
                        placed_idx += [current]
                        idx += [ count_k for count_k,k in enumerate(adj_mat[current]) if k == 1 and count_k not in placed_idx ]
    return mol_count

# Description:
# Rotate Point by an angle, theta, about the vector with an orientation of v1 passing through v2. 
# Performs counter-clockwise rotations (i.e., if the direction vector were pointing
# at the spectator, the rotations would appear counter-clockwise)
# For example, a 90 degree rotation of a 0,0,1 about the canonical 
# y-axis results in 1,0,0.
#
# Point: 1x3 array, coordinates to be rotated
# v1: 1x3 array, rotation direction vector
# v2: 1x3 array, point the rotation vector passes through
# theta: scalar, magnitude of the rotation (defined by default in degrees)
def axis_rot(Point,v1,v2,theta,mode='angle'):

    # Temporary variable for performing the transformation
    rotated=array([Point[0],Point[1],Point[2]])

    # If mode is set to 'angle' then theta needs to be converted to radians to be compatible with the
    # definition of the rotation vectors
    if mode == 'angle':
        theta = theta*pi/180.0

    # Rotation carried out using formulae defined here (11/22/13) http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/)
    # Adapted for the assumption that v1 is the direction vector and v2 is a point that v1 passes through
    a = v2[0]
    b = v2[1]
    c = v2[2]
    u = v1[0]
    v = v1[1]
    w = v1[2]
    L = u**2 + v**2 + w**2

    # Rotate Point
    x=rotated[0]
    y=rotated[1]
    z=rotated[2]

    # x-transformation
    rotated[0] = ( a * ( v**2 + w**2 ) - u*(b*v + c*w - u*x - v*y - w*z) )\
             * ( 1.0 - cos(theta) ) + L*x*cos(theta) + L**(0.5)*( -c*v + b*w - w*y + v*z )*sin(theta)

    # y-transformation
    rotated[1] = ( b * ( u**2 + w**2 ) - v*(a*u + c*w - u*x - v*y - w*z) )\
             * ( 1.0 - cos(theta) ) + L*y*cos(theta) + L**(0.5)*(  c*u - a*w + w*x - u*z )*sin(theta)

    # z-transformation
    rotated[2] = ( c * ( u**2 + v**2 ) - w*(a*u + b*v - u*x - v*y - w*z) )\
             * ( 1.0 - cos(theta) ) + L*z*cos(theta) + L**(0.5)*( -b*u + a*v - v*x + u*y )*sin(theta)

    rotated = rotated/L
    return rotated

# A simple wrapper function for generating a fix statement for scaled charges
def scale_charges(charge_scale,Atom_type_dict,Atom_types,Charges):

    fixes = []
    if charge_scale != 1.0:
        fixes = '# Setting up fixes to scale the charges by {}\n'.format(charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):
            fixes += 'group charge_{} type {}\n'.format(i,i)        
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):
            index = Atom_types.index([ k for k in list(Atom_type_dict.keys()) if Atom_type_dict[k] == i ][0])
            fixes += 'variable charge_{} equal {}\n'.format(i,Charges[index]*charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):        
            index = Atom_types.index([ k for k in list(Atom_type_dict.keys()) if Atom_type_dict[k] == i ][0])
            fixes += 'fix {} charge_{} adapt 0 atom charge v_charge_{} reset yes\n'.format(i,i,i)
        fixes += '\n'
    return fixes

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/gen_mol_traj.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

if __name__ == "__main__":
   main(sys.argv[1:])
