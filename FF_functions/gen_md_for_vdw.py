#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *

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

    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file (usually outputted from an intramolecular mode parse) '+\
                                                 'and writes inputs for a LAMMPS job to generate configurations for the parsing VDW parameters.')


    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')


    #optional arguments    
    parser.add_argument('-FF',dest='FF_files', default="",
                         help = 'A quoted list of the force-field files. Formatting of the force-field files is as produced by the FF_gen.py/FF_extract.py programs)')

    parser.add_argument('-N', dest='N_mol', default="25",
                        help = 'Controls the number of molecules put in the simulation cell. The program expects a list of integers. If less integers '+\
                               'are supplied than the number of coord files, then the list is automatically expanded using the last supplied (or default) '+\
                               'element. A single integer rather than a list is also accepted. (default: 25, i.e. 25 of each molecule supplied are placed in the simulation cell)')

    parser.add_argument('-T', dest='T_equil', default=400,
                        help = 'Controls the temperature of the MD simulation. (in Kelvin; default: 400)')

    parser.add_argument('-t', dest='t_equil', default=1E6,
                        help = 'Controls the length of the MD equilibration segment. (default: 1E6)')

    parser.add_argument('-T_A', dest='T_anneal', default=400,
                        help = 'Controls the temperature of the MD simulation. (in Kelvin; default: 400)')

    parser.add_argument('-t_A', dest='t_anneal', default=1E6,
                        help = 'Controls the length of the MD equilibration segment. (default: 1E6)')

    parser.add_argument('-P', dest='pressure', default=1,
                        help = 'Controls the pressure during the MD simulation. (in ATM; default: 1)')

    parser.add_argument('-d', dest='density', default=0,
                        help = 'Sets the density for the MD simulation. When this option is set, an NVT simulation is performed instead '
                               'of an NPT simulation, and the value supplied for pressure (-P) is disregarded. '+\
                               'This option is incompatible with the -d_N option, only one can be supplied at a time (default: 0, NPT simulation; units of g/cubic cm)')

    parser.add_argument('-d_N', dest='N_density', default=0,
                        help = 'Sets the number density for the MD simulation. When this option is set, an NVT simulation is performed instead '
                               'of an NPT simulation, and the value supplied for pressure (-P) is disregarded. '+\
                               'This option is incompatible with the -d option, only one can be supplied at a time (default: 0, NPT simulation; units of atoms/cubic Angstrom)')

    parser.add_argument('-fix', dest='fixes', default="",
                        help = 'Determines if bonds and/or angles are contstrained in the simulation. "bonds" and "angles" are values options that can be supplied as a space-delimited string. By default the rattle '+\
                               'algorithm is used for constraining the modes. (default: no constrained modes)')

    parser.add_argument('-f', dest='frequency', default=1000,
                        help = 'Controls the sampling frequency during the MD simulation. (default: 1000)')

    parser.add_argument('-o', dest='outputname', default='',
                        help = 'Sets the output filename prefix (default: VDW_selfterms)')

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

    parser.add_argument('--UFF', dest='force_UFF', default=0, action='store_const', const=1,
                        help = 'Forces the use of UFF parameters (regardless of the presence of corresponding parameters in the read FF files). (default: off)')

    parser.add_argument('-pair_styles', dest='pair_styles', default='lj/cut/coul/long 10.0 10.0',
                        help = 'Supplies the information for the pair style setting(s). Supplied as a string and should be formatted for comply '+\
                               'with the LAMMPS "pair_style hybrid" command (default: "lj/cut/coul/long 10.0 10.0")')

    # Make parse inputs
    args=parser.parse_args(argv)

    # Converting each list input argument. 
    args.coord_files = args.coord_files.split()
    args.FF_files = args.FF_files.split()
    args.N_mol = [ int(i) for i in args.N_mol.split() ]

    # Check if a directory was supplied to the coord_files list.
    bools = [ os.path.isdir(i) for i in args.coord_files ]
    if True in bools: 
        if len([ True for i in bools if i == True ]) > 1:
            print "ERROR: multiple directories were supplied to the script. Only one directory can\nbe supplied if the user intends to start the job from a previous run. Exiting..."
            quit()        
        args.coord_files = [args.coord_files[bools.index(True)]]
        read_flag = 1
    else:
        read_flag = 0

    # Extend N_mol to match the length of coord_file
    while len(args.N_mol) < len(args.coord_files): args.N_mol = args.N_mol + [args.N_mol[-1]] 

    # Parse input and make N values global 
    if type(args.N_mol) != type([]): args.N_mol = ast.literal_eval(args.N_mol)

    args.T_equil      = float(args.T_equil)
    args.T_anneal     = float(args.T_anneal)
    args.t_equil      = int(float(args.t_equil))
    args.t_anneal     = int(float(args.t_anneal))
    args.pressure     = float(args.pressure)
    args.frequency    = int(float(args.frequency))
    args.eps_scale    = float(args.eps_scale)
    args.sigma_scale  = float(args.sigma_scale)
    args.charge_scale = float(args.charge_scale)
    args.density      = float(args.density)
    args.N_density    = float(args.N_density)
    args.gens         = int(args.gens)
    if args.density != 0 and args.N_density != 0: print "ERROR: Setting both the -d and -d_N options is inconsistent (both specify density). Please revise. Exiting..."; quit()

    # Parse output filename 
    if args.outputname != '':

        # Handle the special case where the script is run from without the output directory
        if args.outputname == ".":
            print "ERROR: the output folder must differ from the current folder."
            quit()
            Filename = os.getcwd().split('/')[-1]
            args.coord_files = [ Filename+'/'+i if i[0] != '/' else i for i in args.coord_files ]  # Prepend the current folder to the path of relatively specified *.xyz files
            args.FF_files = [ Filename+'/'+i if i[0] != '/' else i for i in args.FF_files ]  # Prepend the current folder to the path of relatively specified *.xyz files
            os.chdir('..')

        # Else, directly assign the output directory based on the -o argument
        else:
            Filename = args.outputname        

    # If the output folder is not specified then the program defaults to the root of the first supplied *xyz file
    else:
        Filename = args.coord_files[0].split('.')[0]
    
    # Check that the input is an .xyz file. 
    if read_flag == 0:
        for i in args.coord_files:
            if i.split('.')[-1] != 'xyz':
                print "ERROR: Check to ensure that the input coordinate file(s) are in .xyz format."
                quit()
            elif os.path.isfile(i) != True:
                print "ERROR: Specified *.xyz file ({}) does not exit. Please check the path. Exiting...".format(i)
                quit()

    # Check that the supplied FF files exist
    for i in args.FF_files:
        if os.path.isfile(i) != True:
            print "ERROR: Specified FF_file ({}) does not exist. Please check the path. Exiting...".format(i)
            quit()

    # If the output directory already exists, check for md cycle information to avoid overwriting data
    
    if os.path.isdir(Filename):

        # Check if the specified folder contains md cycles
        os.chdir(Filename)
        md_folders = [d for dp, dn, filenames in os.walk('.') for d in dn if fnmatch.fnmatch(d,"md-*") ]        

        # Initialize subfolder based on the number of previous cycles present in this folder
        if len(md_folders) == 0:
            current_cycle = 1
        else:
            current_cycle = max([ int(i.split('-')[-1]) for i in md_folders ]) + 1
        Filename= 'md-'+str(current_cycle)
        os.makedirs(Filename)            
        args.coord_files = [ '../'+i if i[0] != '/' else i for i in args.coord_files ]  # Prepend an up-folder call to the path of relatively specified *.xyz files
        args.FF_files = [ '../'+i if i[0] != '/' else i for i in args.FF_files ]  # Prepend an up-folder call to the path of relatively specified *.xyz files
        sys.stdout = Logger(Filename)
        print "PROGRAM CALL: python gen_md_for_vdw.py {}\n".format(' '.join([ i for i in argv]))            

        # If the configuration is being read, then things are simplified since there is no geometry initialization and the
        # subsequent few lines of code are all that are required for generating the job. 
        if read_flag == 1:

            # Check for the existence of the necessary files from the previous run and update the vdw parameters
            cycle_read,Dihedral_styles = prepare_data("md-{}".format(current_cycle-1),"md-{}".format(current_cycle))
            Write_input(Filename,args.T_equil,args.t_equil,args.T_anneal,args.t_anneal,args.pressure,args.frequency,args.onefourscale,args.pair_styles,Dihedral_styles,\
                        args.density,args.N_density,[],args.fixes,[],[])

            # Print banner
            print "Generated a new molecular dynamics cycle ({}) using the last configuraiton from run md-{} and vdw parameters from cycle {}.".format(Filename,current_cycle-1,cycle_read)    
            quit()

        else:
            print "Generating a new molecular dynamics cycle ({}) using supplied *.db and *.xyz files".format(Filename)    

    # Else, create a bare directory
    else:
        os.makedirs(Filename)
        sys.stdout = Logger(Filename)
        print "PROGRAM CALL: python gen_md_for_vdw.py {}\n".format(' '.join([ i for i in argv]))

    # Catenate the FF files and copy them to the output folder
    FF_all = Filename+"/"+Filename+'.db'
    with open(Filename+"/"+Filename+'.db','w') as f:
        for i in args.FF_files:
            with open(i,'r') as ff_file:
                for j in ff_file:
                    f.write(j)
            f.write("\n")

    # Grab data on each molecule being added to the MD simulation.
    Data = get_data(FF_all,args.coord_files,args.N_mol,args.gens)

    # Pack the simulation cell
    Elements,Atom_types,Geometry,Bonds,Bond_types,Angles,Angle_types,Dihedrals,Dihedral_types,Charges,Molecule,Molecule_names,Adj_mat,Sim_Box = \
    Pack_Box(Data,Box_offset=0,Density=args.density,N_Density=args.N_density)

    # Generate VDW parameters    
    VDW_params = initialize_VDW(sorted(set(Atom_types)),sigma_scale=args.sigma_scale,eps_scale=args.eps_scale,VDW_type="lj/cut/coul/long",VDW_FF=Data[args.coord_files[0]]["VDW_params"],Force_UFF=args.force_UFF)

    # Generate Simulation Dictionaries
    # The bond, angle, and diehdral parameters for each molecule are combined into one dictionary
    Bond_params = {}; Angle_params = {}; Dihedral_params = {}; Masses = {}
    for i in Data.keys():
        for j in Data[i]["Bond_params"].keys(): Bond_params[j] = Data[i]["Bond_params"][j]
        for j in Data[i]["Angle_params"].keys(): Angle_params[j] = Data[i]["Angle_params"][j]
        for j in Data[i]["Dihedral_params"].keys(): Dihedral_params[j] = Data[i]["Dihedral_params"][j]
        for j in Data[i]["Masses"].keys(): Masses[j] = Data[i]["Masses"][j]
        
    # Write the lammps datafile (the function returns a dictionary, Atom_type_dict, holding the mapping between atom_types and the lammps id numbers; this mapping is needed for setting fixes)
    print "Writing LAMMPS datafile ({})...".format(Filename+'.data')
    print "Writing LAMMPS settings file ({})...".format(Filename+'.in.settings')
    Atom_type_dict = Write_data(Filename,Atom_types,Sim_Box,Elements,Geometry,Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,\
                                Dihedrals,Dihedral_types,Dihedral_params,Charges,VDW_params,Masses,Molecule)

    # Generate fix statement for scaled charges
    fixes = scale_charges(args.charge_scale,Atom_type_dict,Atom_types,Charges)

    # Gather the different dihedral_styles
    Dihedral_styles = set([ Dihedral_params[i][0] for i in Dihedral_types ])

    # Write the lammps input files
    print "Writing LAMMPS input file ({})...".format(Filename+'.in.init')
    Write_input(Filename,args.T_equil,args.t_equil,args.T_anneal,args.t_anneal,args.pressure,args.frequency,args.onefourscale,args.pair_styles,Dihedral_styles,args.density,args.N_density,fixes,args.fixes,Bond_types,Angle_types)

    # Write map file to more easily post-process the trajectory
    print "Writing mapfile ({})...".format(Filename+'.map')
    write_map(Filename,Elements,Atom_types,Charges,Masses,Adj_mat,zeros([len(Elements)]),args.N_mol)

    # Writing molecule map file
    print "Writing molfile ({})...".format(Filename+'.mol.txt')
    write_molecule(Filename,Molecule_names)

    # Write *.pairs file used by gen_jobs_for_vdw.py for determining the initial guess
    print "Writing pair file ({})...".format(Filename+'.pairs')
    Write_pairs(Filename,Atom_types,VDW_params,args.eps_scale,args.sigma_scale)
    
    # Remove concatenated FF file from run folder
    os.remove(FF_all)

    # Print banner
    print "\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Success! Have a Nice Day!","*"*167)

# Description: A wrapper for the write commands for generating the lammps input file
def Write_input(Filename,equil_temp,equil_time,anneal_temp,anneal_time,pressure,frequency,Onefourscale,Pair_styles,Dihedral_styles,Density,N_Density,fixes,ba_fixes,Bond_types,Angle_types):
    
    # Initialize rigid bond/angle fix commands
    ba_fix_cmd = ''
    if 'bonds' and 'angles' in ba_fixes:
        ba_fix_cmd="fix rigid all rattle 0.0001 20 ${coords_freq} "+"b {} a {}\n".format(" ".join([ str(count_i+1) for count_i,i in enumerate(set(Bond_types))])," ".join([ str(count_i+1) for count_i,i in enumerate(set(Angle_types))]))
    elif 'bonds' in ba_fixes:
        ba_fix_cmd="fix rigid all rattle 0.0001 20 ${coords_freq} "+"b {}\n".format(" ".join([ str(count_i+1) for count_i,i in enumerate(set(Bond_types))]))
    elif 'angles' in ba_fixes:
        ba_fix_cmd="fix rigid all rattle 0.0001 20 ${coords_freq} "+"a {}\n".format(" ".join([ str(count_i+1) for count_i,i in enumerate(set(Angle_types))]))

    with open(Filename+'/'+Filename+'.in.init','w') as f:
        f.write("# lammps input file for polymer simulation with dilute ions\n\n"+\
                "# VARIABLES\n"+\
                "variable        data_name       index   {}\n".format(Filename+'.data')+\
                "variable        settings_name   index   {}\n".format(Filename+'.in.settings')+\
                "variable        log_name        index   {}\n".format(Filename+'.log')+\
                "variable        nSteps_ramp     index   {} # number of data steps for the ramped anneal\n".format(anneal_time/2)+\
                "variable        nSteps_equil    index   {} # number of data steps for the ramped anneal\n".format(equil_time)+\
                "variable        avg_freq        index   1000\n".format(frequency)+\
                "variable        coords_freq     index   1000\n".format(frequency)+\
                "variable        thermo_freq     index   1000\n".format(frequency)+\
                "variable        dump4avg        index   100\n"+\
                "variable        vseed           index   {: <6d}\n".format(int(random.random()*100000))+\
                "variable        ANNEAL_TEMP     index   {} # Temperature during the initial anneal\n".format(anneal_temp)+\
                "variable        FINAL_TEMP      index   {} # Temperature ramped to during the final anneal\n".format(equil_temp)+\
                "variable        pressure        index   {} # Pressure during the simulations\n\n".format(pressure)+\

                "# Change the name of the log output #\n"+\
                "log ${log_name}\n\n"+\
                
                "#===========================================================\n"+\
                "# GENERAL PROCEDURES\n"+\
                "#===========================================================\n"+\
                "units		real	# g/mol, angstroms, fs, kcal/mol, K, atm, charge*angstrom\n"+\
                "dimension	3	# 3 dimensional simulation\n"+\
                "newton		off	# use Newton's 3rd law\n"+\
                "boundary	p p p	# periodic boundary conditions \n"+\
                "atom_style	full    # molecular + charge\n\n"+\

                "#===========================================================\n"+\
                "# FORCE FIELD DEFINITION\n"+\
                "#===========================================================\n"+\
                'special_bonds  lj   0.0 0.0 0.0     # NO     1-4 LJ interactions\n'+\
                'special_bonds  coul 0.0 0.0 {}     # REDUCE 1-4 electrostatics\n'.format(Onefourscale)+\
                'pair_style     hybrid {} # outer_LJ outer_Coul (cutoff values, see LAMMPS Doc)\n'.format(Pair_styles)+\
                'bond_style     harmonic             # parameters needed: k_bond, r0\n'+\
                'angle_style    harmonic             # parameters needed: k_theta, theta0\n')
        if len(Dihedral_styles) >= 1:
            f.write('dihedral_style hybrid {}            # parameters needed: V1, V2, V3, V4\n'.format(' '.join(Dihedral_styles)))            
        f.write('kspace_style   pppm 0.0001          # long-range electrostatics sum method\n'+\
                'pair_modify    shift yes mix arithmetic       # using Lorenz-Berthelot mixing rules\n\n'+\

                "#===========================================================\n"+\
                "# SETUP SIMULATIONS\n"+\
                "#===========================================================\n\n"+\
                "# READ IN COEFFICIENTS/COORDINATES/TOPOLOGY\n"+\
                'read_data ${data_name}\n'+\
                'include ${settings_name}\n\n'+\

                "# SET RUN PARAMETERS\n"+\
                "timestep	1.0		# fs\n"+\
                "run_style	verlet 		# Velocity-Verlet integrator\n"+\
                "neigh_modify every 1 delay 0 check no # More relaxed rebuild criteria can be used\n\n")
        
        if fixes != []:
            f.write("{}".format(fixes))

        f.write("#===========================================================\n"+\
                "# RUN CONSTRAINED RELAXATION\n"+\
                "#===========================================================\n\n"+

                "# SET OUTPUTS\n"+\
                "thermo_style    custom step temp vol density etotal pe ebond eangle edihed ecoul elong evdwl enthalpy press\n"+\
                "thermo_modify   format float %14.6f\n"+\
                "thermo ${thermo_freq}\n\n"+\

                "# DECLARE RELEVANT OUTPUT VARIABLES\n"+\
                "variable        my_step equal   step\n"+\
                "variable        my_temp equal   temp\n"+\
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

                "fix  averages all ave/time ${dump4avg} $(v_avg_freq/v_dump4avg) ${avg_freq} v_my_temp v_my_rho v_my_vol v_my_pe v_my_edih v_my_evdw v_my_eel v_my_ent v_my_P file thermo.avg\n\n"+\

                "# INITIALIZE VELOCITIES AND CREATE THE CONSTRAINED RELAXATION FIX\n"+\
                "# NOTE: MULTIPLE NVE CYCLES WITH REINITIALIZED VELOCITIES IS A PRECAUTION FOR HIGHLY DISTORTED INITIAL CONDITIONS\n"+\
                "velocity        all create ${ANNEAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n\n"+\

                "fix relax all nve/limit 0.1\n"+\
                "run             10000\n"+\

                "# INITIALIZE VELOCITIES AND CREATE THE CONSTRAINED RELAXATION FIX\n"+\
                "velocity        all create ${ANNEAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n\n"+\

                "fix relax all nve/limit 0.1\n"+\
                "run             10000\n"+\

                "# INITIALIZE VELOCITIES AND CREATE THE CONSTRAINED RELAXATION FIX\n"+\
                "velocity        all create ${ANNEAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n\n"+\

                "fix relax all nve/limit 0.1\n"+\
                "run             10000\n"+\

                "unfix relax\n\n"+\

                "#===========================================================\n"+\
                "# RUN RAMPED ANNEAL\n"+\
                "#===========================================================\n\n"+\

                "# REINITIALIZE THE VELOCITIES AND CREATE THE ANNEALING FIX\n"+\
                "velocity        all create ${ANNEAL_TEMP} ${vseed} mom yes rot yes     # DRAW VELOCITIES\n\n")

        if Density == 0.0 and N_Density == 0.0:
            f.write("fix anneal all npt temp ${ANNEAL_TEMP} ${FINAL_TEMP} 100.0 iso ${pressure} ${pressure} 1000.0 # NPT, nose-hoover 100 fs T relaxation\n\n")
        else:
            f.write("fix anneal all nvt temp ${ANNEAL_TEMP} ${FINAL_TEMP} 100.0 # NVT, nose-hoover 100 fs T relaxation\n\n")

        f.write("# CREATE COORDINATE DUMPS FOR ANNEAL\n"+\
                "dump anneal all custom ${coords_freq} anneal.lammpstrj id type x y z mol\n"+\
                "dump_modify anneal sort  id\n\n"+\

                "# RUN RAMPED ANNEAL\n")
        f.write("run	        ${nSteps_ramp}\n")
        f.write("unfix anneal\n\n"+\

                "# RUN EQUILIBRATION PHASE AT FINAL TEMP\n")
        if Density == 0.0 and N_Density == 0.0:
            f.write("fix anneal all npt temp ${FINAL_TEMP} ${FINAL_TEMP} 100.0 iso ${pressure} ${pressure} 1000.0 # NPT, nose-hoover 100 fs T relaxation\n")
        else:
            f.write("fix anneal all nvt temp ${FINAL_TEMP} ${FINAL_TEMP} 100.0 # NVT, nose-hoover 100 fs T relaxation\n")

        if len(ba_fix_cmd) != 0:
            f.write(ba_fix_cmd)
        f.write("run ${nSteps_ramp}\n")
        if len(ba_fix_cmd) != 0:
            f.write("unfix rigid\n")
        f.write("unfix anneal\n"+\
                "undump anneal\n\n"+\

                "#===========================================================\n"+\
                "# RUN EQUILIBRIUM SIM\n"+\
                "#===========================================================\n\n"+

                "# UPDATE RUN PARAMETERS AND CREATE FIX\n")
        if Density == 0.0 and N_Density == 0.0:
            f.write("fix equil all npt temp ${FINAL_TEMP} ${FINAL_TEMP} 100.0 iso ${pressure} ${pressure} 1000.0 # NPT, nose-hoover 100 fs T relaxation\n\n")
        else:
            f.write("fix equil all nvt temp ${FINAL_TEMP} ${FINAL_TEMP} 100.0 # NVT, nose-hoover 100 fs T relaxation\n\n")

        f.write("# CREATE COORDINATE DUMPS FOR EQUILIBRIUM\n"+\
                "dump equil all custom ${coords_freq} equil.lammpstrj id type x y z mol\n"+\
                "dump_modify equil sort  id\n\n"+\

                "# RUN EQUIL\n")
        if len(ba_fix_cmd) != 0:
            f.write(ba_fix_cmd)
        f.write("run		${nSteps_equil}\n")
        if len(ba_fix_cmd) != 0:
            f.write("unfix rigid\n")
        f.write("unfix equil\n"+\
                "undump equil\n\n"+\
                    
                "# WRITE RESTART FILES, CLEANUP, AND EXIT\n"+\
                "write_restart   {}\n".format(Filename.split('/')[-1]+'.end.restart')+\
                "write_data      {} pair ii\n".format(Filename.split('/')[-1]+'.end.data')+\
                "unfix		averages\n")
    f.close()

    return

# Description: A wrapper for the commands to generate a cubic box and array of molecules for the lammps run
def Pack_Box(Data,Box_offset=0,Density=0.0,N_Density=0.0):

    
    print "\nPacking simulation cell..."
    # Center each molecule at the origin
    for i in Data.keys():
        Data[i]["Geometry"] -= (mean(Data[i]["Geometry"][:,0]),mean(Data[i]["Geometry"][:,1]),mean(Data[i]["Geometry"][:,2]))

    # Define box circumscribing the molecule
    for i in Data.keys():
        Data[i]["Mol_Box"] = \
        ( min(Data[i]["Geometry"][:,0]),max(Data[i]["Geometry"][:,0]),min(Data[i]["Geometry"][:,1]),max(Data[i]["Geometry"][:,1]),min(Data[i]["Geometry"][:,2]),max(Data[i]["Geometry"][:,2]) )

    # Find the largest step_size
    # Use the geometric norm of the circumscribing cube as the step size for tiling
    Step_size=0.0
    for i in Data.keys():
        Current_step = ( (Data[i]["Mol_Box"][1]-Data[i]["Mol_Box"][0])**2 +\
                         (Data[i]["Mol_Box"][3]-Data[i]["Mol_Box"][2])**2 +\
                         (Data[i]["Mol_Box"][5]-Data[i]["Mol_Box"][4])**2 )**(0.5)+3  # +3 is just to be safe in the case of perfect alignment
        if Step_size < Current_step:
            Step_size = Current_step


    # Find the smallest N^3 cubic lattice that is greater than N_tot
    N_tot = sum([ Data[i]["N_mol"] for i in Data.keys()]) # Define the total number of molecules to be placed
    N_lat = 1
    while(N_lat**3 < N_tot):
        N_lat = N_lat + 1

    # Find molecular centers
    Centers = zeros([N_tot,3])
    count = 0
    for i in range(N_lat):

        if count == N_tot:
            break

        for j in range(N_lat):

            if count == N_tot:
                    break

            for k in range(N_lat):

                if count == N_tot:
                    break
                
                Centers[i*N_lat**2 + j*N_lat + k] = array([Step_size*i,Step_size*j,Step_size*k])
                count = count + 1    
    
    # Randomize the order of the centers so that molecules are randomly placed
    random.shuffle(Centers[:])

    # Find the box extrema. The same value is used for x y and z so that a cubic box is obtained (Note, the extra "Step_size" builds in padding even without Box_offset)
    low = min([min(Centers[:,0])-Box_offset-Step_size,min(Centers[:,1])-Box_offset-Step_size,min(Centers[:,2])-Box_offset-Step_size])
    high = max([max(Centers[:,0])+Box_offset+Step_size,max(Centers[:,1])+Box_offset+Step_size,max(Centers[:,2])+Box_offset+Step_size])

    # Center the box at the origin
    disp = -1*mean([low,high])
    low += disp
    high += disp
    for count_i,i in enumerate(Centers):
        Centers[count_i] = i + disp

    # Define sim box
    Sim_Box = array([low,high,low,high,low,high])

    # Intialize lists for iterating over molecules and keeping track of how many have been placed
    keys = Data.keys()             # list of unique molecule keys
    placed_num = [0]*len(keys)     # list for keeping track of how many of each molecule have been placed
    atom_index = 0                 # an index for keeping track of how many atoms have been placed
    mol_index = 0                  # an index for keeping track of how many molecules have been placed

    # Initialize various lists to hold the elements, atomtype, molid labels etc.
    # Create a list of molecule ids for each atom        
    Geometry_sim       = zeros([sum([ Data[i]["N_mol"]*len(Data[i]["Geometry"]) for i in Data.keys() ]),3])
    Adj_mat_sim        = zeros([sum([ Data[i]["N_mol"]*len(Data[i]["Geometry"]) for i in Data.keys() ]),sum([ Data[i]["N_mol"]*len(Data[i]["Geometry"]) for i in Data.keys() ])])
    Molecule_sim       = []
    Molecule_files     = []
    Elements_sim       = [] 
    Atom_types_sim     = []
    Bonds_sim          = []
    Bond_types_sim     = []
    Angles_sim         = []
    Angle_types_sim    = []
    Dihedrals_sim      = []
    Dihedral_types_sim = []
    Charges_sim        = []
    Masses_sim         = []


    # Place molecules in the simulation box and extend simulation lists 
    while (sum(placed_num) < N_tot):

        # Place the molecules in round-robin fashion to promote mixing
        # Each molecule key is iterated over, if all molecules of this type have been
        # placed then the molecule type is skipped
        for count_k,k in enumerate(keys):
            
            # If all the current molecules types have been placed continue
            # else: increment counter
            if placed_num[count_k] >= Data[k]["N_mol"]:
                continue
            else:
                placed_num[count_k]+=1
            print "placing atoms {}:{} with molecule {}".format(atom_index,len(Data[k]["Geometry"]),k)

            # perform x rotations
            angle = random.random()*360
            for count_j,j in enumerate(Data[k]["Geometry"]):
                Data[k]["Geometry"][count_j,:] = axis_rot(j,array([1.0,0.0,0.0]),array([0.0,0.0,0.0]),angle,mode='angle')

            # perform y rotations
            angle = random.random()*360
            for count_j,j in enumerate(Data[k]["Geometry"]):
                Data[k]["Geometry"][count_j,:] = axis_rot(j,array([0.0,1.0,0.0]),array([0.0,0.0,0.0]),angle,mode='angle')

            # perform z rotations
            angle = random.random()*360
            for count_j,j in enumerate(Data[k]["Geometry"]):
                Data[k]["Geometry"][count_j,:] = axis_rot(j,array([0.0,0.0,1.0]),array([0.0,0.0,0.0]),angle,mode='angle')

            # Move the current molecule to its box, append to Geometry_sim, return to the molecule to the origin
            Data[k]["Geometry"] += Centers[mol_index]
            Geometry_sim[atom_index:(atom_index+len(Data[k]["Geometry"])),:] = Data[k]["Geometry"]
            Data[k]["Geometry"] -= Centers[mol_index]

            # Extend various lists (total elements lists, atomtypes lists, etc)
            # Note: the lammps input expects bonds,angles,dihedrals, etc to be defined in terms of atom
            #       id so, the atom_index is employed to keep track of how many atoms have been placed.
            Adj_mat_sim[atom_index:(atom_index+len(Data[k]["Geometry"])),atom_index:(atom_index+len(Data[k]["Geometry"]))] = Data[k]["Adj_mat"]
            Molecule_sim       = Molecule_sim + [mol_index]*len(Data[k]["Elements"])
            Molecule_files     = Molecule_files + [k]
            Elements_sim       = Elements_sim + Data[k]["Elements"]
            Atom_types_sim     = Atom_types_sim + Data[k]["Atom_types"]
            Bonds_sim          = Bonds_sim + [ (j[0]+atom_index,j[1]+atom_index) for j in Data[k]["Bonds"] ]
            Bond_types_sim     = Bond_types_sim + Data[k]["Bond_types"]
            Angles_sim         = Angles_sim + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index) for j in Data[k]["Angles"] ]
            Angle_types_sim    = Angle_types_sim + Data[k]["Angle_types"]
            Dihedrals_sim      = Dihedrals_sim + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index,j[3]+atom_index) for j in Data[k]["Dihedrals"] ]
            Dihedral_types_sim = Dihedral_types_sim + Data[k]["Dihedral_types"]
            Charges_sim        = Charges_sim + Data[k]["Charges"]
            Masses_sim         = Masses_sim + [ Data[k]["Masses"][j] for j in Data[k]["Atom_types"] ] 

            # Increment atom_index based on the number of atoms in the current geometry
            atom_index += len(Data[k]["Geometry"])
            mol_index += 1

    # If a mass density is set, then the coordinates of the molecules and simulation box are rescaled to match the requested density.
    if Density != 0.0:
        print "\nRescaling box size and coordinates to match requested mass density ({} g/cc):\n".format(Density)
        A = 6.0221413e23
        mass_in_box = sum(Masses_sim) / A                                 # sum atomic masses then convert to g
        box_vol = ((Sim_Box[1]-Sim_Box[0])*10.0**(-8))**3                 # Standard density units are per cm^3 in lammps
        current_density = mass_in_box/box_vol                             # save current density to variable
        rescale_factor = (current_density/Density)**(1./3.)               # calculate rescale_factor based on the density ratio. cubed root is owing to the 3-dimensionality of the sim box.
        print "\tmass_in_box:     {:< 12.6f} g".format(mass_in_box)
        print "\tbox_vol:         {:< 12.6f} cm^3".format(box_vol)
        print "\tcurrent_density: {:< 12.6f} g/cm^3".format(current_density)
        print "\trescale_factor:  {:< 12.6f} ".format(rescale_factor)

        # If the requested density is less than the current density, then only the box is expanded
        if rescale_factor > 1:
            Sim_Box *= rescale_factor
        # For a contraction, rescale coordinates and box to ensure there are no clashes.
        else:
            Sim_Box *= rescale_factor
            Geometry_sim *= rescale_factor
        print "\nFixing the mass density by rescaling the coordinates by {}\n".format(rescale_factor)

    # If a number density is set, then the coordinates of the molecules and simulation box are rescaled to match the requested density.
    if N_Density != 0.0:
        print "\nRescaling box size and coordinates to match requested number density ({} atoms/A^3):\n".format(N_Density)
        A = 6.0221413e23
        atoms_in_box = len(Elements_sim)                                  # sum atomic masses then convert to g
        box_vol = (Sim_Box[1]-Sim_Box[0])**3                              # number density is specified in 1/A^3
        current_N_density = atoms_in_box/box_vol                          # save current density to variable
        rescale_factor = (current_N_density/N_Density)**(1./3.)           # calculate rescale_factor based on the density ratio. cubed root is owing to the 3-dimensionality of the sim box.
        print "\tatoms_in_box:      {:< 16d} atoms".format(atoms_in_box)
        print "\tbox_vol:           {:< 16.6f} A^3".format(box_vol)
        print "\tcurrent_N_density: {:< 16.6f} atoms/A^3".format(current_N_density)
        print "\trescale_factor:    {:< 16.6f} ".format(rescale_factor)

        # If the requested density is less than the current density, then only the box is expanded
        if rescale_factor > 1:
            Sim_Box *= rescale_factor
        # For a contraction, rescale coordinates and box to ensure there are no clashes.
        else:
            Sim_Box *= rescale_factor
            Geometry_sim *= rescale_factor
        print "\nFixing the number density by rescaling the coordinates by {}\n".format(rescale_factor)
        
    print "The simulation cell has {} molecules and dimensions of {:<3.2f} x {:<3.2f} x {:<3.2f}".format(N_tot,Sim_Box[1]-Sim_Box[0],Sim_Box[3]-Sim_Box[2],Sim_Box[5]-Sim_Box[4])
    return Elements_sim,Atom_types_sim,Geometry_sim,Bonds_sim,Bond_types_sim,Angles_sim,Angle_types_sim,Dihedrals_sim,Dihedral_types_sim,Charges_sim,Molecule_sim,Molecule_files,Adj_mat_sim,Sim_Box
                            
    f.close()

# Wrapper function for writing the lammps data file which is read by the .in.init file during run initialization
def Write_data(Filename,Atom_types,Sim_Box,Elements,Geometry,Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,Dihedrals,Dihedral_types,Dihedral_params,Charges,VDW_params,Masses,Molecule):

    # Write an xyz for easy viewing
    with open(Filename+'/'+Filename+'.xyz','w') as f:
        f.write('{}\n\n'.format(len(Geometry)))
        for count_i,i in enumerate(Geometry):
            f.write('{:20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(Elements[count_i],i[0],i[1],i[2]))

    # Create type dictionaries (needed to convert each atom,bond,angle, and dihedral type to consecutive numbers as per LAMMPS convention)
    # Note: LAMMPS orders atomtypes, bonds, angles, dihedrals, etc as integer types. Each of the following dictionaries holds the mapping between
    #       the true type (held in the various types lists) and the lammps type_id, obtained by enumerated iteration over the respective set(types).
    Atom_type_dict = {}
    for count_i,i in enumerate(sorted(set(Atom_types))):
        for j in Atom_types:
            if i == j:
                Atom_type_dict[i]=count_i+1
            if i in Atom_type_dict.keys():
                break
    Bond_type_dict = {}
    for count_i,i in enumerate(sorted(set(Bond_types))):
        for j in Bond_types:
            if i == j:
                Bond_type_dict[i]=count_i+1
            if i in Bond_type_dict.keys():
                break
    Angle_type_dict = {}
    for count_i,i in enumerate(sorted(set(Angle_types))):
        for j in Angle_types:
            if i == j:
                Angle_type_dict[i]=count_i+1
            if i in Angle_type_dict.keys():
                break
    Dihedral_type_dict = {}
    for count_i,i in enumerate(sorted(set(Dihedral_types))):
        for j in Dihedral_types:
            if i == j:
                Dihedral_type_dict[i]=count_i+1
            if i in Dihedral_type_dict.keys():
                break

    # Write the data file
    with open(Filename+'/'+Filename+'.data','w') as f:
        
        # Write system properties
        f.write("LAMMPS data file via vdw_self_gen.py, on {}\n\n".format(datetime.datetime.now()))

        f.write("{} atoms\n".format(len(Elements)))
        f.write("{} atom types\n".format(len(set(Atom_types))))
        if len(Bonds) > 0:
            f.write("{} bonds\n".format(len(Bonds)))
            f.write("{} bond types\n".format(len(set(Bond_types))))
        if len(Angles) > 0:
            f.write("{} angles\n".format(len(Angles)))
            f.write("{} angle types\n".format(len(set(Angle_types))))
        if len(Dihedrals) > 0:
            f.write("{} dihedrals\n".format(len(Dihedrals)))
            f.write("{} dihedral types\n\n".format(len(set(Dihedral_types))))

        # Write box dimensions
        f.write("{:< 20.16f} {:< 20.16f} xlo xhi\n".format(Sim_Box[0],Sim_Box[1]))
        f.write("{:< 20.16f} {:< 20.16f} ylo yhi\n".format(Sim_Box[2],Sim_Box[3]))
        f.write("{:< 20.16f} {:< 20.16f} zlo zhi\n\n".format(Sim_Box[4],Sim_Box[5]))

        # Write Masses
        f.write("Masses\n\n")
        for count_i,i in enumerate(sorted(set(Atom_types))):
            for j in set(Atom_types):
                if Atom_type_dict[j] == count_i+1:
                    f.write("{} {:< 8.6f}\n".format(count_i+1,Masses[str(j)])) # count_i+1 bc of LAMMPS 1-indexing
        f.write("\n")

        # Write Bond Coeffs
        f.write("Bond Coeffs\n\n")
        for count_i,i in enumerate(set(Bond_types)):
            for j in set(Bond_types):
                if Bond_type_dict[j] == count_i+1:
                    f.write("{} {:< 12.6f} {:< 12.6f}\n".format(count_i+1,Bond_params[j][0],Bond_params[j][1])) # count_i+1 bc of LAMMPS 1-indexing
        f.write("\n")

        # Write Angle Coeffs
        f.write("Angle Coeffs\n\n")
        for count_i,i in enumerate(set(Angle_types)):
            for j in set(Angle_types):
                if Angle_type_dict[j] == count_i+1:
                    f.write("{} {:< 12.6f} {:< 12.6f}\n".format(count_i+1,Angle_params[j][0],Angle_params[j][1])) # count_i+1 bc of LAMMPS 1-indexing
        f.write("\n")

        # Write Atoms
        f.write("Atoms\n\n")
        for count_i,i in enumerate(Atom_types):
            f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
            .format(count_i+1,Molecule[count_i],Atom_type_dict[i],Charges[count_i],Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],0,0,0))

        # Write Bonds
        if len(Bonds) > 0:
            f.write("\nBonds\n\n")
            for count_i,i in enumerate(Bonds):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Bond_type_dict[Bond_types[count_i]],i[0]+1,i[1]+1))

        # Write Angles
        if len(Angles) > 0:
            f.write("\nAngles\n\n")
            for count_i,i in enumerate(Angles):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Angle_type_dict[Angle_types[count_i]],i[0]+1,i[1]+1,i[2]+1))

        # Write Dihedrals
        if len(Dihedrals) > 0: 
            f.write("\nDihedrals\n\n")
            for count_i,i in enumerate(Dihedrals):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Dihedral_type_dict[Dihedral_types[count_i]],i[0]+1,i[1]+1,i[2]+1,i[3]+1))

    # Write the settings file
    with open(Filename+'/'+Filename+'.in.settings','w') as f:

        # Write non-bonded interactions (the complicated form of this loop is owed
        # to desire to form a nicely sorted list in terms of the lammps atom types
        # Note: Atom_type_dict was initialize according to sorted(set(Atom_types) 
        #       so iterating over this list (twice) orders the pairs, too.
        f.write("     {}\n".format("# Non-bonded interactions (pair-wise)"))
        for count_i,i in enumerate(sorted(set(Atom_types))):     
            for count_j,j in enumerate(sorted(set(Atom_types))): 

                # Skip duplicates
                if count_j < count_i:
                    continue

                # Conform to LAMMPS i <= j formatting
                if Atom_type_dict[i] <= Atom_type_dict[j]:
                    f.write("     {:20s} {:<10d} {:<10d} ".format("pair_coeff",Atom_type_dict[i],Atom_type_dict[j]))
                else:
                    f.write("     {:20s} {:<10d} {:<10d} ".format("pair_coeff",Atom_type_dict[j],Atom_type_dict[i]))

                # Determine key (ordered by initialize_VDW function such that i > j)
                if i > j:
                    key = (i,j)
                else:
                    key = (j,i)

                # Write the parameters
                for k in VDW_params[key]:
                    if type(k) is str:
                        f.write("{:20s} ".format(k))
                    if type(k) is float:
                        f.write("{:< 20.6f} ".format(k))
                f.write("\n")                        

        # Write stretching interactions
        # Note: Bond_type_dict was initialized by looping over sorted(set(Bond_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Stretching interactions"))
        for i in sorted(set(Bond_types)):
            f.write("     {:20s} {:<10d} ".format("bond_coeff",Bond_type_dict[i]))
            for j in Bond_params[i]:
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
            f.write("\n")

        # Write bending interactions
        # Note: Angle_type_dict was initialized by looping over sorted(set(Angle_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Bending interactions"))
        for i in sorted(set(Angle_types)):
            f.write("     {:20s} {:<10d} ".format("angle_coeff",Angle_type_dict[i]))
            for j in Angle_params[i]:
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
            f.write("\n")

        # Write dihedral interactions
        # Note: Dihedral_type_dict was initialized by looping over sorted(set(Dihedral_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Dihedral interactions"))
        for i in sorted(set(Dihedral_types)):
            f.write("     {:20s} {:<10d} ".format("dihedral_coeff",Dihedral_type_dict[i]))
            for j in Dihedral_params[i]:
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
                if type(j) is int:
                    f.write("{:< 20d} ".format(j))
            f.write("\n")

    return Atom_type_dict

# A wrapper for the commands to parse the bonds, angles, and dihedrals from the adjacency matrix.
# Returns:   list of atomtypes, bond_types, bond instances, angle_types, angle instances, dihedral_types,
#            diehdral instances, charges, and VDW parameters.
def Find_parameters(Adj_mat,Geometry,Atom_types,FF_db="FF_file"):

    # Initialize lists of each instance and type of FF object.
    # instances are stored as tuples of the atoms involved 
    # (e.g., bonds between atoms 1 and 13 and 17 and 5 would be stored as [(1,13),(17,5)] 
    # Similarly, types are stored as tuples of atom types.
    Bonds = []
    Bond_types = []
    Angles = []
    Angle_types = []
    Dihedrals = []
    Dihedral_types = []

    # Find bonds #
    print "Parsing bonds..."
    for count_i,i in enumerate(Adj_mat):        
        Tmp_Bonds = [ (count_i,count_j) for count_j,j in enumerate(i) if j == 1 and count_j > count_i ]

        # Store bond tuple so that lowest atom *type* between the first and the second atom is placed first
        # and avoid redundant placements
        for j in Tmp_Bonds:
            if Atom_types[j[1]] < Atom_types[j[0]] and (j[1],j[0]) not in Bonds and (j[0],j[1]) not in Bonds:
                Bonds = Bonds + [ (j[1],j[0]) ]
                Bond_types = Bond_types + [ (Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1]) not in Bonds and (j[1],j[0]) not in Bonds:
                Bonds = Bonds + [ (j[0],j[1]) ]
                Bond_types = Bond_types + [ (Atom_types[j[0]],Atom_types[j[1]]) ]


    # Find angles #
    print "Parsing angles..."
    for i in Bonds:        

        # Find angles based on connections to first index of Bonds
        Tmp_Angles = [ (count_j,i[0],i[1]) for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j != i[1] ]

        # Store angle tuple so that lowest atom *type* between the first and the third is placed first
        # and avoid redundant placements
        for j in Tmp_Angles:
            if Atom_types[j[2]] < Atom_types[j[0]] and (j[2],j[1],j[0]) not in Angles and (j[0],j[1],j[2]) not in Angles:
                Angles = Angles + [(j[2],j[1],j[0])]
                Angle_types = Angle_types + [ (Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2]) not in Angles and (j[2],j[1],j[0]) not in Angles:
                Angles = Angles + [(j[0],j[1],j[2])]
                Angle_types = Angle_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]]) ]

        # Find angles based on connections to second index of Bonds
        Tmp_Angles = [ (i[0],i[1],count_j) for count_j,j in enumerate(Adj_mat[i[1]]) if j == 1 and count_j != i[0] ]

        # Store angle tuple so that lowest atom *type* between the first and the third is placed first
        # and avoid redundant placements
        for j in Tmp_Angles:
            if Atom_types[j[2]] < Atom_types[j[0]] and (j[2],j[1],j[0]) not in Angles and (j[0],j[1],j[2]) not in Angles:
                Angles = Angles + [(j[2],j[1],j[0])]
                Angle_types = Angle_types + [ (Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2]) not in Angles and (j[2],j[1],j[0]) not in Angles:
                Angles = Angles + [(j[0],j[1],j[2])]
                Angle_types = Angle_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]]) ]

        
    # Find dihedrals #
    print "Parsing dihedrals..."
    for i in Angles:
        
        # Find atoms attached to first atom of each angle
        Tmp_Dihedrals = [ (count_j,i[0],i[1],i[2]) for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in [i[1],i[2]] ]
        
        # Store dihedral tuple so that the lowest atom *type* between the first and fourth is placed first
        # and avoid redundant placements        
        for j in Tmp_Dihedrals:

            # If the first and fourth atoms are equal, then sorting is based on the second and third
            if Atom_types[j[3]] == Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                if Atom_types[j[2]] < Atom_types[j[1]]:
                    Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
                else:
                    Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

            elif Atom_types[j[3]] < Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2],j[3]) not in Dihedrals and (j[3],j[2],j[1],j[0]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

        # Find atoms attached to the third atom of each angle
        Tmp_Dihedrals = [ (i[0],i[1],i[2],count_j) for count_j,j in enumerate(Adj_mat[i[2]]) if j == 1 and count_j not in [i[0],i[1]] ]
        
        # Store dihedral tuple so that the lowest atom *type* between the first and fourth is placed first
        # and avoid redundant placements        
        for j in Tmp_Dihedrals:

            # If the first and fourth atoms are equal, then sorting is based on the second and third
            if Atom_types[j[3]] == Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                if Atom_types[j[2]] < Atom_types[j[1]]:
                    Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
                else:
                    Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                    Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

            elif Atom_types[j[3]] < Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
                Dihedral_types = Dihedral_types+ [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
            elif (j[0],j[1],j[2],j[3]) not in Dihedrals and (j[3],j[2],j[1],j[0]) not in Dihedrals:
                Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
                Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]
                    
    ##############################################################
    # Read in parameters: Here the stretching, bending, dihedral #
    # and non-bonding interaction parameters are read in from    #
    # parameters file. Mass and charge data is also included.    #
    # The program looks for a simple match for the first entry   #
    # in each line with one of the bond or angle types.          #
    # INPUT: param_file, BOND_TYPES_LIST, ANGLE_TYPES_LIST       #
    #        DIHEDRAL_TYPES_LIST, ELEMENTS                       #
    # OUTPUT: CHARGES, MASSES, BOND_PARAMS, ANGLE_PARAMS,        #
    #         DIHERAL_PARAMS, PW_PARAMS                          #
    ##############################################################

    # Initialize dictionaries

    # Read in masses and charges
    Masses = {}
    with open(FF_db,'r') as f:
        content=f.readlines()
        
    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() == 'atom':
            Masses[fields[1]] = float(fields[3]) 

    # Read in bond parameters
    Bond_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0] == 'bond':
            if fields[3] == "harmonic":
                Bond_params[(fields[1],fields[2])] = [float(fields[4]),float(fields[5])]
            else:
                print "ERROR: only harmonic bond definitions are currently supported by gen_md_for_vdw.py. Exiting..."
                quit()

    # Read in angle parameters
    Angle_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() == 'angle':
            if fields[4] == "harmonic":
                Angle_params[(fields[1],fields[2],fields[3])] = [float(fields[5]),float(fields[6])]
            else:
                print "ERROR: only harmonic angle definitions are currently supported by gen_md_for_vdw.py. Exiting..."
                quit()

    # Read in dihedral parameters
    Dihedral_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() in ['dihedral','torsion']:            
            if fields[5] == "opls":
                Dihedral_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
            elif fields[5] == "harmonic":
#                Dihedral_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(1)]
                Dihedral_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ]
            else:
                print "ERROR: Only opls type dihedral definitions are currently supported by gen_md_for_vdw.py. Exiting..."
                quit()
                

    # Search for charges based on atom type
    with open(FF_db,'r') as f:
        content=f.readlines()

    Charges = zeros(len(Atom_types))
    for i in range(len(Atom_types)):
        for lines in content:
            fields=lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['charge'] and Atom_types[i] == fields[1]:
                Charges[i] = float(fields[2])

    # Search for VDW parameters
    VDW_params = {}
    with open(FF_db,'r') as f:
        for lines in f:
            fields = lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['vdw']:

                if fields[1] > fields[2]:
                    VDW_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                else:
                    VDW_params[(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]

    # Check for missing parameters
    Missing_masses = [ i for i in Atom_types if str(i) not in Masses.keys() ] 
    Missing_charges = [ count_i for count_i,i in enumerate(Charges) if i == -100.0 ]; Missing_charges = [ Atom_types[i] for i in Missing_charges ]
    Missing_bonds = [ i for i in Bond_types if (i[0],i[1]) not in Bond_params.keys() ]
    Missing_angles = [ i for i in Angle_types if (i[0],i[1],i[2]) not in Angle_params.keys() ]
    Missing_dihedrals = [ i for i in Dihedral_types if (i[0],i[1],i[2],i[3]) not in Dihedral_params.keys() ]

    # Print diagnostics on missing parameters and quit if the prerequisites are missing.
    if ( len(Missing_masses) + len(Missing_charges) + len(Missing_bonds) + len(Missing_angles) + len(Missing_dihedrals) ) > 0:
        print "\nUh Oh! There are missing FF parameters...\n"

        if Missing_masses:
            print "Missing masses for the following atom types: {}".format([ i for i in set(Missing_masses) ])
        if Missing_charges:
            print "Missing charges for the following atom types: {}".format([ i for i in set(Missing_charges) ])
        if Missing_bonds:
            print "Missing bond parameters for the following bond types: {}".format([ i for i in set(Missing_bonds) ])
        if Missing_angles:
            print "Missing angle parameters for the following angle types: {}".format([ i for i in set(Missing_angles) ])
        if Missing_dihedrals:
            print "Missing dihedral parameters for the following dihedral types: {}".format([ i for i in set(Missing_dihedrals) ])
        
        print "\nEnsure the specification of the missing parameters. Exiting..."
        quit()

    return Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,Dihedrals,Dihedral_types,Dihedral_params,Charges.tolist(),Masses,VDW_params

# # Function for converting geometry into a UA geometry
# def Sub_UA(Elements,Adj_mat,Atom_types,Geometry):

#     # List of united atom types. By default, Hydrogens are removed for atom types belonging to this list (Carbons). 
#     #UA_types = range(1,51)
    
#     UA_idx = [ count_i for count_i,i in enumerate(Elements) if i == "C" ]
#     Delete_list = []
#     for i in UA_idx:
#         Delete_list = Delete_list + [ count_j for count_j,j in enumerate(Adj_mat[i]) if j == 1 and Elements[count_j] == 'H' ]
        
#     # Perform deletions
#     Elements = [ Elements[m] for m in range(len(Elements)) if m not in Delete_list ]
#     Atom_types = [ Atom_types[m] for m in range(len(Atom_types)) if m not in Delete_list ]
#     Geometry = Geometry[[ m for m in range(len(Geometry)) if m not in Delete_list ],:]
#     Adj_mat = Adj_mat[[ m for m in range(len(Adj_mat[:,0])) if m not in Delete_list ],: ]
#     Adj_mat = Adj_mat[:,[ m for m in range(len(Adj_mat[0,:])) if m not in Delete_list ]]                

#     # Center new geometry at origin
#     Centroid = array([mean(Geometry[:,0]),mean(Geometry[:,1]),mean(Geometry[:,2])])
#     for count_i,i in enumerate(Geometry):
#         Geometry[count_i] = i-Centroid

#     return Elements,Adj_mat,Geometry,Atom_types

# # Generates the adjacency matrix 
# # OLD METHOD BASED ON HEURISTIC RADII
# def Table_generator(Elements,Geometry):
    
#     Radii = {'H' :1.4, 'He':1.5, 
#              'Li':2.0, 'Be':1.6,                                                                                                     'B' :1.7, 'C' :1.7, 'N' :1.7, 'O' :1.6, 'F' :1.7, 'Ne':1.6,
#              'Na':2.2, 'Mg':1.9,                                                                                                     'Al':2.0, 'Si':2.1, 'P' :2.1, 'S' :2.0, 'Cl':1.9, 'Ar':1.7,
#              'K' :2.6, 'Ca':2.2, 'Sc':2.4, 'Ti':2.4, 'V' :2.6, 'Cr':2.5, 'Mn':2.6, 'Fe':2.5, 'Co':2.4, 'Ni':2.4, 'Cu':2.5, 'Zn':2.5, 'Ga':2.1, 'Ge':2.1, 'As':2.2, 'Se':2.1, 'Br':2.1, 'Kr':1.8,
#              'Rb':2.9, 'Sr':2.5, 'Y' :2.7, 'Zr':2.6, 'Nb':2.7, 'Mo':2.6, 'Tc':2.5, 'Ru':2.6, 'Rh':2.5, 'Pd':2.5, 'Ag':2.6, 'Cd':2.7, 'In':2.4, 'Sn':2.4, 'Sb':2.4, 'Te':2.4, 'I' :2.3, 'Xe':1.9,
#              'Cs':3.2, 'Ba':2.7, 'La':2.6, 'Hf':2.6, 'Ta':2.7, 'W' :2.6, 'Re':2.6, 'Os':2.6, 'Ir':2.6, 'Pt':2.6, 'Au':2.5, 'Hg':2.6, 'Tl':2.4, 'Pb':2.5, 'Bi':2.5, 'Po':2.5, 'At':2.5, 'Rn':2.0,
#              'default' : 1.8 }

#     # Print warning for uncoded elements.
#     for i in Elements:
#         if i not in Radii.keys():
#             print "ERROR: The geometry contains an element ({}) that the Table_generator function doesn't have bonding information for. This needs to be directly added to the Radii".format(i)+\
#                   " dictionary before proceeding. Exiting..."
#             quit()


#     # Generate distance matrix holding atom-atom separations (only save upper right)
#     Dist_Mat = triu(cdist(Geometry,Geometry))
    
#     # Find plausible connections
#     x_ind,y_ind = where( (Dist_Mat > 0.0) & (Dist_Mat < max([ Radii[i] for i in Radii.keys() ])) )

#     # Initialize Adjacency Matrix
#     Adj_mat = zeros([len(Geometry),len(Geometry)])

#     # Iterate over plausible connections and determine actual connections
#     for count,i in enumerate(x_ind):
        
#         # Determine connection based on the longer of the two radii
#         if Dist_Mat[i,y_ind[count]] < max(Radii[Elements[i]],Radii[Elements[y_ind[count]]]):
#             Adj_mat[i,y_ind[count]]=1

#     # Hermitize Adj_mat
#     Adj_mat=Adj_mat + Adj_mat.transpose()

#     return Adj_mat

# Generates the adjacency matrix based on UFF bond radii
def Table_generator(Elements,Geometry):

    # Initialize UFF bond radii (Rappe et al. JACS 1992)
    # NOTE: Units of angstroms 
    # NOTE: These radii neglect the bond-order and electronegativity corrections in the original paper. Where several values exist for the same atom, the largest was used. 
    Radii = {  'H':0.354, 'He':0.849,\
              'Li':1.336, 'Be':1.074,                                                                                                                          'B':0.838,  'C':0.757,  'N':0.700,  'O':0.658,  'F':0.668, 'Ne':0.920,\
              'Na':1.539, 'Mg':1.421,                                                                                                                         'Al':1.244, 'Si':1.117,  'P':1.117,  'S':1.064, 'Cl':1.044, 'Ar':1.032,\
               'K':1.953, 'Ca':1.761, 'Sc':1.513, 'Ti':1.412,  'V':1.402, 'Cr':1.345, 'Mn':1.382, 'Fe':1.335, 'Co':1.241, 'Ni':1.164, 'Cu':1.302, 'Zn':1.193, 'Ga':1.260, 'Ge':1.197, 'As':1.211, 'Se':1.190, 'Br':1.192, 'Kr':1.147,\
              'Rb':2.260, 'Sr':2.052,  'Y':1.698, 'Zr':1.564, 'Nb':1.473, 'Mo':1.484, 'Tc':1.322, 'Ru':1.478, 'Rh':1.332, 'Pd':1.338, 'Ag':1.386, 'Cd':1.403, 'In':1.459, 'Sn':1.398, 'Sb':1.407, 'Te':1.386,  'I':1.382, 'Xe':1.267,\
              'Cs':2.570, 'Ba':2.277, 'La':1.943, 'Hf':1.611, 'Ta':1.511,  'W':1.526, 'Re':1.372, 'Os':1.372, 'Ir':1.371, 'Pt':1.364, 'Au':1.262, 'Hg':1.340, 'Tl':1.518, 'Pb':1.459, 'Bi':1.512, 'Po':1.500, 'At':1.545, 'Rn':1.42,\
              'default' : 0.7 }
    
    # Scale factor is used for determining the bonding threshold. 1.2 is a heuristic that give some lattitude in defining bonds since the UFF radii correspond to equilibrium lengths. 
    scale_factor = 1.2

    # Print warning for uncoded elements.
    for i in Elements:
        if i not in Radii.keys():
            print "ERROR: The geometry contains an element ({}) that the Table_generator function doesn't have bonding information for. This needs to be directly added to the Radii".format(i)+\
                  " dictionary before proceeding. Exiting..."
            quit()


    # Generate distance matrix holding atom-atom separations (only save upper right)
    Dist_Mat = triu(cdist(Geometry,Geometry))
    
    # Find plausible connections
    x_ind,y_ind = where( (Dist_Mat > 0.0) & (Dist_Mat < max([ Radii[i]**2.0 for i in Radii.keys() ])) )

    # Initialize Adjacency Matrix
    Adj_mat = zeros([len(Geometry),len(Geometry)])

    # Iterate over plausible connections and determine actual connections
    for count,i in enumerate(x_ind):
        
        # Assign connection if the ij separation is less than the UFF-sigma value times the scaling factor
        if Dist_Mat[i,y_ind[count]] < (Radii[Elements[i]]+Radii[Elements[y_ind[count]]])*scale_factor:            
            Adj_mat[i,y_ind[count]]=1

    # Hermitize Adj_mat
    Adj_mat=Adj_mat + Adj_mat.transpose()

    # Perform some simple checks on bonding to catch errors
    problem_dict = { i:0 for i in Radii.keys() }
    conditions = { "H":1, "C":4, "F":1, "Cl":1, "Br":1, "I":1, "O":2, "N":4, "B":4 }
    for count_i,i in enumerate(Adj_mat):

        if Elements[count_i] in conditions.keys():
            if sum(i) > conditions[Elements[count_i]]:
                problem_dict[Elements[count_i]] += 1

    # Print warning messages for obviously suspicious bonding motifs.
    if sum( [ problem_dict[i] for i in problem_dict.keys() ] ) > 0:
        print "Table Generation Warnings:"
        for i in sorted(problem_dict.keys()):
            if problem_dict[i] > 0:
                if i == "H": print "WARNING: {} hydrogen(s) have more than one bond.".format(problem_dict[i])
                if i == "C": print "WARNING: {} carbon(s) have more than four bonds.".format(problem_dict[i])
                if i == "F": print "WARNING: {} fluorine(s) have more than one bond.".format(problem_dict[i])
                if i == "Cl": print "WARNING: {} chlorine(s) have more than one bond.".format(problem_dict[i])
                if i == "Br": print "WARNING: {} bromine(s) have more than one bond.".format(problem_dict[i])
                if i == "I": print "WARNING: {} iodine(s) have more than one bond.".format(problem_dict[i])
                if i == "O": print "WARNING: {} oxygen(s) have more than two bonds.".format(problem_dict[i])
                if i == "N": print "WARNING: {} nitrogen(s) have more than four bonds.".format(problem_dict[i])
                if i == "B": print "WARNING: {} bromine(s) have more than four bonds.".format(problem_dict[i])
        print ""

    return Adj_mat
        
def xyz_parse(input):

    # Open file and read contents into variable
    with open(input,'r') as f:
            content=f.readlines()

    # Find number of atoms and initialize various matrices
    Atom_Number = int(content[0].split()[0])
    Elements = ['']*Atom_Number
    Geometry = zeros([Atom_Number,3])
    Atom_types = [0]*Atom_Number
    types_flag = 1

    # Iterate over the remainder of contents and read the
    # geometry and elements into variable. Note that the
    # first two lines are considered a header
    count=0
    for lines in content[2:]:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        # Write geometry containing lines to variable
        if len(fields) > 3:
            Elements[count]=fields[0]
            Geometry[count,:]=array([float(fields[1]),float(fields[2]),float(fields[3])])

            # If atom type information is stored in the fifth column it is saved to variable
            if len(fields) >= 5:
                Atom_types[count] = fields[4]
            else:
                types_flag = 0

            # If the nuclear charge is stored in the sixth column it is saved to variable
            if len(fields) >= 6:
                Charges[count]=float(fields[5])
            count = count + 1

    if types_flag == 1:
        return Elements,Geometry,Atom_types
    else:
        return Elements,Geometry,[]

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

# Wrapper function for the write commands for creating the *.map file
def write_map(Filename,Elements,Atom_types,Charges,Masses,Adj_mat,Structure,N_mol):

    # Open file for writing and write header (first two lines of the map file are header)    
    with open(Filename+'/'+Filename+'.map','w') as f:
        f.write('{} {}\n {:<50} {:<10} {:<10} {:<14}  {:<13} {}\n'.format(len(Atom_types),sum(N_mol),'Atom_type','Element','Structure','Mass','Charge','Adj_mat'))
        for count_i,i in enumerate(Atom_types):
            adj_mat_entry = (' ').join([ str(count_j) for count_j,j in enumerate(Adj_mat[count_i,:]) if j == 1 ])
            f.write(' {:<50} {:<10} {:< 9d} {:<14.6f} {:< 14.8f} {}\n'.format(i,Elements[count_i],int(Structure[count_i]),Masses[str(i)],Charges[count_i],adj_mat_entry))
        f.close()

# Wrapper function for the write commands for creating the *.map file
def write_molecule(Filename,Molecules):

    # Open file for writing and write header (first two lines of the map file are header)    
    with open(Filename+'/'+Filename+'.mol.txt','w') as f:
        f.write('{:<30s} {:<30s}\n'.format("Molecule_id","Molecule_Name"))
        for count_i,i in enumerate(Molecules):
            f.write("{:<30} {:<30}\n".format(count_i,i))

        # Write mol clusters
        f.write('\nLists for each molecule:\n')
        for i in sorted(set(Molecules)):
            f.write('{}: '.format(i))
            for count_j,j in enumerate(Molecules):
                if i == j:
                    f.write(' {}'.format(count_j))
            f.write('\n')

# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def initialize_VDW(atomtypes,sigma_scale=1.0,eps_scale=1.0,VDW_type="lj/cut/coul/long",VDW_FF={},Force_UFF=0):

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

    # NEW: Initialize VDW_dict first guess based on element types and Lorentz-Berthelot mixing rules
    # Order of operations: (1) If the parameters are in the supplied FF database then they are used as is
    # (2) If the self-terms are in the supplied FF datebase then they are used to generate the mixed
    # interactions (3) UFF parameters are used. 
    VDW_dict = {}
    for count_i,i in enumerate(atomtypes):
        for count_j,j in enumerate(atomtypes):
            if count_i < count_j:
                continue

            # Check for parameters in the database
            if (i,j) in VDW_FF and Force_UFF != 1:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
            elif (j,i) in VDW_FF and Force_UFF != 1:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1:
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(j,j)][1])**(0.5) * eps_scale
                sigma  = (VDW_FF[(i,i)][2]+VDW_FF[(j,j)][2])/2.0 * sigma_scale
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]

            # Last resort: Use UFF parameters.
            else:
                type_1 = int(i.split('[')[1].split(']')[0])
                type_2 = int(j.split('[')[1].split(']')[0])
                eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) * eps_scale
                sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0 * sigma_scale
                if i > j:
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]

#     # OLD: Initialize VDW_dict first guess based on element types and Lorentz-Berthelot mixing rules
#     VDW_dict = {}
#     for count_i,i in enumerate(atomtypes):
#         for count_j,j in enumerate(atomtypes):
#             if count_i < count_j:
#                 continue

#             type_1 = int(i.split('[')[1].split(']')[0])
#             type_2 = int(j.split('[')[1].split(']')[0])
#             eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) * eps_scale
#             sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0 * sigma_scale
#             if i > j:
#                 if (i,j) in VDW_FF: VDW_dict[(i,j)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
#                 else: VDW_dict[(i,j)] = [VDW_type,eps,sigma]
#             else:
#                 if (j,i) in VDW_FF: VDW_dict[(j,i)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]
#                 else: VDW_dict[(j,i)] = [VDW_type,eps,sigma]

    # Print summary
    print "\n{}".format("*"*167)
    print "* {:^163s} *".format("Initializing VDW parameters for the simulation (those with * were read from the FF file(s))")
    print "*{}*".format("-"*165)
    print "* {:<50s} {:<50s} {:<20s}  {:<18s} {:<18s}   *".format("Type","Type","VDW_type","eps (kcal/mol)","sigma (angstroms)")
    print "{}".format("*"*167)
    for j in VDW_dict.keys():
        if j in VDW_FF and Force_UFF != 1:
            print "  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f} *".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2])
        else:
            print "  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f}".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2])
    print ""


    return VDW_dict

# Description:   A wrapper for the commands to create a data dictionary holding the FF, mode, and geometry information for each molecule
# Inputs:  Filename:
# Returns: Data: a dictionary keyed to the list of coord_files
#                a keyed entry is itself a dictionary holding the geometry, 
#                elements, charges, bonds, bond_types, angles, angle_types,
#                dihedrals, dihedral_types and masses.
#                subdictionary keys (e.g., Data[*]["Geometry"])
#                Geometry: array of the molecule's geometry
#                Elements: list of the molecule's elements (indexed to the rows of "Geometry")
#                Atom_types: list of the molecule's atom_types (indexed to the rows of "Geometry")
#                Charges: list of the molecule's charges (indexed to the rows of "Geometry")
#                Bonds: a list of the molecule's bonds, each element being a tuple of the atom indices that are bonded
#                Angles: a list of the molecule's angles, each element being a tuple of the atom indices that form an angle mode
#                Dihedrals: a list of the molecule's dihedrals, each element being a tuple of the atom indices that form a dihedral
#                Bond_types: a list of the molecule's bond_types, indexed to Bonds
#                Angle_types: a list of the molecule's angle_types, indexed to Angles
#                Dihedral_types: a list of the molecule's dihedral_types, indexed to Dihedrals
#
def get_data(FF_all,coord_files,N_mol,gens):

    # Initialize dictionary to hold the FF, mode, and geometry information about each molecule
    Data = {}

    # Iterate over all molecules being simulated and collect their FF, mode, and geometry information.
    for count_i,i in enumerate(coord_files):

        # Initialize dictionary for this geometry and set number of molecules to place in the simulated system
        Data[i] = {}      
        Data[i]["N_mol"] = N_mol[count_i]

        # Extract Element list and Coord list from the file
        Data[i]["Elements"],Data[i]["Geometry"],Data[i]["Atom_types"] = xyz_parse(i)

        # Generate adjacency table
        Data[i]["Adj_mat"] = Table_generator(Data[i]["Elements"],Data[i]["Geometry"])

        # If the types weren't present in the supplied xyz file, then automatically parse the atom types based on the adjacency matrix and gens argument. 
        if Data[i]["Atom_types"] == []:
            
            # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
            Hybridizations = Hybridization_finder(Data[i]["Elements"],Data[i]["Adj_mat"])
    
            # Find atom types
            print "Atomtypes could not be parsed from file {}, Determining atom types based on a {}-bond deep search...".format(i,gens)
            Data[i]["Atom_types"] = id_types(Data[i]["Elements"],Data[i]["Adj_mat"],gens,Hybridizations,Data[i]["Geometry"])

        # Check the number of molecules
        mol_in_in = mol_count(Data[i]["Adj_mat"])
        if mol_in_in > 1:
            print "ERROR: {} molecules were discovered in geometry {}. Check the geometry of the input file. Exiting...".format(mol_in_in,i)
            quit()

        # Generate list of bonds angles and dihedrals    
        print "\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Parsing Modes and FF Information for Molecule {}".format(i),"*"*167)
        Data[i]["Bonds"],Data[i]["Bond_types"],Data[i]["Bond_params"],Data[i]["Angles"],Data[i]["Angle_types"],Data[i]["Angle_params"],\
        Data[i]["Dihedrals"],Data[i]["Dihedral_types"],Data[i]["Dihedral_params"],Data[i]["Charges"],Data[i]["Masses"],Data[i]["VDW_params"] =\
            Find_parameters(Data[i]["Adj_mat"],Data[i]["Geometry"],Data[i]["Atom_types"],FF_db=FF_all)

        # Print System characteristics
        print "\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Mode Summary for Molecule {}".format(i),"*"*167)
        print "\nAtom_types ({}):\n".format(len(set(Data[i]["Atom_types"])))
        for j in sorted(set(Data[i]["Atom_types"])):
            print "\t{}".format(j)
        print "\nBond types ({}):\n".format(len(set(Data[i]["Bond_types"])))
        for j in sorted(set(Data[i]["Bond_types"])):
            print "\t{}".format(j)
        print "\nAngle types ({}):\n".format(len(set(Data[i]["Angle_types"])))
        for j in sorted(set(Data[i]["Angle_types"])):
            print "\t{}".format(j)
        print "\nDihedral types ({}):\n".format(len(set(Data[i]["Dihedral_types"])))
        for j in sorted(set(Data[i]["Dihedral_types"])):
            print "\t{}".format(j)

        print "\n{}".format("*"*167)
        print "* {:^163s} *".format("System Characteristics")
        print "*{}*".format("-"*165)
        print "* {:<87s} {:<24s} {:<24s} {:<25s} *".format("Type","Element","Mass","Charge")
        print "{}".format("*"*167)
        for j in range(len(Data[i]["Atom_types"])):
            print " {:<88s} {:<23s} {:< 24.6f} {:< 24.6f}".format(Data[i]["Atom_types"][j],Data[i]["Elements"][j],Data[i]["Masses"][str(Data[i]["Atom_types"][j])],Data[i]["Charges"][j])

    return Data

# A simple wrapper function for generating a fix statement for scaled charges
def scale_charges(charge_scale,Atom_type_dict,Atom_types,Charges):

    fixes = []
    if charge_scale != 1.0:
        fixes = '# Setting up fixes to scale the charges by {}\n'.format(charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in Atom_type_dict.keys() ]):
            fixes += 'group charge_{} type {}\n'.format(i,i)        
        for i in sorted([ Atom_type_dict[j] for j in Atom_type_dict.keys() ]):
            index = Atom_types.index([ k for k in Atom_type_dict.keys() if Atom_type_dict[k] == i ][0])
            fixes += 'variable charge_{} equal {}\n'.format(i,Charges[index]*charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in Atom_type_dict.keys() ]):        
            index = Atom_types.index([ k for k in Atom_type_dict.keys() if Atom_type_dict[k] == i ][0])
            fixes += 'fix {} charge_{} adapt 0 atom charge v_charge_{} reset yes\n'.format(i,i,i)
        fixes += '\n'
    return fixes

# A Simple wrapper function for writing the unscaled pair parameters
def Write_pairs(Filename,Atom_types,VDW_params,eps_scale,sigma_scale):

    # EXACT DUPLICATE OF THE LOOP IN THE WRITE_DATA FUNCTION
    # Create type dictionaries (needed to convert each atom,bond,angle, and dihedral type to consecutive numbers as per LAMMPS convention)
    # Note: LAMMPS orders atomtypes, bonds, angles, dihedrals, etc as integer types. Each of the following dictionaries holds the mapping between
    #       the true type (held in the various types lists) and the lammps type_id, obtained by enumerated iteration over the respective set(types).
    Atom_type_dict = {}
    for count_i,i in enumerate(sorted(set(Atom_types))):
        for j in Atom_types:
            if i == j:
                Atom_type_dict[i]=count_i+1
            if i in Atom_type_dict.keys():
                break

    # Write the settings file
    with open(Filename+'/'+Filename+'.pairs','w') as f:

        # Write non-bonded interactions (the complicated form of this loop is owed
        # to desire to form a nicely sorted list in terms of the lammps atom types
        # Note: Atom_type_dict was initialize according to sorted(set(Atom_types) 
        #       so iterating over this list (twice) orders the pairs, too.
        f.write("     {}\n".format("# Unscaled non-bonded interactions used for fitting vdw parameters"))
        for count_i,i in enumerate(sorted(set(Atom_types))):     
            for count_j,j in enumerate(sorted(set(Atom_types))): 

                # Skip duplicates
                if count_j < count_i:
                    continue

                # Conform to LAMMPS i <= j formatting
                if Atom_type_dict[i] <= Atom_type_dict[j]:
                    f.write("     {:20s} {:<10d} {:<10d} ".format("pair_coeff",Atom_type_dict[i],Atom_type_dict[j]))
                else:
                    f.write("     {:20s} {:<10d} {:<10d} ".format("pair_coeff",Atom_type_dict[j],Atom_type_dict[i]))

                # Determine key (ordered by initialize_VDW function such that i > j)
                if i > j:
                    key = (i,j)
                else:
                    key = (j,i)

                # Write the parameters
                f.write("{:20s} {:< 20.6f} {:< 20.6f}\n".format(VDW_params[key][0],VDW_params[key][1]/eps_scale,VDW_params[key][2]/sigma_scale))

# Description: This functions uses some simple heuristics to determine each atom's
#              hybridization based on the elements types and adjacency matrix
def Hybridization_finder(Elements,Adj_mat):

    Hybridizations = ['X']*len(Elements)
    fail_flag = 0
    
    for i in range(len(Hybridizations)):

        # Find connected atoms and coordination number
        Coord_ind = [ count_j for count_j,j in enumerate(Adj_mat[i]) if j == 1 ]
        Coord_number = len(Coord_ind)

        # Based upon element and coordination number determine hybridization
        if Elements[i] in ['B','Al','C','Si','Ge']:
            if Coord_number == 2:
                Hybridizations[i] = 'sp'
            elif Coord_number == 3:
                Hybridizations[i] = 'sp2'
            elif Coord_number == 4:
                Hybridizations[i] = 'sp3'

        elif Elements[i] in ['N','P','As']:
            if Coord_number == 1:
                Hybridizations[i] = 'sp'
            elif Coord_number == 2:
                Hybridizations[i] = 'sp2'
            elif Coord_number == 3:
                Hybridizations[i] = 'sp3'

        elif Elements[i] in ['O','S','Se']:
            if Coord_number == 1:
                Hybridizations[i] = 'sp2'
            elif Coord_number == 2:
                Hybridizations[i] = 'sp3'

        elif Elements[i] in ['H','Li','Na','K','Rb','Cs','Fr','Be','Mg','Ca','Sr','Ba','Ra']:
            if Coord_number == 1:
                Hybridizations[i] = 's' 

        elif Elements[i] in ['F','Cl','Br','I']:
            if Coord_number == 1:
                Hybridizations[i] = 'sp3'

        else: 
            print "ERROR: hybridization of atom {} could not be determined! Exiting...".format(i)
            fail_flag = 1

    if fail_flag == 1:
        quit()

    return Hybridizations

# Description: This function takes in an element list and adjacency matrix and returns 
#              an index matches list of topological atom types. Here "topological" means
#              the bonding configurations around each atom, going out to "gens" bonds deep.
#              The function uses a canonicalization algorithm to ensure that the same 
#              topological name is returned regardless of the element ordering. 
# 
# Inputs:      Elements:   A list of atomic symbols corresponding to the atom types in the geometry
#                          This list is indexed matched to the adjacency matrix.
#              Adj_mat:    A matrix where each column corresponds to an atom in the geometry, 1's 
#                          correspond to bonds between atoms and 0's to no bonds between atoms.
#
# Returns      Atom_types: A list of topologically named atom types                           
def id_types(Elements,Adj_mat,gens=3,hybridizations=[],geo=[]):

    # Initialize periodic table
    periodic = { "h": 1,  "he": 2,\
                 "li":3,  "be":4,                                                                                                      "b":5,    "c":6,    "n":7,    "o":8,    "f":9,    "ne":10,\
                 "na":11, "mg":12,                                                                                                     "al":13,  "si":14,  "p":15,   "s":16,   "cl":17,  "ar":18,\
                 "k":19,  "ca":20,  "sc":21,  "ti":22,  "v":23,  "cr":24,  "mn":25,  "fe":26,  "co":27,  "ni":28,  "cu":29,  "zn":30,  "ga":31,  "ge":32,  "as":33,  "se":34,  "br":35,  "kr":36,\
                 "rb":37, "sr":38,  "y":39,   "zr":40,  "nb":41, "mo":42,  "tc":43,  "ru":44,  "rh":45,  "pd":46,  "ag":47,  "cd":48,  "in":49,  "sn":50,  "sb":51,  "te":52,  "i":53,   "xe":54,\
                 "cs":55, "ba":56,            "hf":72,  "ta":73, "w":74,   "re":75,  "os":76,  "ir":77,  "pt":78,  "au":79,  "hg":80,  "tl":81,  "pb":82,  "bi":83,  "po":84,  "at":85,  "rn":86}

    # Formatted as a list of lists. Each list corresponds to a generation in the tree
    # and each element in the list is a three element tuple whose first index points to its
    # its parent in the previous generation, whose second index is its index in the
    # Adj_mat/Elements list, and whose third element is its atomic number. The zeroth generation is the atom itself
    trees = [ [ [] ] for i in Elements ]

    for count_i,i in enumerate(Adj_mat):

        # Initialize a list of added atoms so that none get added twice (looping topologies)
        added_atoms = [count_i]

        # Fill in first slot with the parent atom's index
        trees[count_i][0] = [(count_i,count_i,periodic[Elements[count_i].lower()])]

        # Create a list for the first generation connections
        trees[count_i].append([])

        # Iterate over connections and add them to the first generation list
        for count_a,a in enumerate(i):
            
            # If a connection exists add the index to the list
            if a == 1 and count_a not in added_atoms:
                trees[count_i][1] += [(count_i,count_a,periodic[Elements[count_a].lower()])]
                added_atoms += [count_a]

        # Recursively add branch information based on the number of generations requested by user
        for g in range(1,gens):
        
            # Create a list for the second generation connections
            trees[count_i].append([])

            # Iterate over the adj_mat rows for the first generation connections
            for count_j,j in enumerate(trees[count_i][g]):
                for count_a,a in enumerate(Adj_mat[j[1]]):

                    # If a connection exists that is not the connection from the parent generation add the index to the list
                    if a == 1 and count_a != j[0] and count_a not in added_atoms:
                        trees[count_i][g+1] += [(j[1],count_a,periodic[Elements[count_a].lower()])]
                        added_atoms += [count_a]

    # Sort trees: each generation is sorted based on atomic number,
    # where equal, sort is based on the number of attached atoms, 
    # where equal, sort is based on the atomic number of attached atoms
    # where equal, sort continues based on these criteria through all generations
    # until a point of difference arises. If none is found, ordering is left as is.
    for count_i,i in enumerate(trees):
        trees[count_i] = sort_tree(i)

    # Find E/Z stereocenters (if they exist and the necessary 
    # hybridizations and geometry have been supplied) 
    if hybridizations != [] and geo != []:
        EZ_labels = find_EZ(trees,hybridizations,geo)

    # Nest the trees according to topology 
    nested_trees = gen_nested_trees(trees)

    # Created nested atom lists according to topology
    # (same structure as nested tree, except tuples are 
    #  replaced with the atomic number only)
    nested_atoms = list(nested_trees)
    for count_i,i in enumerate(nested_atoms):
        nested_atoms[count_i] = strip_trees(i)

    # Create atom_type labels from the nested atom lists
    atom_types = ['X']*len(nested_atoms)        
    for count_i,i in enumerate(nested_atoms):
        atom_types[count_i] = '['.join(str(i).split(', ['))
        atom_types[count_i] = EZ_labels[count_i] + atom_types[count_i]

    return atom_types

# Sort trees: canonicalization algorithm follows cahn-ingold-prelog rules
# each generation is sorted based on atomic number,
# where equal, each atom is replaced by a list of attached atoms and comparisons are made between each sublist
# where equal, each sublist is replaced by a list of attached atoms in the next generation
# this process continues until a difference is found. 
# Note 1: the common multiple bond rule is omitted as it is unneccesary for canonicalization)
# Note 2: there is subtlety when applying the CIP rules to highly branched structures that is associated 
#         with the order in which subshells are compared. Neither Mcarrey, Wade, or Sorrell's books
#         are explicit on the matter. Here I always compare lists based upon highest MW preceeding
#         atom.
def sort_tree(tree):

    # Sort all generations by atomic weight
    for count_i,i in enumerate(tree):
        tree[count_i].sort(key=lambda x:x[2],reverse=True)

    # Iterate over each generation ranking each node's priority
    for count_i,i in enumerate(tree):        

        # If there are empty generations, avoid sorting
        if len(i) == 0: continue

        # Find nodes with identical atomic weights and the number of connections of each node in the next generation
        # Initialize rank list
        rank = [-1]*len(i)
        rank_counter = 0
        
        for j in range(len(i)):

            # If unassigned, assign a rank
            if rank[j] == -1:
                rank[j] = rank_counter
                
            # Check for matches based on atomic number and assign equal ranks to those
            for k in range(len(i)):
                if i[j][2] == i[k][2] and j != k and rank[k] == -1:
                    rank[k]=rank_counter

            # Increment rank_counter based on the number of unresolved pairs that had the same priority
            rank_counter += len([ k for k in rank if k == rank_counter ])

        # Recursively search through the depth of the tree to resolve priority differences
        depth = count_i
        for j in set(rank):

            # Find connections that need a deeper comparison to resolve priority
            idx = [ count_k for count_k,k in enumerate(rank) if k == j ]

            # Skip already resolved ranks
            if len(idx) == 1:
                continue

            # Initialize tmp_ranks (holds the sublist of ranks for the atoms being resolved)
            tmp_ranks = [0]*len(idx)
            
            # Initialize the sorting list (a list of lists holding sublists of connections for each connection)
            # i.e., the elements in the first list correspond to lists holding the attached atoms
            #       the elements in each sublist hold sublists to compare branches between generations. 
            sorting_list = [ [find_connections(tree,[depth,k])] for k in idx ]

            # Recursively search the tree for differences in priority. 
            # The loop is exited once all of the unresolved branches have been assigned a rank
            while len(idx) > 0:
                                
                # rank the sorting_list (list of lists of lists; each element in the first list
                disp = 0
                for k in sorted(set(tmp_ranks)):

                    # This set of loops guarantees that comparisons are only performed between equally ranked subsets of 
                    # the unresolved set (e.g,, there might be 4 branches that are still unresolved, but 2 are ranked 0,0
                    # and 2 are ranked 2,2. It is important in this case that 0,0 and 2,2 get compared only amongst themselves.
                    loc_list = []
                    sub_sorting_list = []
                    for count_m,m in enumerate(tmp_ranks):
                        if m == k:
                            sub_sorting_list += [sorting_list[count_m]]
                            loc_list += [count_m]

                    # rank_nested_lists: returns the priorities/ranks of the sub_sorting_list
                    sub_ranks = rank_nested_lists(sub_sorting_list)
                    for count_m,m in enumerate(loc_list):
                        tmp_ranks[m] = sub_ranks[count_m]+disp 
                        
                    # the minimum ranking of the next atoms being resolved starts at disp + #of_just_resolved
                    disp += len(sub_ranks)
                
                # Update rank for resolved atoms
                for count_k,k in enumerate(tmp_ranks):
                    if tmp_ranks.count(k) == 1:
                        rank[idx[count_k]] += tmp_ranks[count_k]

                # Update the rank of unresolved atoms:
                # the overall rank of each unresolved atom is incremented for each just
                # resolved atom that is ranked in front of it.
                for count_k,k in enumerate(tmp_ranks):
                    if tmp_ranks.count(k) == 1:
                        for count_m,m in enumerate(tmp_ranks):
                            if tmp_ranks.count(m) != 1 and m > k:
                                rank[idx[count_m]] += 1

                # ID the remaining unresolved branches and remove the
                # rest from sorting_list and idx
                tmp_sorting_list = []
                tmp_idx = []                
                tmp_tmp_ranks = []
                for count_k,k in enumerate(tmp_ranks):
                    if tmp_ranks.count(k) != 1:
                        tmp_sorting_list += [sorting_list[count_k]]
                        tmp_idx += [idx[count_k]]
                        tmp_tmp_ranks += [k]
                sorting_list = tmp_sorting_list
                idx = tmp_idx
                tmp_ranks = tmp_tmp_ranks

                # Repopulate sorting_list from the next generation
                tmp_sorting_list = [ [] for k in range(len(idx)) ]
                for count_k,k in enumerate(sorting_list):
                    for count_m,m in enumerate(k):
                        for count_n,n in enumerate(m):

                            # if n is an empty node then each subranch mutually deadended here and it can be skipped
                            if n == []:
                                continue
                            tmp_sorting_list[count_k] += [ find_connections(tree,n[1]) ]
                sorting_list = tmp_sorting_list

        # NEEDS TO BE MORE THOROUGHLY ADDRESSED (ENCOUNTERED ERROR IN FITTING CARBONATES)
        if len(set(rank)) != len(rank):
            missing = [ j for j in range(len(rank)) if j not in rank ]
            for j in missing:
                rank[rank.index(j+1)] = j
        
        # Sort the tree based on the final rankings
        tmp = [[]]*len(tree[count_i])
        for count_j,j in enumerate(rank):
            tmp[j] = tree[count_i][count_j]
        tree[count_i] = tmp

    return tree

def rank_nested_lists(lol):    

    # matrix based sort algorithm. This is knowingly suboptimal, but it 
    # gives much more flexibility during the sort than standard algorithms. 
    order=zeros([len(lol),len(lol)])
    rank = [0]*len(lol)
    terminal_list = []
    for count_i,i in enumerate(lol):
        
        # If there were no connections in the current generation then comparisons have 
        # been exhausted and ranking should be set in this round.
        if len(flatten(i)) == 0:
            terminal_list += [count_i]
            order[count_i,:] = 0
            continue

        for count_j,j in enumerate(lol):

            if count_j > count_i:

                flag=0

                # iterate over sublists in i while performing
                # comparisons with corresponding sublists in j
                for count_k,k in enumerate(i):

                    if flag == 1:
                        break

                    # If there are no sublists left in j and difference has been found, i has priority
                    elif count_k > len(j)-1:
                        order[count_i,count_j]=1
                        order[count_j,count_i]=0
                        break

                    for count_m,m in enumerate(k):

                        # If the sublist in j has run out of values, i has priority
                        if count_m > len(j[count_k])-1:
                            order[count_i,count_j]=1
                            order[count_j,count_i]=0
                            flag = 1
                            break                       

                        # If the sublist element is an empty node (with branching empty nodes can be mixed with full nodes)
                        if i[count_k][count_m] == [] and j[count_k][count_m] != []:
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break

                        # If the sublist element is an empty node (with branching empty nodes can be mixed with full nodes) (inverse)
                        if i[count_k][count_m] != [] and j[count_k][count_m] == []:
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break

                        # If both sublist elements are empty nodes (with branching empty nodes can be mixed with full nodes)
                        if i[count_k][count_m] == [] and j[count_k][count_m] == []:
                            continue

                        # If there is a difference in the sublist elements assign priority
                        elif i[count_k][count_m][0] > j[count_k][count_m][0]:
                            order[count_i,count_j]=1
                            order[count_j,count_i]=0
                            flag = 1
                            break

                        # If there is a difference in the sublist elements assign priority (inverse)
                        elif i[count_k][count_m][0] < j[count_k][count_m][0]:
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break

                        # If the sublist in i has run out of values and j still has attachments, j has priority
                        elif count_m == len(i[count_k])-1 and len(j[count_k]) > len(i[count_k]):
                            order[count_i,count_j]=0
                            order[count_j,count_i]=1
                            flag = 1
                            break
                            
                    # i has run out of sublists without finding a difference and j still has lists, j has priority
                    if count_k == len(i)-1 and len(j) > len(i) and flag == 0:
                        order[count_i,count_j]=0
                        order[count_j,count_i]=1                        
                        break

                    # if both i and j have run out of sublists without finding a difference then they are identical up to this point
                    elif count_k == len(i)-1 and len(i) == len(j) and flag == 0:
                        order[count_i,count_j]=1
                        order[count_j,count_i]=1
                        break
                    
    # Assign ranks based on entries in the order matrix
    for count_i,i in enumerate(order):
        rank[count_i]=int(len(order)-sum(i)-1)

    # Assign unique ranks for branches that are at their terminus
    # Reverse order shouldn't matter, but it guarrantees that the
    # first discovered gets higher priority
    for count_i,i in enumerate(terminal_list[::-1]):
        rank[i]=rank[i]-count_i

    return rank

def find_connections(tree,loc):
    if len(tree) > loc[0]+1:
        
        # Find connections
        connections = [ [i[2],[loc[0]+1,count_i]] for count_i,i in enumerate(tree[loc[0]+1]) if i[0] == tree[loc[0]][loc[1]][1] ]

        # If no connections are found then a 0-entry is returned
        if connections == []:
            connections = [[]]
        # If connections are found then proceed with sort
        else:
            connections.sort(key=lambda x:x[0],reverse=True)

    # If a next generation doesn't exist then a 0-entry is returned
    else:
        connections = [[]]

    return connections

def flatten(lol):
    flat = []
    for i in lol:
        if type(i) == type([]):
            flat += flatten(i)
        else:
            flat += [i]
    return flat
    
# lol: nested list of lists, where each element is a 3-element tuple
# This function returns a list of lists with the same structure as lol, except 
# that each tuple is replaced by its third element (atomic number in this case)
def strip_trees(lol):
    tmp = []
    for count_i,i in enumerate(lol):
        if type(i) == type([]):
            tmp += [strip_trees(i)]
        else:
            tmp += [i[2]]
    return tmp

# Description: iterate over the nodes in each generation and search the previous generation for matches (using find_loc function)
# insert all matches into the label as a list in the index following the previous generation node.
# Operates on the total tree structure (so, trees is a list of trees)
def gen_nested_trees(trees):
    nested_trees = [ [] for i in range(len(trees)) ]
    for count_i,i in enumerate(trees):

        # Count the number of occupied generations
        num_gens = len([ count_j for count_j,j in enumerate(i) if len(j)>0])

        # Intiialize the nested list with the first generation atom
        nest = [i[0][0]]

        # Loop over generations and unwind the topology
        for j in range(1,num_gens):

            # Reversed list keeps the ordering so that highest priority ends up at the front
            for k in i[j][::-1]:         
                loc = find_loc(k[0],nest)
                if len(loc) == 1:
                    nest.insert(loc[0]+1,[k])
                elif len(loc) == 2:
                    nest[loc[0]].insert(loc[1]+1,[k])
                elif len(loc) == 3:
                    nest[loc[0]][loc[1]].insert(loc[2]+1,[k])
                elif len(loc) == 4:
                    nest[loc[0]][loc[1]][loc[2]].insert(loc[3]+1,[k])
                elif len(loc) == 5:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]].insert(loc[4]+1,[k])
                elif len(loc) == 6:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]].insert(loc[5]+1,[k])
                elif len(loc) == 7:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]].insert(loc[6]+1,[k])
                elif len(loc) == 8:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]][loc[6]].insert(loc[7]+1,[k])
                elif len(loc) == 9:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]][loc[6]][loc[7]].insert(loc[8]+1,[k])
                elif len(loc) == 10:
                    nest[loc[0]][loc[1]][loc[2]][loc[3]][loc[4]][loc[5]][loc[6]][loc[7]][loc[8]].insert(loc[9]+1,[k])
                else:
                    print "Nesting is too deep, ensure that the requested generations do not exceed 10"
                    quit()
        nested_trees[count_i] = nest

    return nested_trees

# Description: Recursive search to id location of the matching node in the list of lists.
#              Returns a list of indices corresponding to the position within the list of
#              lists of the matching element. (e.g., loc=[2,4,3] means lol[2][4][3] is the match)
#
# Inputs:      idx: scalar, sought match for the second element in the tuples held within the list of lists.
#              lol: list of lists (arbitrary nesting) whose lowest level elements are 3-tuples. 
#
# Returns:     loc: a list of indices specifying the location in lol of the matching tuple
def find_loc(idx,lol):
    loc = []                                   # Initialize list to hold the location indices
    for count_i,i in enumerate(lol):           # Iterate over the elements in the list
        if type(i) == type([]):                # If a list is encountered a recursion is initiated
            tmp = find_loc(idx,i)              # recursive search for a match (yields [] is no match is found)
            if tmp != []:                      # If the recursion yields a match then we have a winner
                loc = loc + [count_i] + tmp    # append the list index (count_i) and the loc list returned by the recursion
        elif idx == i[1]:                      # If a match is found in the element-wise search append its location
            loc += [count_i]      
    return loc

def find_EZ(trees,hybridizations,geo):
    
    EZ_labels = ['']*len(trees)
    for count_i,i in enumerate(trees):

        # Check if the anchor atom is a potential stereocenter
        # NOTE: in the current implementation only C,N,Si,and P are considered potential EZ centers)
        if hybridizations[count_i] != 'sp2' or i[0][0][2] not in [6,7,14,15]:
            continue

        # Find indices of potential EZ bonds
        idx = [ j[1] for j in i[1] if hybridizations[j[1]] == 'sp2' and j[2] in [6,7,14,15] ]

        # Iterate over all of the EZ bonds and determine stereochemistry
        for j in idx:

            # Find anchor atom's highest priority connection (excluding the sp2 bonded group)      
            ind_1 = [ k[1] for k in i[1] if k[1] != j ][0]
            
            # Find bonded atom's highest priority connection (excluding the sp2 bonded group)
            ind_4 = [ k[1] for k in trees[j][1] if k[1] != count_i ][0] 

            # Intialize the coordinate array of the relevant dihedral
            xyz = zeros([4,3])
            xyz[0] = geo[ind_1,:]
            xyz[1] = geo[count_i,:]
            xyz[2] = geo[j,:]
            xyz[3] = geo[ind_4,:]

            # Calculate dihedral using call to dihedral_calc
            angle = 180.0*dihedral_calc(xyz)/pi         # Calculate dihedral
            angle = int(180*round(abs(angle)/180.0))    # Round to 0 or 180

            # ID E/Z based on angle
            if angle == 0:
                EZ_labels[count_i] += "Z"
            elif angle == 180:
                EZ_labels[count_i] += "E"

    return EZ_labels

# Description: Calculates the dihedral angle (in radians) for
#              a quadruplet of atoms
#
# Inputs:      xyz: a 4x3 numpy array, where each row holds the 
#                   cartesian position of each atom and row indices
#                   determine placement in the dihedral ( i.e, which
#                   atoms correspond to 1-2-3-4 in the dihedral )
# 
# Returns:    angle: dihedral angle in radians
#
def dihedral_calc(xyz):
    
    # Calculate the 2-1 vector           
    v1 = (xyz[1]-xyz[0]) 
                                                             
    # Calculate the 3-2 vector           
    v2 = (xyz[2]-xyz[1]) 
                                                             
    # Calculate the 4-3 bond vector      
    v3 = (xyz[3]-xyz[2]) 

    # Calculate dihedral (funny arctan2 formula avoids the use of arccos and retains the sign of the angle)
    angle = arctan2( dot(v1,cross(v2,v3))*(dot(v2,v2))**(0.5) , dot(cross(v1,v2),cross(v2,v3)) )
    
    return angle

def prepare_data(old_folder,new_folder):

    # Check for the presence of a cycle trajectory
    cycle_dir = natural_sort([ i for i in next(os.walk('.'))[1] if'cycle-' in i ])
    if len(cycle_dir) > 0:
        cycle_dir = cycle_dir[-1]
        new_params = [ "{}/{}".format(cycle_dir,i) for i in os.listdir(cycle_dir) if os.path.isfile("{}/{}".format(cycle_dir,i)) == True and "AA.db" in i ]
        if len(new_params) > 0:
            new_params = new_params[-1]

    # check several existence conditions before starting the copying
    if os.path.isfile("{}/{}.end.data".format(old_folder,old_folder)) == False:
        print "ERROR: the file {}.end.data is missing from the previous run in folder {}. Cannot restart the run. Exiting...".format(old_folder,old_folder)
        quit()
    elif os.path.isfile("{}/{}.map".format(old_folder,old_folder)) == False:
        print "ERROR: the file {}.map is missing from the previous run in folder {}. Cannot restart the run. Exiting...".format(old_folder,old_folder)
        quit()
    elif os.path.isfile("{}/{}.in.settings".format(old_folder,old_folder)) == False:
        print "ERROR: the file {}.in.settings is missing from the previous run in folder {}. Cannot restart the run. Exiting...".format(old_folder,old_folder)
        quit()
    elif len(cycle_dir) == 0:
        print "ERROR: No 'cycle-' diretories were found in the folder {}. Cannot start the run because there are no VDW parameters. Exiting...".format(old_folder)
    elif len(new_params) == 0:
        print "ERROR: No file ending in *AA.db was found in the folder {}. Cannot start the run because there are no VDW parameters. Exiting...".format(cycle_dir)

    # Copy files and generate taffy2lammps mapping
    copyfile("{}/{}.end.data".format(old_folder,old_folder),"{}/{}.data".format(new_folder,new_folder))
    copyfile("{}/{}.map".format(old_folder,old_folder),"{}/{}.map".format(new_folder,new_folder))
#    copyfile("{}/{}.pairs".format(old_folder,old_folder),"{}/{}.pairs".format(new_folder,new_folder))
    taffy2lammps_dict, lammps2taffy_dict = taffy2lammps("{}/{}.map".format(new_folder,new_folder),"{}/{}.data".format(new_folder,new_folder))

    # Parse the vdw parameters that need to be updated in the *in.settings file
    new_vdw = {}
    with open(new_params,'r') as f:
        for lines in f:
            fields = lines.split()
            if len(fields) > 0 and fields[0] == "vdw" and fields[1] in taffy2lammps_dict.keys() and fields[2] in taffy2lammps_dict.keys():
                new_vdw[(taffy2lammps_dict[fields[1]],taffy2lammps_dict[fields[2]])] = [fields[3]] + [ float(i) for i in fields[4:] ]
                new_vdw[(taffy2lammps_dict[fields[2]],taffy2lammps_dict[fields[1]])] = [fields[3]] + [ float(i) for i in fields[4:] ]

    # Write new in.settings file with the updated parameters
    dihedral_styles = []
    with open("{}/{}.in.settings".format(old_folder,old_folder),'r') as o:
        with open("{}/{}.in.settings".format(new_folder,new_folder),'w') as f:
            for lines in o:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "pair_coeff" and (fields[1],fields[2]) in new_vdw:
                    if new_vdw[(fields[1],fields[2])][0] == 'lj':                        
                        f.write("     {:20s} {:<10s} {:<10s} {:20s} {: <20.6f} {: <20.6f}\n".format(fields[0],fields[1],fields[2],fields[3],new_vdw[(fields[1],fields[2])][1],new_vdw[(fields[1],fields[2])][2]))
                else:
                    f.write(lines)
                    if len(fields) > 0 and fields[0] == "dihedral_coeff": dihedral_styles += [fields[2]]

    # Write new pairs file with the updated parameters
    with open("{}/{}.pairs".format(old_folder,old_folder),'r') as o:
        with open("{}/{}.pairs".format(new_folder,new_folder),'w') as f:
            for lines in o:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "pair_coeff" and (fields[1],fields[2]) in new_vdw:
                    if new_vdw[(fields[1],fields[2])][0] == 'lj':                        
                        f.write("     {:20s} {:<10s} {:<10s} {:20s} {: <20.6f} {: <20.6f}\n".format(fields[0],fields[1],fields[2],fields[3],new_vdw[(fields[1],fields[2])][1],new_vdw[(fields[1],fields[2])][2]))
                else:
                    f.write(lines)

    return cycle_dir.split('-')[1],list(set(dihedral_styles))
                    
# Generate dictionaries for mapping between taffy atomtypes and lammps atomtypes
def taffy2lammps(mapfile,datafile):

    # Check existence conditions for the input files
    if os.path.isfile(mapfile) == False:
        print "ERROR in taffy2lammps: there is no mapfile named {}. Exiting...".format(mapfile)
        quit()
    elif os.path.isfile(datafile) == False:
        print "ERROR in taffy2lammps: there is no datafile named {}. Exiting...".format(datafile)
        quit()

    # Get taffy types
    with open(mapfile,'r') as f:
        for lc,lines in enumerate(f):
            if lc == 0:
                N_atoms = int(lines.split()[0])
                taffytypes = ["X"]*N_atoms
                lammpstypes = ["X"]*N_atoms
            if lc > 1:
                taffytypes[lc-2] = lines.split()[0]
            if lc == N_atoms+1:
                break
    
    # Get LAMMPS types
    flag = 0
    count = 0
    with open(datafile,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            if flag == 0 and len(fields) > 0 and fields[0] == "Atoms":
                flag = 1
                continue
            if flag > 0: 
                if len(fields) == 0:
                    flag += 1            
                    if flag == 3:
                        break
                else:
                    lammpstypes[int(fields[0])-1] = fields[2]
                    count += 1

    # Return dictionaries for mapping taffy onto lammps types and the reverse
    return { i:lammpstypes[taffytypes.index(i)] for i in set(taffytypes) }, { j:taffytypes[lammpstypes.index(j)] for j in set(lammpstypes) } 

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/gen_md_for_vdw.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

if __name__ == "__main__":
   main(sys.argv[1:])
