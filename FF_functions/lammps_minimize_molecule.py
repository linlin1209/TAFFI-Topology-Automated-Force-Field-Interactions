#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime
from numpy import *
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from extract_charges import *
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *
import random
from gen_md_for_sampling import get_data, Find_parameters, Write_data

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file and writes inputs for a LAMMPS job for a single molecule shrink-wrapped simulation.'
                                                 'These simulations are meant to be used in conjunction with free-energy trajectories for subtracting out the intramolecular '+\
                                                 'electrostatics from the final results.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_file', help = 'A string with the name of the  *.xyz file to be geometry optimized.')

    parser.add_argument('FF_files', help = 'A quoted list of the force-field files. '+\
                                           'Formatting of the force-field files is as produced by the FF_gen.py/FF_extract.py programs)')

    #optional arguments    
    parser.add_argument('-f', dest='frequency', default=1000,
                        help = 'Controls the sampling frequency during the MD simulation. (default: 1000)')

    parser.add_argument('-o', dest='outputname', default='minimize',
                        help = 'Sets the output filename prefix (default: minimize)')

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

    parser.add_argument('-pair_styles', dest='pair_styles', default='lj/cut/coul/cut 14.0 14.0',
                        help = 'Supplies the information for the pair style setting(s). Supplied as a string and should be formatted for comply '+\
                               'with the LAMMPS "pair_style hybrid" command and should match the cutoffs for the soft potential used in the free energy trajectory '+\
                               '(default: "lj/cut/coul/cut 14.0 14.0")')

    parser.add_argument('--keep_intermediates', dest='keep_intermediates', default=False, const=True, action='store_const',
                        help = 'When this flag is enabled the intermediate files are kept. (Default: off)')

    parser.add_argument('-lammps_exe', dest='lammps_exe', default='/home/bsavoie/LAMMPS/exe/lmp_fep_tachus_openmpi_sigfig_rigid',
                        help = 'This holds the location of the lammps executable on your system. LAMMPS is called to perform the MD-constrained minimizations for dihedral fits. (default: sys_dep)')

    parser.add_argument('-q', dest='q', default='none',
                        help = 'Controls the total charge on the molecule. By default, the rounded integer charge is used for each molecule (round). The program expects a list of integers (or none, or round for the default). '+\
                               'If less integers are supplied than the number of coord files, then the list is automatically expanded using the last supplied (or default) element. (default: None)') 

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'When atomtypes are being automatically determined (i.e. when the supplied *.xyz files do not have atomtype information) this variable controls the bond depth used to define unique atomtypes. (default: 2)')

    parser.add_argument('--impropers', dest='improper_flag', default=False, action='store_const', const=True,
                        help = 'When present, the simulation will use improper dihedral terms. (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)

    # Converting each list input argument. 
    args.coord_file = args.coord_file.split()
    args.FF_files = args.FF_files.split()
    args.frequency = int(float(args.frequency))
    args.eps_scale = float(args.eps_scale)
    args.sigma_scale = float(args.sigma_scale)
    args.charge_scale = float(args.charge_scale)
    args.gens = int(args.gens)
    args.q = [str(args.q)]

    # Check that the input is an .xyz file. 
    if len(args.coord_file) > 1:
        print("ERROR in lammps_minimize_molecule: Only one geometry can be supplied. Exiting...")
        quit()
    elif args.coord_file[0].split('.')[-1] != 'xyz':
        print("ERROR in lammps_minimize_molecule: Check to ensure that the input coordinate file is in .xyz format. Exiting...")
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
        print("PROGRAM CALL: python gen_md_for_vdw.py {}\n".format(' '.join([ i for i in argv])))
        
    # Check that the supplied FF files exist
    for i in args.FF_files:
        if os.path.isfile(i) != True:
            print("ERROR: Specified FF_file ({}) doesn't exist. Please check the path. Exiting...".format(i))
            quit()
            
    # Catenate the FF files and copy them to the output folder
    FF_all = Filename+"/"+Filename+'.db'
    with open(Filename+"/"+Filename+'.db','w') as f:
        for i in args.FF_files:
            with open(i,'r') as ff_file:
                for j in ff_file:
                    f.write(j)
            f.write("\n")

    # Grab data for the molecule being added to the MD simulation.
    Data = get_data(FF_all,args.coord_file,[1],args.q,args.gens,Improper_flag=args.improper_flag)
    Data = Data[args.coord_file[0]]

    # Generate VDW parameters
    VDW_params = initialize_VDW(sorted(set(Data["Atom_types"])),sigma_scale=args.sigma_scale,eps_scale=args.eps_scale,VDW_type=args.pair_styles.split()[0],VDW_FF=Data["VDW_params"])

    # Generate Simulation Dictionaries
    # The bond, angle, and diehdral parameters for each molecule are combined into one dictionary
    Bond_params = {}; Angle_params = {}; Dihedral_params = {}; Masses = {}
    for j in list(Data["Bond_params"].keys()): Bond_params[j] = Data["Bond_params"][j]
    for j in list(Data["Angle_params"].keys()): Angle_params[j] = Data["Angle_params"][j]
    for j in list(Data["Dihedral_params"].keys()): Dihedral_params[j] = Data["Dihedral_params"][j]
    for j in list(Data["Masses"].keys()): Masses[j] = Data["Masses"][j]

    # # Write the lammps datafile (the function returns a dictionary, Atom_type_dict, holding the mapping between atom_types and the lammps id numbers; this mapping is needed for setting fixes)
    # Atom_type_dict = Write_data(Filename,Atom_types,Sim_Box,Elements,Geometry,Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,\
    #                             Dihedrals,Dihedral_types,Dihedral_params,Charges,VDW_params,Masses,Molecule)

    # # Write the lammps input files
    # print "Writing LAMMPS input file ({})...".format(Filename+'.in.init')
    # Write_input(Filename,args.T_equil,args.t_equil,args.T_anneal,args.t_anneal,args.frequency,args.onefourscale,args.pair_styles,Dihedral_styles,args.density,args.N_density,fixes)

    #### NEW COMMMANDS ####
    print("Writing LAMMPS datafile ({})...".format(Filename+'.data'))
    print("Writing LAMMPS settings file ({})...".format(Filename+'.in.settings'))

    # Write the lammps datafile (the function returns a dictionary, Atom_type_dict, holding the mapping between
    # atom_types and the lammps id numbers; this mapping is needed for setting fixes)                    

    Sim_Box=array([min(Data['Geometry'][:,0])-0.5,max(Data['Geometry'][:,0])+0.5,\
                   min(Data['Geometry'][:,1])-0.5,max(Data['Geometry'][:,1])+0.5,\
                   min(Data['Geometry'][:,2])-0.5,max(Data['Geometry'][:,2])+0.5])

    Atom_type_dict,fixed_modes = Write_data(Filename,Data["Atom_types"],Sim_Box,Data["Elements"],Data["Geometry"],Data["Bonds"],Data["Bond_types"],Data["Bond_params"],Data["Angles"],Data["Angle_types"],Data["Angle_params"],\
                                            Data["Dihedrals"],Data["Dihedral_types"],Data["Dihedral_params"],Data["Impropers"],Data["Improper_types"],Data["Improper_params"],Data["Charges"],VDW_params,Data["Masses"],\
                                            [0]*len(Data["Elements"]),Improper_flag = args.improper_flag)
    os.chdir(Filename)
#    Atom_type_dict = write_lammps_data('minimize',Data["Atom_types"],Data["Elements"],Data["Geometry"],Data["Bonds"],Data["Bond_types"],Data["Bond_params"],Data["Angles"],Data["Angle_types"],Data["Angle_params"],\
#                                       Data["Dihedrals"],Data["Dihedral_types"],Data["Dihedral_params"],{},[],Data["Charges"],VDW_params,Data["Masses"])

    # Generate fix statement for scaled charges
    fixes = scale_charges(args.charge_scale,Atom_type_dict,Data["Atom_types"],Data["Charges"])

    # Gather the different dihedral_styles
    Dihedral_styles = set([ Dihedral_params[i][0] for i in Data["Dihedral_types"] ])

    # Write the input file for the forward constrained scan starting from configuration z
    write_minimize_input(Filename,100,0.0,"lj/cut/coul/cut 100.0 100.0",Dihedral_styles,fixes)    

    # Run the LAMMPS minimization
    lammps_call = "{} -in {}.in.init -screen {}.screen -log {}.log".format(args.lammps_exe,Filename,"minimize","minimize")
    subprocess.call(lammps_call.split())

    # Parse geometry
    final_geo = parse_lammpstrj('minimize.lammpstrj')

    # Write the final geometry
    xyz_write("minimized.xyz",Data["Elements"],final_geo)

    # Clean up lammps files
    if args.keep_intermediates is False:
        os.remove("{}.xyz".format(Filename))
        os.remove("{}.data".format(Filename))
        os.remove("{}.in.init".format(Filename))
        os.remove("{}.in.settings".format(Filename))
        os.remove("minimize.lammpstrj")
        os.remove("minimize.log")
        os.remove("minimize.screen")
        os.remove(Filename+'.db')

    # return to the original directory
    os.chdir('..')
    
    # Print banner
    print("\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Success! Have a Nice Day!","*"*167))
    return

# Wrapper function for parsing the lammpstrj file and returning the last geometry
def parse_lammpstrj(filename):

    # Parse geometry
    N_atoms   = None
    atom_flag = 0
    counter   = -1
    with open('minimize.lammpstrj','r') as f:
        for lines in f:
            fields = lines.split()

            # skip empty
            if len(fields) == 0: continue

            # parse the number of atoms from the first frame
            if N_atoms is None and len(fields) == 4 and fields[0] == "ITEM:" and fields[1] == "NUMBER" and fields[2] == "OF" and fields[3] == "ATOMS":
                atom_flag = 1
                continue

            # parse the number of atoms
            if atom_flag == 1:
                N_atoms = int(fields[0])
                geo = zeros([N_atoms,3])
                atom_flag = 0
                continue

            # reset counter at the start of the frame
            if fields[0] == "ITEM:" and fields[1] == "TIMESTEP":
                counter = -1

            # trigger the geometry parse
            if fields[0] == "ITEM:" and fields[1] == "ATOMS":
                counter = 0
                continue

            # parse the geometry
            if counter > -1:
                geo[counter] = array([float(fields[2]),float(fields[3]),float(fields[4])])
                counter += 1

    # return geometry
    return geo

# Wrapper function for writing the lammps data file which is read by the .in.init file during run initialization. This function is 
# nearly identical to the function of the same name within the gen_md_for_vdw.py script, except the qc option has been added to be compatible
# with the FF_dict structure and the Molecule variable has been removed because only one molecule is present in each simulation.
def write_lammps_data(Filename,Atom_types,Elements,Geometry,Bonds,Bond_types,Bond_params={},Angles=[],Angle_types=[],Angle_params={},Dihedrals=[],Dihedral_types=[],Dihedral_params={},Dihedral_harmonic_params={},One_fives={},\
                          Charges={},VDW_params={},Masses={},qc="dft",one_five_opt=0,tables=[],Sim_Box=None):
    
    # Write an xyz for easy viewing
    with open(Filename+'.xyz','w') as f:
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
            if i in list(Atom_type_dict.keys()):
                break
    Bond_type_dict = {}
    for count_i,i in enumerate(sorted(set(Bond_types))):
        for j in Bond_types:
            if i == j:
                Bond_type_dict[i]=count_i+1
            if i in list(Bond_type_dict.keys()):
                break
    
    # Add 1-5 bonds as the last type (used for turning off 1-5 non-bonded interactions
    if one_five_opt == 1:
        Bond_type_dict["one_fives"] = count_i+2   

    Angle_type_dict = {}
    for count_i,i in enumerate(sorted(set(Angle_types))):
        for j in Angle_types:
            if i == j:
                Angle_type_dict[i]=count_i+1
            if i in list(Angle_type_dict.keys()):
                break
    Dihedral_type_dict = {}
    for count_i,i in enumerate(sorted(set(Dihedral_types))):
        for j in Dihedral_types:
            if i == j:
                Dihedral_type_dict[i]=count_i+1
            if i in list(Dihedral_type_dict.keys()):
                break

    # Write the data file
    with open(Filename+'.data','w') as f:
        
        # Write system properties
        f.write("LAMMPS data file via vdw_self_gen.py, on {}\n\n".format(datetime.datetime.now()))

        f.write("{} atoms\n".format(len(Elements)))
        f.write("{} atom types\n".format(len(set(Atom_types))))
        if len(Bonds) > 0:
            if one_five_opt == 1:
                f.write("{} bonds\n".format(len(Bonds)+len(One_fives)))
                f.write("{} bond types\n".format(len(set(Bond_types))+1))
            else:
                f.write("{} bonds\n".format(len(Bonds)))
                f.write("{} bond types\n".format(len(set(Bond_types))))
        if len(Angles) > 0:
            f.write("{} angles\n".format(len(Angles)))
            f.write("{} angle types\n".format(len(set(Angle_types))))
        if len(Dihedrals) > 0:
            f.write("{} dihedrals\n".format(len(Dihedrals)))
            f.write("{} dihedral types\n\n".format(len(list(Dihedral_type_dict.keys()))))

        # Write box dimensions
        if Sim_Box is None:
            f.write("{:< 20.16f} {:< 20.16f} xlo xhi\n".format(min(Geometry[:,0])-0.5,max(Geometry[:,0])+0.5))
            f.write("{:< 20.16f} {:< 20.16f} ylo yhi\n".format(min(Geometry[:,1])-0.5,max(Geometry[:,1])+0.5))
            f.write("{:< 20.16f} {:< 20.16f} zlo zhi\n\n".format(min(Geometry[:,2])-0.5,max(Geometry[:,2])+0.5))
        else:
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
                    f.write("{} {:< 12.6f} {:< 12.6f}\n".format(count_i+1,Bond_params[j][1],Bond_params[j][2])) # count_i+1 bc of LAMMPS 1-indexing

        # Add one_five ghost bonds
        if one_five_opt == 1:
            f.write("{} {:< 12.6f} {:< 12.6f}\n".format(count_i+2,0.0,2.0)) # dummy bond for 1-5 interactions
        f.write("\n")

        # Write Angle Coeffs
        f.write("Angle Coeffs\n\n")
        for count_i,i in enumerate(set(Angle_types)):
            for j in set(Angle_types):
                if Angle_type_dict[j] == count_i+1:
                    f.write("{} {:< 12.6f} {:< 12.6f}\n".format(count_i+1,Angle_params[j][1],Angle_params[j][2])) # count_i+1 bc of LAMMPS 1-indexing
        f.write("\n")

        # Write Atoms
        f.write("Atoms\n\n")
        for count_i,i in enumerate(Atom_types):
            f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
            .format(count_i+1,0,Atom_type_dict[i],Charges[count_i],Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],0,0,0))

        # Write Bonds
        bond_count = 1
        if len(Bonds) > 0:
            f.write("\nBonds\n\n")
            for count_i,i in enumerate(Bonds):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(bond_count,Bond_type_dict[Bond_types[count_i]],i[0]+1,i[1]+1))
                bond_count += 1

        # Write 1-5 ghost bonds
        if one_five_opt == 1: 
            for i in One_fives:
                f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(bond_count,Bond_type_dict["one_fives"],i[0]+1,i[-1]+1))
                bond_count += 1

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
    with open(Filename+'.in.settings','w') as f:

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

                # Write the parameters (note: the internal type (e.g., lj) is converted into the corresponding lammps type (e.g., lj/cut/coul/long))
                for k in VDW_params[key]:
                    if k == "lj": k = "lj/cut/coul/cut"
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
            for j in Bond_params[i][1:]:
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
                if type(j) is double:
                    f.write("{:< 20.6f} ".format(j))
                if type(j) is int:
                    f.write("{:< 20d} ".format(j))
            f.write("\n")

        # Write 1-5 ghost bond parameter
        if one_five_opt == 1:
            f.write("     {:20s} {:<10d} {:< 20.6f} {:< 20.6f}\n".format("bond_coeff",Bond_type_dict["one_fives"],0.0,2.0))

        # Write bending interactions
        # Note: Angle_type_dict was initialized by looping over sorted(set(Angle_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Bending interactions"))
        for i in sorted(set(Angle_types)):
            f.write("     {:20s} {:<10d} ".format("angle_coeff",Angle_type_dict[i]))
            for j in Angle_params[i][1:]:
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
                if type(j) is double:
                    f.write("{:< 20.6f} ".format(j))
                if type(j) is int:
                    f.write("{:< 20d} ".format(j))
            f.write("\n")

        # Write dihedral interactions
        # Note: Dihedral_type_dict was initialized by looping over sorted(set(Dihedral_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Dihedral interactions"))
        for i in sorted(set(Dihedral_types)):
            
                # The potential for the scanned dihedrals are set to zero
                if i in tables:
                    f.write("     {:20s} {:<10d} {:20s} {:20s} {:20s}".format("dihedral_coeff",Dihedral_type_dict[i],"table",tables[i][0],tables[i][1]))

                # This is clunky. The harmonic/opls type parameters reside in separate dictionaries so a check needs to be performed as to which. *** NOTE TRUE ANYMORE, AFTER TESTING THE LATER ELIF CAN BE REMOVED ***
                elif i in list(Dihedral_params.keys()):
                    f.write("     {:20s} {:<10d} ".format("dihedral_coeff",Dihedral_type_dict[i]))

                    # This check is performed in case a V0 normalization constant is included in the dihedral definition
                    if len(Dihedral_params[i]) == 6:
                        tmp = Dihedral_params[i][0] + Dihedral_params[i][2:]
                    else:
                        tmp = Dihedral_params[i]

                    # Write the parameters (with trimmed normalization constant)
                    for j in tmp:
                        if type(j) is str:
                            f.write("{:20s} ".format(j))
                        if type(j) is float:
                            f.write("{:< 20.6f} ".format(j))
                        if type(j) is double:
                            f.write("{:< 20.6f} ".format(j))
                        if type(j) is int:
                            f.write("{:< 20d} ".format(j))

                # Write the dihedral params
                elif i in list(Dihedral_harmonic_params.keys()):

                    f.write("     {:20s} {:<10d} ".format("dihedral_coeff",Dihedral_type_dict[i]))
                    for j in Dihedral_harmonic_params[i]:
                        if type(j) is str:
                            f.write("{:20s} ".format(j))
                        if type(j) is float:
                            f.write("{:< 20.6f} ".format(j))
                        if type(j) is double:
                            f.write("{:< 20.6f} ".format(j))
                        if type(j) is int:
                            f.write("{:< 20d} ".format(j))
                    if Dihedral_harmonic_params[i][0] == "harmonic":
                        f.write(" {:<10d}".format(2))   # Harmonic dihedrals in LAMMPS have an optional "n" parameter that multiplies the cos argument. Fits in taffi all assume 2.
                    
                # Print error and exit if the dihedral parameters are not available
                else:
                    print("ERROR: Dihedral parameters for {} could not be found. Exiting...".format(i))
                    quit()

                f.write("\n")

    return Atom_type_dict

# Description: A wrapper for the commands to write a lammps input file to optimize a geometry
def write_minimize_input(Filename,frequency,Onefourscale,Pair_styles,Dihedral_styles,fixes=[]):

    with open(Filename+'.in.init','w') as f:
        f.write("# lammps input file for polymer simulation with dilute ions\n\n"+\
                "# VARIABLES\n"+\
                "variable        data_name       index   {}\n".format(Filename+'.data')+\
                "variable        settings_name   index   {}\n".format(Filename+'.in.settings')+\
                "variable        coords_freq     index   {}\n".format(frequency)+\
                "variable        vseed           index   {: <6d}\n".format(int(random.random()*100000))+\
                
                "#===========================================================\n"+\
                "# GENERAL PROCEDURES\n"+\
                "#===========================================================\n"+\
                "units		real	# g/mol, angstroms, fs, kcal/mol, K, atm, charge*angstrom\n"+\
                "dimension	3	# 3 dimensional simulation\n"+\
                "newton		off	# use Newton's 3rd law\n"+\
                "boundary	f f f	# fixed boundary conditions \n"+\
                "atom_style	full    # molecular + charge\n\n"+\

                "#===========================================================\n"+\
                "# FORCE FIELD DEFINITION\n"+\
                "#===========================================================\n"+\
                'special_bonds  lj/coul  0.0 0.0 0.0     # NO 1-4 LJ/Coul interactions\n'+\
                'pair_style     hybrid {} # outer_LJ outer_Coul (cutoff values, see LAMMPS Doc)\n'.format(Pair_styles)+\
                'pair_modify    shift yes\n'+\
                'bond_style     harmonic             # parameters needed: k_bond, r0\n'+\
                'angle_style    harmonic             # parameters needed: k_theta, theta0\n')
        if len(Dihedral_styles) >= 1:
            f.write('dihedral_style hybrid {}            # parameters depend on the kinds of dihedrals that are present\n'.format(' '.join(Dihedral_styles)))                    
        f.write('pair_modify    shift yes mix arithmetic table 0       # using Lorenz-Berthelot mixing rules\n\n'+\

                "#===========================================================\n"+\
                "# SETUP SIMULATIONS\n"+\
                "#===========================================================\n\n"+\
                "# READ IN COEFFICIENTS/COORDINATES/TOPOLOGY\n"+\
                'read_data ${data_name}\n'+\
                'include ${settings_name}\n\n'+\

                "# SET RUN PARAMETERS\n"+\
                "timestep	1.0		# fs\n"+\
                "run_style	verlet 		# Velocity-Verlet integrator\n"+\
                "neigh_modify every 1 delay 0 check no # More relaxed rebuild criteria can be used\n\n"+\
        
                "# Define Theromstyle variables for print\n"+\
                "thermo_style    custom step temp vol density etotal pe ebond eangle edihed ecoul elong evdwl\n"+\
                "variable        my_pe   equal   pe\n"+\
                "variable        my_ebon equal   ebond\n"+\
                "variable        my_eang equal   eangle\n"+\
                "variable        my_edih equal   edihed\n"+\
                "variable        my_evdw equal   evdwl\n"+\
                "variable        my_ecol equal   ecoul\n\n")

        # Write fixes if they are supplied
        if fixes != []:
            f.write("{}".format(fixes))

        f.write("#===========================================================\n"+\
                "# RUN MINIMIZATION\n"+\
                "#===========================================================\n\n")

        f.write('# Set minimization style (hftn seems to be the most robust)\n')
        f.write('min_style hftn\n\n')

        # Run minimization on the initial structure
        max_evals = 100000
        f.write("# minimize molecule energy with restraints\n"+\
                "dump           minimization all atom ${coords_freq}"+" minimize.lammpstrj\n"+\
                "dump_modify    minimization scale no\n"+\
                "dump_modify    minimization sort  id\n"+\
                "minimize 0.0 0.0 {} 100000\n".format(max_evals)+\
                "undump minimization\n")
    return

# # A wrapper for the commands to parse the bonds, angles, and dihedrals from the adjacency matrix.
# # Returns:   list of atomtypes, bond_types, bond instances, angle_types, angle instances, dihedral_types,
# #            diehdral instances, charges, and VDW parameters.
# def Find_parameters(Adj_mat,Geometry,Atom_types,FF_db="FF_file"):

#     # Initialize lists of each instance and type of FF object.
#     # instances are stored as tuples of the atoms involved 
#     # (e.g., bonds between atoms 1 and 13 and 17 and 5 would be stored as [(1,13),(17,5)] 
#     # Similarly, types are stored as tuples of atom types.
#     Bonds = []
#     Bond_types = []
#     Angles = []
#     Angle_types = []
#     Dihedrals = []
#     Dihedral_types = []

#     # Find bonds #
#     print "Parsing bonds..."
#     for count_i,i in enumerate(Adj_mat):        
#         Tmp_Bonds = [ (count_i,count_j) for count_j,j in enumerate(i) if j == 1 and count_j > count_i ]

#         # Store bond tuple so that lowest atom *type* between the first and the second atom is placed first
#         # and avoid redundant placements
#         for j in Tmp_Bonds:
#             if Atom_types[j[1]] < Atom_types[j[0]] and (j[1],j[0]) not in Bonds and (j[0],j[1]) not in Bonds:
#                 Bonds = Bonds + [ (j[1],j[0]) ]
#                 Bond_types = Bond_types + [ (Atom_types[j[1]].split('-UA')[0],Atom_types[j[0]].split('-UA')[0]) ]
#             elif (j[0],j[1]) not in Bonds and (j[1],j[0]) not in Bonds:
#                 Bonds = Bonds + [ (j[0],j[1]) ]
#                 Bond_types = Bond_types + [ (Atom_types[j[0]].split('-UA')[0],Atom_types[j[1]].split('-UA')[0]) ]


#     # Find angles #
#     print "Parsing angles..."
#     for i in Bonds:        

#         # Find angles based on connections to first index of Bonds
#         Tmp_Angles = [ (count_j,i[0],i[1]) for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j != i[1] ]

#         # Store angle tuple so that lowest atom *type* between the first and the third is placed first
#         # and avoid redundant placements
#         for j in Tmp_Angles:
#             if Atom_types[j[2]] < Atom_types[j[0]] and (j[2],j[1],j[0]) not in Angles and (j[0],j[1],j[2]) not in Angles:
#                 Angles = Angles + [(j[2],j[1],j[0])]
#                 Angle_types = Angle_types + [ (Atom_types[j[2]].split('-UA')[0],Atom_types[j[1]].split('-UA')[0],Atom_types[j[0]].split('-UA')[0]) ]
#             elif (j[0],j[1],j[2]) not in Angles and (j[2],j[1],j[0]) not in Angles:
#                 Angles = Angles + [(j[0],j[1],j[2])]
#                 Angle_types = Angle_types + [ (Atom_types[j[0]].split('-UA')[0],Atom_types[j[1]].split('-UA')[0],Atom_types[j[2]].split('-UA')[0]) ]

#         # Find angles based on connections to second index of Bonds
#         Tmp_Angles = [ (i[0],i[1],count_j) for count_j,j in enumerate(Adj_mat[i[1]]) if j == 1 and count_j != i[0] ]

#         # Store angle tuple so that lowest atom *type* between the first and the third is placed first
#         # and avoid redundant placements
#         for j in Tmp_Angles:
#             if Atom_types[j[2]] < Atom_types[j[0]] and (j[2],j[1],j[0]) not in Angles and (j[0],j[1],j[2]) not in Angles:
#                 Angles = Angles + [(j[2],j[1],j[0])]
#                 Angle_types = Angle_types + [ (Atom_types[j[2]].split('-UA')[0],Atom_types[j[1]].split('-UA')[0],Atom_types[j[0]].split('-UA')[0]) ]
#             elif (j[0],j[1],j[2]) not in Angles and (j[2],j[1],j[0]) not in Angles:
#                 Angles = Angles + [(j[0],j[1],j[2])]
#                 Angle_types = Angle_types + [ (Atom_types[j[0]].split('-UA')[0],Atom_types[j[1]].split('-UA')[0],Atom_types[j[2]].split('-UA')[0]) ]

        
#     # Find dihedrals #
#     print "Parsing dihedrals..."
#     for i in Angles:
        
#         # Find atoms attached to first atom of each angle
#         Tmp_Dihedrals = [ (count_j,i[0],i[1],i[2]) for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in [i[1],i[2]] ]
        
#         # Store dihedral tuple so that the lowest atom *type* between the first and fourth is placed first
#         # and avoid redundant placements        
#         for j in Tmp_Dihedrals:

#             # If the first and fourth atoms are equal, then sorting is based on the second and third
#             if Atom_types[j[3]] == Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
#                 if Atom_types[j[2]] < Atom_types[j[1]]:
#                     Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
#                     Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
#                 else:
#                     Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
#                     Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

#             elif Atom_types[j[3]] < Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
#                 Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
#                 Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
#             elif (j[0],j[1],j[2],j[3]) not in Dihedrals and (j[3],j[2],j[1],j[0]) not in Dihedrals:
#                 Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
#                 Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

#         # Find atoms attached to the third atom of each angle
#         Tmp_Dihedrals = [ (i[0],i[1],i[2],count_j) for count_j,j in enumerate(Adj_mat[i[2]]) if j == 1 and count_j not in [i[0],i[1]] ]
        
#         # Store dihedral tuple so that the lowest atom *type* between the first and fourth is placed first
#         # and avoid redundant placements        
#         for j in Tmp_Dihedrals:

#             # If the first and fourth atoms are equal, then sorting is based on the second and third
#             if Atom_types[j[3]] == Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
#                 if Atom_types[j[2]] < Atom_types[j[1]]:
#                     Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
#                     Dihedral_types = Dihedral_types + [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
#                 else:
#                     Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
#                     Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]

#             elif Atom_types[j[3]] < Atom_types[j[0]] and (j[3],j[2],j[1],j[0]) not in Dihedrals and (j[0],j[1],j[2],j[3]) not in Dihedrals:
#                 Dihedrals = Dihedrals + [(j[3],j[2],j[1],j[0])]
#                 Dihedral_types = Dihedral_types+ [ (Atom_types[j[3]],Atom_types[j[2]],Atom_types[j[1]],Atom_types[j[0]]) ]
#             elif (j[0],j[1],j[2],j[3]) not in Dihedrals and (j[3],j[2],j[1],j[0]) not in Dihedrals:
#                 Dihedrals = Dihedrals + [(j[0],j[1],j[2],j[3])]
#                 Dihedral_types = Dihedral_types + [ (Atom_types[j[0]],Atom_types[j[1]],Atom_types[j[2]],Atom_types[j[3]]) ]
                    
#     ##############################################################
#     # Read in parameters: Here the stretching, bending, dihedral #
#     # and non-bonding interaction parameters are read in from    #
#     # parameters file. Mass and charge data is also included.    #
#     # The program looks for a simple match for the first entry   #
#     # in each line with one of the bond or angle types.          #
#     # INPUT: param_file, BOND_TYPES_LIST, ANGLE_TYPES_LIST       #
#     #        DIHEDRAL_TYPES_LIST, ELEMENTS                       #
#     # OUTPUT: CHARGES, MASSES, BOND_PARAMS, ANGLE_PARAMS,        #
#     #         DIHERAL_PARAMS, PW_PARAMS                          #
#     ##############################################################

#     # Initialize dictionaries

#     # Read in masses and charges
#     Masses = {}
#     with open(FF_db,'r') as f:
#         content=f.readlines()
        
#     for lines in content:
#         fields=lines.split()

#         # Skip empty lines
#         if len(fields) == 0:
#             continue

#         if fields[0].lower() == 'atom':
#             Masses[fields[1]] = float(fields[3]) 

#     # Read in bond parameters
#     Bond_params = {}
#     with open(FF_db,'r') as f:
#         content=f.readlines()

#     for lines in content:
#         fields=lines.split()

#         # Skip empty lines
#         if len(fields) == 0:
#             continue

#         if fields[0] == 'bond':
#             if fields[3] == "harmonic":
#                 Bond_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
#             else:
#                 print "ERROR: only harmonic bond definitions are currently supported by gen_md_for_vdw.py. Exiting..."
#                 quit()

#     # Read in angle parameters
#     Angle_params = {}
#     with open(FF_db,'r') as f:
#         content=f.readlines()

#     for lines in content:
#         fields=lines.split()

#         # Skip empty lines
#         if len(fields) == 0:
#             continue

#         if fields[0].lower() == 'angle':
#             if fields[4] == "harmonic":
#                 Angle_params[(fields[1],fields[2],fields[3])] = [fields[4],float(fields[5]),float(fields[6])]
#             else:
#                 print "ERROR: only harmonic angle definitions are currently supported by gen_md_for_vdw.py. Exiting..."
#                 quit()

#     # Read in dihedral parameters
#     Dihedral_params = {}
#     with open(FF_db,'r') as f:
#         content=f.readlines()

#     for lines in content:
#         fields=lines.split()

#         # Skip empty lines
#         if len(fields) == 0:
#             continue

#         if fields[0].lower() in ['dihedral','torsion']:            
#             if fields[5] == "opls":
#                 Dihedral_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
#             elif fields[5] == "harmonic":
#                 Dihedral_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ]
#             else:
#                 print "ERROR: Only opls and harmonic type dihedral definitions are currently supported by gen_md_for_vdw.py. Exiting..."
#                 quit()
                

#     # Search for charges based on atom type
#     with open(FF_db,'r') as f:
#         content=f.readlines()

#     Charges = zeros(len(Atom_types))
#     for i in range(len(Atom_types)):
#         for lines in content:
#             fields=lines.split()

#             # Skip empty lines
#             if len(fields) == 0:
#                 continue
                    
#             if fields[0].lower() in ['charge'] and Atom_types[i] == fields[1]:
#                 Charges[i] = float(fields[2])

#     # Search for VDW parameters
#     VDW_params = {}
#     with open(FF_db,'r') as f:
#         for lines in f:
#             fields = lines.split()

#             # Skip empty lines
#             if len(fields) == 0:
#                 continue
                    
#             if fields[0].lower() in ['vdw']:

#                 if fields[1] > fields[2]:
#                     VDW_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
#                 else:
#                     VDW_params[(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]

#     # Check for missing parameters
#     Missing_masses = [ i for i in Atom_types if str(i) not in Masses.keys() ] 
#     Missing_charges = [ count_i for count_i,i in enumerate(Charges) if i == -100.0 ]; Missing_charges = [ Atom_types[i] for i in Missing_charges ]
#     Missing_bonds = [ i for i in Bond_types if (i[0],i[1]) not in Bond_params.keys() ]
#     Missing_angles = [ i for i in Angle_types if (i[0],i[1],i[2]) not in Angle_params.keys() ]
#     Missing_dihedrals = [ i for i in Dihedral_types if (i[0],i[1],i[2],i[3]) not in Dihedral_params.keys() ]

#     # Print diagnostics on missing parameters and quit if the prerequisites are missing.
#     if ( len(Missing_masses) + len(Missing_charges) + len(Missing_bonds) + len(Missing_angles) + len(Missing_dihedrals) ) > 0:
#         print "\nUh Oh! There are missing FF parameters...\n"

#         if Missing_masses:
#             print "Missing masses for the following atom types: {}".format([ i for i in set(Missing_masses) ])
#         if Missing_charges:
#             print "Missing charges for the following atom types: {}".format([ i for i in set(Missing_charges) ])
#         if Missing_bonds:
#             print "Missing bond parameters for the following bond types: {}".format([ i for i in set(Missing_bonds) ])
#         if Missing_angles:
#             print "Missing angle parameters for the following angle types: {}".format([ i for i in set(Missing_angles) ])
#         if Missing_dihedrals:
#             print "Missing dihedral parameters for the following dihedral types: {}".format([ i for i in set(Missing_dihedrals) ])
        
#         print "\nEnsure the specification of the missing parameters. Exiting..."
#         quit()

#     return Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,Dihedrals,Dihedral_types,Dihedral_params,Charges.tolist(),Masses,VDW_params

# # Generates the adjacency matrix 
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
        
# def xyz_parse(input):

#     # Open file and read contents into variable
#     with open(input,'r') as f:
#             content=f.readlines()

#     # Find number of atoms and initialize various matrices
#     Atom_Number = int(content[0].split()[0])
#     Elements = ['']*Atom_Number
#     Geometry = zeros([Atom_Number,3])
#     Atom_types = [0]*Atom_Number

#     # Iterate over the remainder of contents and read the
#     # geometry and elements into variable. Note that the
#     # first two lines are considered a header
#     count=0
#     for lines in content[2:]:
#         fields=lines.split()

#         # Skip empty lines
#         if len(fields) == 0:
#             continue

#         # Write geometry containing lines to variable
#         if len(fields) > 3:
#             Elements[count]=fields[0]
#             Geometry[count,:]=array([float(fields[1]),float(fields[2]),float(fields[3])])

#             # If atom type information is stored in the fifth column it is saved to variable
#             if len(fields) >= 5:
#                 Atom_types[count] = fields[4]

#             # If the nuclear charge is stored in the sixth column it is saved to variable
#             if len(fields) >= 6:
#                 Charges[count]=float(fields[5])
#             count = count + 1

#     return Elements,Geometry,Atom_types

# # Description: finds the number of disconnected subnetworks in the 
# #              adjacency matrix, which corresponds to the number of 
# #              separate molecules.
# #
# # Inputs:      adj_mat: numpy array holding a 1 in the indices of bonded
# #                        atom types. 
# #
# # Returns:     mol_count: scalar, the number of molecules in the adj_mat
# def mol_count(adj_mat):
    
#     # Initialize list of atoms assigned to molecules and a counter for molecules
#     placed_idx = []    
#     mol_count = 0

#     # Continue the search until all the atoms have been assigned to molecules
#     while len(placed_idx)<len(adj_mat):

#         # Use sequential elements of the adj_mat as seeds for the spanning network search
#         for count_i,i in enumerate(adj_mat):

#             # Only proceed with search if the current atom hasn't been placed in a molecule
#             if count_i not in placed_idx:

#                 # Increment mol_count for every new seed and add the seed to the list of placed atoms
#                 mol_count += 1               
#                 placed_idx += [count_i]
                
#                 # Find connections
#                 idx = [ count_j for count_j,j in enumerate(i) if j==1 and count_j not in placed_idx ]
                
#                 # Continue until no new atoms are found
#                 while len(idx) > 0:
#                     current = idx.pop(0)
#                     if current not in placed_idx:
#                         placed_idx += [current]
#                         idx += [ count_k for count_k,k in enumerate(adj_mat[current]) if k == 1 and count_k not in placed_idx ]
#     return mol_count

# # Description:
# # Rotate Point by an angle, theta, about the vector with an orientation of v1 passing through v2. 
# # Performs counter-clockwise rotations (i.e., if the direction vector were pointing
# # at the spectator, the rotations would appear counter-clockwise)
# # For example, a 90 degree rotation of a 0,0,1 about the canonical 
# # y-axis results in 1,0,0.
# #
# # Point: 1x3 array, coordinates to be rotated
# # v1: 1x3 array, rotation direction vector
# # v2: 1x3 array, point the rotation vector passes through
# # theta: scalar, magnitude of the rotation (defined by default in degrees)
# def axis_rot(Point,v1,v2,theta,mode='angle'):

#     # Temporary variable for performing the transformation
#     rotated=array([Point[0],Point[1],Point[2]])

#     # If mode is set to 'angle' then theta needs to be converted to radians to be compatible with the
#     # definition of the rotation vectors
#     if mode == 'angle':
#         theta = theta*pi/180.0

#     # Rotation carried out using formulae defined here (11/22/13) http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/)
#     # Adapted for the assumption that v1 is the direction vector and v2 is a point that v1 passes through
#     a = v2[0]
#     b = v2[1]
#     c = v2[2]
#     u = v1[0]
#     v = v1[1]
#     w = v1[2]
#     L = u**2 + v**2 + w**2

#     # Rotate Point
#     x=rotated[0]
#     y=rotated[1]
#     z=rotated[2]

#     # x-transformation
#     rotated[0] = ( a * ( v**2 + w**2 ) - u*(b*v + c*w - u*x - v*y - w*z) )\
#              * ( 1.0 - cos(theta) ) + L*x*cos(theta) + L**(0.5)*( -c*v + b*w - w*y + v*z )*sin(theta)

#     # y-transformation
#     rotated[1] = ( b * ( u**2 + w**2 ) - v*(a*u + c*w - u*x - v*y - w*z) )\
#              * ( 1.0 - cos(theta) ) + L*y*cos(theta) + L**(0.5)*(  c*u - a*w + w*x - u*z )*sin(theta)

#     # z-transformation
#     rotated[2] = ( c * ( u**2 + v**2 ) - w*(a*u + b*v - u*x - v*y - w*z) )\
#              * ( 1.0 - cos(theta) ) + L*z*cos(theta) + L**(0.5)*( -b*u + a*v - v*x + u*y )*sin(theta)

#     rotated = rotated/L
#     return rotated

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
def initialize_VDW(atomtypes,sigma_scale=1.0,eps_scale=1.0,VDW_type="lj/cut/coul/long",VDW_FF={}):

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
            if (i,j) in VDW_FF:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
            elif (j,i) in VDW_FF:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF:
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
    print("\n{}".format("*"*167))
    print("* {:^163s} *".format("Initializing VDW parameters for the simulation (those with * were read from the FF file(s))"))
    print("*{}*".format("-"*165))
    print("* {:<50s} {:<50s} {:<20s}  {:<18s} {:<18s}   *".format("Type","Type","VDW_type","eps (kcal/mol)","sigma (angstroms)"))
    print("{}".format("*"*167))
    for j in list(VDW_dict.keys()):
        if j in VDW_FF:
            print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f} *".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2]))
        else:
            print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f}".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2]))
    print("")


    return VDW_dict

# # Description:   A wrapper for the commands to create a data dictionary holding the FF, mode, and geometry information for each molecule
# # Inputs:  Filename:
# # Returns: Data: a dictionary keyed to the list of coord_files
# #                a keyed entry is itself a dictionary holding the geometry, 
# #                elements, charges, bonds, bond_types, angles, angle_types,
# #                dihedrals, dihedral_types and masses.
# #                subdictionary keys (e.g., Data[*]["Geometry"])
# #                Geometry: array of the molecule's geometry
# #                Elements: list of the molecule's elements (indexed to the rows of "Geometry")
# #                Atom_types: list of the molecule's atom_types (indexed to the rows of "Geometry")
# #                Charges: list of the molecule's charges (indexed to the rows of "Geometry")
# #                Bonds: a list of the molecule's bonds, each element being a tuple of the atom indices that are bonded
# #                Angles: a list of the molecule's angles, each element being a tuple of the atom indices that form an angle mode
# #                Dihedrals: a list of the molecule's dihedrals, each element being a tuple of the atom indices that form a dihedral
# #                Bond_types: a list of the molecule's bond_types, indexed to Bonds
# #                Angle_types: a list of the molecule's angle_types, indexed to Angles
# #                Dihedral_types: a list of the molecule's dihedral_types, indexed to Dihedrals
# #
# def get_data(FF_all,coord_files):

#     # Initialize dictionary to hold the FF, mode, and geometry information about each molecule
#     Data = {}

#     # Iterate over all molecules being simulated and collect their FF, mode, and geometry information.
#     for count_i,i in enumerate(coord_files):

#         # Initialize dictionary for this geometry and set number of molecules to place in the simulated system
#         Data[i] = {}      

#         # Extract Element list and Coord list from the file
#         Data[i]["Elements"],Data[i]["Geometry"] = xyz_parse(i)

#         # Generate adjacency table
#         Data[i]["Adj_mat"] = Table_generator(Data[i]["Elements"],Data[i]["Geometry"])

#         # Check the number of molecules
#         mol_in_in = mol_count(Data[i]["Adj_mat"])
#         if mol_in_in > 1:
#             print "ERROR: {} molecules were discovered in geometry {}. Check the geometry of the input file. Exiting...".format(mol_in_in,i)
#             quit()

#         # Parse the atomtypes
#         Data[i]["Atom_types"] = id_types(Data[i]["Elements"],Data[i]["Adj_mat"],geo=Data[i]["Geometry"])

#         # Generate list of bonds angles and dihedrals    
#         print "\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Parsing Modes and FF Information for Molecule {}".format(i),"*"*167)
#         Data[i]["Bonds"],Data[i]["Bond_types"],Data[i]["Bond_params"],Data[i]["Angles"],Data[i]["Angle_types"],Data[i]["Angle_params"],\
#         Data[i]["Dihedrals"],Data[i]["Dihedral_types"],Data[i]["Dihedral_params"],Data[i]["Charges"],Data[i]["Masses"],Data[i]["VDW_params"] =\
#             Find_parameters(Data[i]["Adj_mat"],Data[i]["Geometry"],Data[i]["Atom_types"],FF_db=FF_all)

#         # Print System characteristics
#         print "\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Mode Summary for Molecule {}".format(i),"*"*167)
#         print "\nAtom_types ({}):\n".format(len(set(Data[i]["Atom_types"])))
#         for j in sorted(set(Data[i]["Atom_types"])):
#             print "\t{}".format(j)
#         print "\nBond types ({}):\n".format(len(set(Data[i]["Bond_types"])))
#         for j in sorted(set(Data[i]["Bond_types"])):
#             print "\t{}".format(j)
#         print "\nAngle types ({}):\n".format(len(set(Data[i]["Angle_types"])))
#         for j in sorted(set(Data[i]["Angle_types"])):
#             print "\t{}".format(j)
#         print "\nDihedral types ({}):\n".format(len(set(Data[i]["Dihedral_types"])))
#         for j in sorted(set(Data[i]["Dihedral_types"])):

#             print "\t{}".format(j)

#         print "\n{}".format("*"*167)
#         print "* {:^163s} *".format("System Characteristics")
#         print "*{}*".format("-"*165)
#         print "* {:<87s} {:<24s} {:<24s} {:<25s} *".format("Type","Element","Mass","Charge")
#         print "{}".format("*"*167)
#         for j in range(len(Data[i]["Atom_types"])):
#             print " {:<88s} {:<23s} {:< 24.6f} {:< 24.6f}".format(Data[i]["Atom_types"][j],Data[i]["Elements"][j],Data[i]["Masses"][str(Data[i]["Atom_types"][j])],Data[i]["Charges"][j])

#     return Data

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

# Logger class for redirecting stdout to a logfile
class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/lammps_minimize_molecule.log", "a",buffering = 1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass

# main sentinel
if __name__ == "__main__":
   main(sys.argv[1:])
