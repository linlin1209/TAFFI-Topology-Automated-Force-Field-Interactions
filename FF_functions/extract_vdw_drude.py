#!/bin/env python                                                                                                                                                              
import sys,argparse,os,fnmatch,random,time,matplotlib
import numpy as np
from numpy.linalg import matrix_rank
from scipy import stats
from scipy.spatial.distance import cdist
from scipy.optimize import minimize
from copy import deepcopy
matplotlib.use('Agg') # Needed for cluster image generation
from matplotlib import pyplot as plt
# Add TAFFI Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from adjacency import Table_generator
from id_types import id_types


def main(argv):

    parser = argparse.ArgumentParser(description='This program collects force field parameters from the output of Orca calculations. '+\
                                                 'This script should be used in conjunction with the gen_jobs_for_vdw.py program, which automatically generates '\
                                                 'the Orca jobs for submission. After the quantum chemistry completes, this script collates the information and '\
                                                 ' generates a database file that can be used in conjunction with polygen to generate the MD job.')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-f', dest='base_name', default='*',
                        help = 'The program operates on all files discovered during a directory walk from the working directory whose name matches this variable.')

    parser.add_argument('-o', dest='output_folder', default='vdw_fits',
                        help = 'All outputs from the fitting are saved to a subfolder of the working directory (determined by the -f argument). '+\
                               'This option determines the folder name. (default: vdw_fits)') 

    parser.add_argument('-gens', dest='gens',type=int, default=2,
                        help = 'When atomtypes are being automatically determined (i.e. when the supplied *.xyz files do not have atomtype information) this variable controls the bond depth used to define unique atomtypes. (default: 2)')

    parser.add_argument('-FF_alpha', dest='FF_alpha', default='',
                        help = 'db file that contains alphas, only alphas in this file will be used, other parameters are obsolete')
    
    parser.add_argument('-c', dest='cycles', default='all',
                        help = 'This variable determines which configurations are included in the fit set. The configs folder which holds the sampled configurations can possess cycle labels ( i.e. "cycle_num-config_num" ). The '+\
                               'user can include all configurations by supplying "all" (the default), or supply individual cycle labels as a space-delimited string (i.e. "1 3 5") or use a colon syntax for the start:end of the '+\
                               'cycles to include (i.e. "1:4" would included cycles 1 through 4, inclusive). The last two options can be used in combination. (default: "all")')
                        
    parser.add_argument('-fit', dest='fit', default="lj",
                        help = 'Specifies the function used to fit the VDW potential (options: buck, lj, born ; default: lj)')
    
    parser.add_argument('-max_cycles', dest='max_cycles',type=int,default=10000,
                        help = 'Specifies the maximum number of iterations to perform before terminating the fit. (default: 1000)')

    parser.add_argument('-min_cycles', dest='min_cycles',type=int, default=10,
                        help = 'Specifies the minimum number of iterations to perform before terminating the fit. (default: 100)')

    parser.add_argument('-xhi2_thresh', dest='xhi2_thresh',type=float, default=1E-6,
                        help = 'Specifies the termination threshold for differences in xhi2 between cycles. '+\
                               'Only breaks once min_cycles of iterations have been performed (kcal2/mol2; default: 1E-6)')

    parser.add_argument('-weight', dest='weight',type=float, default=1.0,
                        help = 'Specifies the weighting factor for updating the VDW parameters at each iteration. For example, '+\
                               '1.0 means the old value is completely placed by the new value, and 0.5 means the updated value '+\
                               'is the average of the new and old. (default: 1.0)')

    parser.add_argument('-N_config', dest='N_config',type=int, default=0,
                        help = 'Specifies the number of configurations to use for the fit. By default (0) all discovered configurations '+\
                               'are used. When set to a value a random subset of the discovered configurations is used. This option is '+\
                               'intended to be used for convergence tests with respect to the number of configurations. If N_config is set '+\
                               'to a value larger than the number of discovered configurations, then a warning is printed and the program '+\
                               'defaults to using all configurations. (default: 0, all configurations)') 

    parser.add_argument('-FF_DFT', dest='DFT_db', default=[],
                        help = 'This variable holds the filename of a DFT parameters database. When supplied, the program will avoid parametrizing any '+\
                               'vdw pairs at the DFT level that already exist in the database (default: None)')

    parser.add_argument('-FF_MP2', dest='MP2_db', default=[],
                        help = 'This variable holds the filename of a MP2 parameters database. When supplied, the program will avoid parametrizing any '+\
                               'vdw pairs at the MP2 level that already exist in the database (default: None)')

    parser.add_argument('-E_max', dest='E_max',type=float, default=10.0,
                        help = 'This variable holds the maximum energy threshold for using a configuration in the fitting procedure. For example, settings '+\
                               'E_max to 10.0 will only use configurations whose interaction energy is less than 10.0 kcal/mol in the fit. (default: 10.0)')

    parser.add_argument('-E_min', dest='E_min',type=float, default=-10000.0,
                        help = 'This variable holds the minimum energy threshold for using a configuration in the fitting procedure. For example, settings '+\
                               'E_max to 0.0 will only use configurations whose interaction energy is greater than 0.0 kcal/mol in the fit. This isn\'t usually '+\
                               'useful since the low energy portion of the interaction potential is most relevant to the fit; it is included for completeness (default: -10000.0)')

    parser.add_argument('-seed', dest='seed',type=int, default=444,
                        help = 'The seed for the random number generator. (default: 444)')

    parser.add_argument('-q_a', dest='q_a', default=None,
                        help = 'The molecular charge on molecule a. By default the charge is automatically determined by rounding the total charge to the nearest integer. (default: None)')

    parser.add_argument('-q_b', dest='q_b', default=None,
                        help = 'The molecular charge on molecule b. By default the charge is automatically determined by rounding the total charge to the nearest integer. (default: None)')

    parser.add_argument('-r_min_scale', dest='r_min_scale',type=float, default=1.0,
                        help = 'Scales the r_min value for the VDW fit. By default, the minima of fit potetials are retrained to lie within the sampled region of pair separations. '+\
                               'This behavior can be modified by setting this variable to <1 (more relaxed) or >1 (forces larger r_min). (default: 1.0)')

    parser.add_argument('-UA_vol_scale', dest='UA_vol_scale',type=float, default=1.0,
                        help = 'Scale factor for UA multiplying the interaction volume of the UA pseudo-atoms. (default: 1.0)')

    parser.add_argument('-mixing_rule', dest='mixing_rule', default="None",
                        help = 'Mixing rule for heteroatomic lennard-jones interactions. Valid options are none (no mixing rules imposed), lb (Lorentz-Berthelot), and wh (Waldman-Hagler). The options are '+\
                               'case insensitive (default: None)')

    parser.add_argument('-L2_sigma', dest='L2_s',type=float, default=0.1,
                        help = 'Sets the weight of L2 regularization for sigma values in the fit penalty function. When set to zero, the LJ parameters are fit using least-squares criteria, '+\
                               'non-zero values result in L2 regularization of the sigma values (default: 0.1)')

    parser.add_argument('-L2_eps', dest='L2_e',type=float, default=0.01,
                        help = 'Sets the weight of L2 regularization for the epsilon values in the fit penalty function. When set to zero, the LJ parameters are fit using least-squares criteria, '+\
                               'non-zero values result in L2 regularization of the epsilon values (default: 0.01)')                               

    parser.add_argument('-UA_method', dest='UA_method', default='AA_fit',
                        help = 'Determines the method for calculating the UA values. '+\
                               '"radius" uses the volume of the sphere that encompasses all UA atoms in the group, scaled by a factor determined by the -UA_vol_scale variable. '+\
                               '"volume" uses the monte-carlo calculated volume of the UA group, with overlaps naturally subtracted. '+\
                               'Both, radius and volume methods base the UA eps values upon scaled AA-eps values in order to reproduce the intermolecular '+\
                               'LJ interaction energy. "fit" performs an iterative fit of the eps and sigma values to reproduce the DFT data, independent of the AA '+\
                               'parameters. radius and volume are meant for use fitting self interactions (where A and B molecules are identical), the fit method is meant '+\
                               'for fitting heteromolecular interactions (e.g., ion polymer or solute solvent) (default: radius)')

    parser.add_argument('-tag',dest='tag',default='',
                        help = 'tag for FF folder and relax drude, this assumes FF name is FF_tag and output drude: _drude_tag.xyz, default:None')

    parser.add_argument('-fun', dest='functional', default='B3LYP',
                        help = 'functional for get_data to identify job_type. (default: B3LYP; other typical options are wB97X-D3, and M062X)')

    parser.add_argument('--extract_charges', dest='extract_charges', default=0, const=1, action='store_const',
                        help = 'If this flag is enabled then the charges are parsed and averaged over all QC calculations (default: off)')

    parser.add_argument('--avoid_read', dest='avoid_read', default=0, const=1, action='store_const',
                        help = 'If this flag is enabled then the force-field potentials found in any present "initial_params.db" files are not used. (default: use initial_params.db to initialize LJ guess)')

    parser.add_argument('--keep_outlier', dest='outlier', default=True, const=False, action='store_const',
                        help = 'If this flag is enabled then the outliers will be kept, default: remove outlier (have modified z-score>3.5')

    parser.add_argument('--polar', dest='polar_flag', default=False, const=True, action='store_const',
                        help = 'If this flag is enabled will do Drude fitting, default: regular fitting')


    # Declare global variables and save the working directory to variable
    global Htokcalmol,Ecoul_Const,Data,min_pairs
    working_dir = os.getcwd()       # Not sure if this is useful anymore
    Htokcalmol = 627.509            # Converts hartree to kcal/mol
    Ecoul_Const = 332.0702108431797 # Converts elementary_charge^2/ang to kcal/mol 

    # Convert inputs to the proper data types
    args=parser.parse_args(argv)
    args.fit = str(args.fit).lower()
    args.mixing_rule = str(args.mixing_rule).lower()
    if args.DFT_db != []:
        args.DFT_db = args.DFT_db.split()
    if args.MP2_db != []:
        args.MP2_db = args.MP2_db.split()    

    # Consistency checks
    if args.fit not in [ 'lj', 'buck', 'born' ]: print("ERROR: an unrecognized fit type was requested by the user. Exiting..."); quit()
    if args.UA_method not in [ 'radius', 'volume', 'AA_fit', 'fit' ]: print("ERROR: an unrecognized -UA_method was requested by the user. Exiting..."); quit()
    if args.mixing_rule not in ['none','lb','wh']: print("ERROR: an unrecognized -mixing_rule was requested by the user. Exiting..."); quit()
    if args.min_cycles > args.max_cycles: print("ERROR: min_cycles must be less than or equal to max_cycles. Exiting..."); quit()
    if args.L2_s < 0.0: print("ERROR: -L2_sigma must be set to a float greater than 0. Exiting..."); quit()
    if args.L2_e < 0.0: print("ERROR: -L2_eps must be set to a float greater than 0. Exiting..."); quit()
    
    # Parse the cycles input. Additional complexity is due to the use of a colon minilanguage that allows the user to supply supply a range cycles.
    args.cycles = args.cycles.split()
    if "all" in args.cycles: args.cycles = ["*"]
    else:
        expanded = [ i for i in args.cycles if ":" in i ]
        args.cycles = [ i for i in args.cycles if i not in expanded ]   
        added = []
        # Expand the colon containing arguments
        for count_i,i in enumerate(expanded): added += [ str(j) for j in range(int(i.split(":")[0]),int(i.split(":")[1])+1) ]
        args.cycles += added
        args.cycles = [ i+'-*' for i in args.cycles ]

    # Use a fixed seed value to ensure reproducibility between runs
    random.seed(int(args.seed))

    # Check that the job directory exists
    if os.path.isdir(args.base_name) == False:
        print("ERROR: Requested directory could not be found. Please check the name. Exiting...")
        quit()
    else:
        os.chdir(args.base_name)
        if os.path.isdir(args.output_folder) == True:
            print("ERROR: Requested output folder ({}/{}) already exists. Exiting to avoid overwriting data...".format(args.base_name,args.output_folder))
            quit()
        else:
            os.makedirs(args.output_folder)
            os.makedirs(args.output_folder+'/figures')
            sys.stdout = Logger(args.output_folder)
            print("PROGRAM CALL: python extract_vdw.py {}\n".format(' '.join([ i for i in argv])))
        if args.DFT_db != []:
            args.DFT_db = [ '../'+i if i[0] != '/' else i for i in args.DFT_db ]
        if args.MP2_db != []:
            args.MP2_db = [ '../'+i if i[0] != '/' else i for i in args.MP2_db ]

    # Find all output files
    Base_Filename = '*out'
    Base_dir = os.getcwd()

    # Use fnmatch to use wildcards and os.walk to perform a directory walk. Names holds the path and filename of each 
    # file to be processed by the program. Name contains the unique thermo filenames that are parsed by the program. 
    # The program expects each trajectory being sampled to have the same number of thermo.avg files, indexed according
    # to the lambda step in the thermodynamic integration simulation. 
    Names = [os.path.join(dp,f) for dp, dn, filenames in os.walk('.') for f in filenames if fnmatch.fnmatch(f,Base_Filename) ]
    Names = keep_complete(Names)
    if len(Names) == 0: print("Sorry, no complete vdw jobs were discovered in the specified directory. Exiting..."); quit()
    Names = [ i for i in Names if True in [ fnmatch.fnmatch(i.split('/')[-1],j) for j in args.cycles ] ] # Only keep files that belong to cycles specifid by the args.cycles argument

    #################
    # PARSE QC Data #
    #################

    # Print header
    print("{}".format("*"*144))
    print("* {:^140s} *".format("Parsing QC Interaction Energies, Charges, and Configurations in {}".format(args.base_name)))
    print("{}".format("*"*144))

    # Returns a dictionary keyed to the filenames, with the information needed for the fits keyed
    # to subdictionaries. See the get_data function for more info on the keys that are populated.
    # NOTE: this is function parses both DFT and MP2 data. The subsequent fit function only makes use
    #       of a subset of the data depending on whether the MP2 or DFT information is being used. 
    Data, job_type = get_data(Names,args.functional,args.extract_charges,args.DFT_db,q_a=args.q_a,q_b=args.q_b,tag=args.tag)
   
    with open('fit_AA.txt','w') as f:
      for i in sorted([ Data[key]['e_dft_fit_AA_ori'] for key in Data]):
         f.write('{}\n'.format(i))
    with open('fit_EC.txt','w') as f:
      for i in sorted([E_C_Tot(i,'AA') for i in Data]):
         f.write('{}\n'.format(i))
    #quit()
    ########################
    # PARSE DFT PARAMETERS #
    ########################

    # If DFT data was discovered in the QC output files then the DFT vdw parameters are parsed
    if "dft" in job_type:
        print("in")

        # Print header
        print("\n{}".format("*"*144))
        print("* {:^140s} *".format("Parsing DFT-Interaction Energy Calculations for Configurations within {}".format(args.base_name)))
        print("{}".format("*"*144))

        # Remove entries based on E_min/E_max criteria
        for i in list(Data.keys()):
            if Data[i]["e_dft"] > args.E_max: Data.pop(i,None)
            elif Data[i]["e_dft"] < args.E_min: Data.pop(i,None)

        
        # Print diagnostic
        print("\n{:70s} {}".format("Number of Configurations Discovered:",len(Names)))    
        print("{:70s} {}".format("Number of Configurations Satisfying E_min/E_max criteria:",len(list(Data.keys()))))

        # Check that at least 1 configuration has been retained based on the min/max criteria
        if len(list(Data.keys())) == 0:
            print("\nERROR: No configurations were discovered that matched the E_max criteria ({:<2.2f} kcal/mol). Exiting...".format(args.E_max))
            quit()

        # Print a warning if more configurations were requested than were discovered.
        if args.N_config > len(list(Data.keys())):
            print("\n                                 ****************************************************************************")
            print("                                 * WARNING: Less configurations were discovered than requested by N_config. *")
            print("                                 *          Continuing with fit using all configurations.                   *")
            print("                                 ****************************************************************************\n")

        # Use the number of configurations requested by user (by default, args.N_configs=0, all are used)
        if args.N_config != 0 and args.N_config <= len(list(Data.keys())):
            keys = list(Data.keys())
            random.shuffle(keys)
            for i in keys[args.N_config:]:
                Data.pop(i,None)

        print("{:70s} {}".format("Number of Configurations Used for Fit:",len(list(Data.keys()))))

        # Parse the unique molecules types in the simulation
        parse_molecules(args.output_folder,Data)

        # Fit all possible pair-wise interactions
        Fits = Fit_pairs(args.output_folder,Data,args.fit,args.weight,args.xhi2_thresh,args.min_cycles,args.max_cycles,args.DFT_db,\
                         QC_type="DFT",charges_opt=args.extract_charges,UA_vol_scale=args.UA_vol_scale,UA_method=args.UA_method,r_min_scale=args.r_min_scale,\
                         avoid_read_flag=args.avoid_read,mixing_rule=args.mixing_rule,L2_s=args.L2_s,L2_e=args.L2_e,outlier_option=args.outlier)

    ########################
    # PARSE MP2 PARAMETERS #
    ########################

    # Deprecated

    print("\n{}".format("*"*144))
    print("* {:^140s} *".format("VDW Parameter Parse Complete!"))
    print("{}".format("*"*144))

    quit()
# Checks for job completion in the discovered orca output files
def keep_complete(Names):
    complete = []
    for i in Names:
        with open(i,'r') as f:
            for lines in f:
                if "SCF NOT CONVERGED AFTER" in lines:
                    break
                if "****ORCA TERMINATED NORMALLY****" in lines:
                    complete += [i]
    return complete

# Scrapes the DFT-D and BSSE-corrected MP2 interaction energies from the discovered orca outputs
# Inputs:      Names      : a list of relative filenames from the working directory
# Returns:     Data       : a dictionary with keys for each filename in "Names", and subkeys holding
#                           each molecule's geometry (geo_a, geo_b), list of types (types_a,types_b)
#                           charge (charges_a, charges_b), and dft-d and BSSE-corrected MP2 interaction
#                           energies (e_dft,e_mp2)
# Dict_info:  "geo_a"     : array holding the geometry of the first molecule
#             "geo_b"     : array holding the geometry of the second molecule
#             "types_a"   : array holding the atom types of the first molecule
#             "types_b"   : array holding the atom types of the second molecule
#             "charges_a" : array holding the atomic partial charges of the first molecule
#             "charges_b" : array holding the atomic partial charges of the second molecule
#             "e_dft"     : array holding the DFT-D interaction energies for each configuration
#             "e_MP2"     : array holding the counterpoise-corrected MP2 interaction energies for each configuration
#             "e_dft_fit" : array holding the DFT-D interaction energies less the electrostatic contribution for each configuration
#             "e_MP2_fit" : array holding the counterpoise-corrected MP2 interaction energies less the electrostatic contribution for each configuration
#             "r_dist"    : array holding the pair-separations for each configurations (rows are indexed to the atoms in geo_a and columns to the atoms in geo_b).
#
def get_data(Names,functional,parse_charges=1,db_files=[],q_a=None,q_b=None,tag=''):

    # Initialize dictionary for all configurations
    global Data
    Data = {}

    # Initialize dictionary to hold charges being read from the database file
    FF_charges = {}
    if db_files != []:
        FF_charges = grab_db_charges(db_files)

    if parse_charges == 1:

        # Initialize a list of unique atom types
        unique_types = []

        # Iterate over all input files and parse atom types and charges for each molecule
        for count_i,i in enumerate(Names):

            # Initialize dictionary for this configuration
            Data[i] = { 'geo_a': np.array([]), 'geo_b': np.array([]), 'types_a': [], 'types_b': [], 'charges_a':np.array([]), 'charges_b':np.array([]), 'e_dft':0.0, 'e_mp2':0.0, 
                        'e_dft_fit_AA':0.0, 'e_mp2_fit_AA':0.0, 'e_dft_fit_UA':0.0, 'e_mp2_fit_UA':0.0, 'charges_a_UA':np.array([]), 'charges_b_UA':np.array([]) }

            # Grab geometry and datatypes from xyz file
            XYZ_file = '.'.join([ j for j in i.split('.')[:-1] ])+'.xyz'
            Data[i]['elem_a'],Data[i]['elem_b'],Data[i]['geo_a'],Data[i]['geo_b'],Data[i]['types_a'],Data[i]['types_b'] = scrape_xyz(XYZ_file,parse_charges=0)

            # Add atom types to unique list
            for j in Data[i]['types_a']:
                if j not in unique_types: unique_types += [j]
            for j in Data[i]['types_b']:
                if j not in unique_types: unique_types += [j]
                    
            # Get charges from output files
            Data[i]['charges_a'],Data[i]['charges_b'] = scrape_charges(i,Data[i]['types_a'],Data[i]['types_b'])

        # Calculate average of each charge type over all configurations
        charges = np.zeros(len(unique_types))
        counts  = np.zeros(len(unique_types))
        for count_i,i in enumerate(unique_types):
            for count_j,j in enumerate(Names):
                for count_k,k in enumerate(Data[j]["types_a"]):
                    if i == k:
                        charges[count_i] += Data[j]["charges_a"][count_k]
                        counts[count_i]  += 1
                for count_k,k in enumerate(Data[j]["types_b"]):
                    if i == k:
                        charges[count_i] += Data[j]["charges_b"][count_k]
                        counts[count_i]  += 1
                        
        # Perform average 
        charges = charges/counts
        
        # Assign charges to each sub-dictionary
        for count_i,i in enumerate(unique_types):
            for count_j,j in enumerate(Names):
                for count_k,k in enumerate(Data[j]["types_a"]):
                    if i == k:
                        Data[j]["charges_a"][count_k] = charges[count_i]
                for count_k,k in enumerate(Data[j]["types_b"]):
                    if i == k:
                        Data[j]["charges_b"][count_k] = charges[count_i]
                
    # Parse configurational energies
    flag_a_AA_charges = 0
    flag_b_AA_charges = 0
    flag_a_UA_charges = 0
    flag_b_UA_charges = 0
    residual_flag_a = 0
    residual_flag_b = 0
    for count_i,i in enumerate(Names):

        # Initialize dictionary for this configuration
        if parse_charges == 0:
            Data[i] = { 'geo_a': np.array([]), 'geo_b': np.array([]), 'types_a': [], 'types_b': [], 'e_dft':0.0, 'e_mp2':0.0, 
                        'e_dft_fit_AA':0.0, 'e_mp2_fit_AA':0.0, 'e_dft_fit_UA':0.0, 'e_mp2_fit_UA':0.0, 'charges_a_UA':np.array([]), 'charges_b_UA':np.array([]) }

        # Grab geometry and datatypes from xyz file
        XYZ_file = '.'.join([ j for j in i.split('.')[:-1] ])+'.xyz'
        XYZdrude_file = '.'.join([ j for j in i.split('.')[:-1] ])+'_drude.xyz'
        if tag != '':
            XYZdrude_file = '.'.join([ j for j in i.split('.')[:-1] ])+'_drude_'+tag+'.xyz'
        if parse_charges == 1:
            Data[i]['elem_a'],Data[i]['elem_b'],Data[i]['geo_a'],Data[i]['geo_b'],Data[i]['types_a'],Data[i]['types_b'] = scrape_xyz(XYZ_file,parse_charges=0)
        else:
            Data[i]['elem_a'],Data[i]['elem_b'],Data[i]['geo_a'],Data[i]['geo_b'],Data[i]['types_a'],Data[i]['types_b'],Data[i]['charges_a'],Data[i]['charges_b'] = scrape_xyz(XYZ_file,parse_charges=1)
            Data[i]['elem_a_drude'],Data[i]['elem_b_drude'],Data[i]['geo_a_drude'],Data[i]['geo_b_drude'],Data[i]['charges_a_drude'],Data[i]['charges_b_drude'] = scrape_xyz_drude(XYZdrude_file,parse_charges=1)


        # Overwrite values for UA-molecules already in the FF_charges dictionary (a check is performed to ensure that all of the atom types
        # are in the FF parameters read from file for each molecule. If they are then the charges are read from file.)
        # NOTE: values are first initialized from a copy of the charges* arrays
        Data[i]["charges_a_UA"] = np.copy(Data[i]["charges_a"])
        Data[i]["charges_b_UA"] = np.copy(Data[i]["charges_b"])
        if False not in [ j+'-UA' in list(FF_charges.keys()) for j in Data[i]['types_a'] if int(j.split('[')[1].split(']')[0]) != 1 ]:
            flag_a_UA_charges = 1
            for count_j,j in enumerate(Data[i]['types_a']):
                if int(j.split('[')[1].split(']')[0]) == 1: continue
                Data[i]['charges_a_UA'][count_j] = FF_charges[j+'-UA']
        if False not in [ j+'-UA' in list(FF_charges.keys()) for j in Data[i]['types_b'] if int(j.split('[')[1].split(']')[0]) != 1 ]:
            flag_b_UA_charges = 1
            for count_j,j in enumerate(Data[i]['types_b']):
                if int(j.split('[')[1].split(']')[0]) == 1: continue
                Data[i]['charges_b_UA'][count_j] = FF_charges[j+'-UA']

        # Assign the total charge on molecule a
        if q_a is None:
            q_a = round(sum(Data[i]['charges_a']))
        else:
            q_a = float(q_a)

        # Assign the total charge on molecule b
        if q_b is None:
            q_b = round(sum(Data[i]['charges_b']))
        else:
            q_b = float(q_b)

        # Balance the charges to ensure an integer number of electrons on each molecule
        # this can potentially happen if the atom types being read from the database file
        # were derived from a different molecule
        residual   = sum(Data[i]['charges_a']) - q_a
        for count_j,j in enumerate(Data[i]['types_a']):

            # Check to see that the amount of charge being redistributed is small            
            if abs(residual/float(len(Data[i]['types_a']))) > 0.05: 
                residual_flag_a = 1

            # Equally redistribute the residual charge
            Data[i]['charges_a'][count_j] -= residual/float(len(Data[i]['types_a']))

        residual   = sum(Data[i]['charges_b']) - q_b
        for count_j,j in enumerate(Data[i]['types_b']):

            # Check to see that the amount of charge being redistributed is small
            if abs(residual/float(len(Data[i]['types_b']))) > 0.05: 
                residual_flag_b = 1

            # Equally redistribute the residual charge
            Data[i]['charges_b'][count_j] -= residual/float(len(Data[i]['types_b']))

        # Calculate UA charges
        Data[i]['adj_mat_a'] = Table_generator(Data[i]['elem_a'],Data[i]['geo_a'])            
        Data[i]['adj_mat_b'] = Table_generator(Data[i]['elem_b'],Data[i]['geo_b'])            

        # Add hydrogen charges to bonded carbons - for molecule a
        Data[i]["charges_a_UA"] = np.copy(Data[i]["charges_a"])
        for count_j,j in enumerate(Data[i]['adj_mat_a']):
            if int(Data[i]["types_a"][count_j].split('[')[1].split(']')[0]) != 6: continue
            H_ind = [ count_k for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_a'][count_k].split('[')[1].split(']')[0]) == 1 ]
            Data[i]["charges_a_UA"][count_j] += sum([ Data[i]["charges_a"][k] for k in H_ind ])
            Data[i]["charges_a_UA"][H_ind] = 0.0

        # Add hydrogen charges to bonded carbons - for molecule b
        Data[i]["charges_b_UA"] = np.copy(Data[i]["charges_b"])
        for count_j,j in enumerate(Data[i]['adj_mat_b']):
            if int(Data[i]["types_b"][count_j].split('[')[1].split(']')[0]) != 6: continue
            H_ind = [ count_k for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_b'][count_k].split('[')[1].split(']')[0]) == 1 ]
            Data[i]["charges_b_UA"][count_j] += sum([ Data[i]["charges_b"][k] for k in H_ind ])
            Data[i]["charges_b_UA"][H_ind] = 0.0
        
        # Determine job type (i.e. 'dft' 'mp2' or combined 'dft' 'mp2')
        if count_i == 0:
            job_type = find_jobtype(i,functional)

        # Parse the orca output files 
        flag = 0
        DFT_energy_count = 0
        MP2_energy_count = 0
        DFT_AB = 0.0
        DFT_A  = 0.0
        DFT_B  = 0.0
        MP2_AB = 0.0
        MP2_A  = 0.0
        MP2_B  = 0.0
        with open(i,'r') as output:
            for lines in output:
                fields = lines.split()

                # The ordering of jobs for 'dft mp2': 1=dimer_DFT, 2=dimer_MP2, 3=A_DFT, 4=A_MP2, 5=B_DFT, 6=B_MP2
                # The ordering of jobs for 'dft':     1=dimer_DFT, 2=A_DFT, 3=B_DFT
                # The ordering of jobs for 'dft mp2': 1=dimer_MP2, 2=A_MP2, 3=B_MP2
                if len(fields) == 5 and fields[0] == "FINAL" and fields[1] == "SINGLE" and fields[2] == "POINT" and fields[3] == "ENERGY":
                    DFT_energy_count += 1
                    # If mp2 jobs are interspersed then the HF jobs will also print "FINAL SINGLE POINT ENERGY"
                    if "dft" in job_type and "mp2" in job_type:
                        if DFT_energy_count == 1:
                            DFT_AB = float(fields[4])
                        elif DFT_energy_count == 3:
                            DFT_A  = float(fields[4])
                        elif DFT_energy_count == 5:
                            DFT_B  = float(fields[4])
                    # If only dft jobs are run then each FINAL SINGLE POINT ENERGY corresponds to one of the DFT jobs
                    elif "dft" in job_type:
                        if DFT_energy_count == 1:
                            DFT_AB = float(fields[4])
                        elif DFT_energy_count == 2:
                            DFT_A  = float(fields[4])
                        elif DFT_energy_count == 3:
                            DFT_B  = float(fields[4])

                if len(fields) == 5 and fields[0] == "MP2" and fields[1] == "TOTAL" and fields[2] == "ENERGY:":
                    MP2_energy_count += 1
                    if MP2_energy_count == 1:
                        MP2_AB = float(fields[3])
                    elif MP2_energy_count == 2:
                        MP2_A  = float(fields[3])
                    elif MP2_energy_count == 3:
                        MP2_B  = float(fields[3])
                        break

        if "dft" in job_type and "mp2" in job_type:
            if DFT_energy_count < 5:
                print("ERROR: Not enough DFT calculations were discovered in file {}".format(i))
                quit()
            elif MP2_energy_count < 3:
                print("ERROR: Not enough MP2 calculations were discovered in file {}".format(i))
                quit()
        elif "dft" in job_type and DFT_energy_count < 3:
            print("ERROR: Not enough DFT calculations were discovered in file {}".format(i))
            quit()
        elif "mp2" in job_type and MP2_energy_count < 3:
            print("ERROR: Not enough MP2 calculations were discovered in file {}".format(i))
            quit()
            
        # Calculate counterpoise-corrected interaction energy
        if "dft" in job_type:
            Data[i]["e_dft"]        = (DFT_AB-DFT_A-DFT_B)*Htokcalmol          # BSSE-corrected DFT interaction energy
        if "mp2" in job_type:
            Data[i]["e_mp2"]        = (MP2_AB-MP2_A-MP2_B)*Htokcalmol          # BSSE-corrected MP2 interaction energy

        # Calculate r_dist matrix, holding the pair-wise separations between atoms on molecule a and b
        Data[i]["r_dist"] = cdist(Data[i]['geo_a'],Data[i]['geo_b'])

        # Calculate pairs lists, holding the pair-wise separations for all instances of each pair-type in this configuration
        Data[i]["pairs"] = {}
        for count_j,j in enumerate(Data[i]["types_a"]):
            for count_k,k in enumerate(Data[i]["types_b"]):
                if (j,k) not in list(Data[i]["pairs"].keys()): Data[i]["pairs"][(j,k)]=[]
                Data[i]["pairs"][(j,k)]+=[Data[i]["r_dist"][count_j,count_k]]

        # Convert lists to arrays
        for j in list(Data[i]["pairs"].keys()):
            Data[i]["pairs"][j] = np.array(Data[i]["pairs"][j])

        # Calculate arrays for vectorized calculation of pairwise interactions
        # "pair_vector" holds the sum of 1/r elements for each pair type indexed to 
        # the "pair_type_vector" pair_types. That is, if element j in pair_type_vector
        # is type (x,y), then element j in "pair_vector" is sum(1/r) instances for 
        # type (x,y) in configuration i. if/else statements are there to remove redundancies.
        Data[i]["pair_vector-6"] = []
        Data[i]["pair_vector-12"] = []
        Data[i]["pair_type_vector"] = []
        added_types = []
        for j in list(Data[i]["pairs"].keys()):
            if j in added_types:
                continue
            if j[0]==j[1]:
                Data[i]["pair_vector-6"]  += [sum(1.0/Data[i]["pairs"][j]**6.0)]
                Data[i]["pair_vector-12"] += [sum(1.0/Data[i]["pairs"][j]**12.0)]
            elif (j[1],j[0]) in list(Data[i]["pairs"].keys()):
                Data[i]["pair_vector-6"]  += [sum(1.0/Data[i]["pairs"][j]**6.0)  + sum(1.0/Data[i]["pairs"][(j[1],j[0])]**6.0)]
                Data[i]["pair_vector-12"] += [sum(1.0/Data[i]["pairs"][j]**12.0) + sum(1.0/Data[i]["pairs"][(j[1],j[0])]**12.0)]                
            else:
                Data[i]["pair_vector-6"]  += [sum(1.0/Data[i]["pairs"][j]**6.0)]
                Data[i]["pair_vector-12"] += [sum(1.0/Data[i]["pairs"][j]**12.0)]
            if j[0] >= j[1]:
                Data[i]["pair_type_vector"]+=[j]
            else:
                Data[i]["pair_type_vector"]+=[(j[1],j[0])]
            added_types += [j,(j[1],j[0])]
        Data[i]["pair_vector-6"]  = np.array(Data[i]["pair_vector-6"])
        Data[i]["pair_vector-12"] = np.array(Data[i]["pair_vector-12"])

    # Print relevant diagnostics
    if residual_flag_a != 0 or residual_flag_b != 0 or flag_a_AA_charges == 1 or flag_b_AA_charges == 1 or flag_a_UA_charges == 1 or flag_b_UA_charges == 1:
        print(" ")
    if residual_flag_a != 0 or residual_flag_b != 0:
        print("{}\n! {:^140} !\n{}\n".format("!"*144,"More than abs(0.05) electrons per atom were required to reach an integer number of electrons","!"*144))
    if flag_a_AA_charges == 1:
        print("AA-charges were read for the a-molecule(s) from file(s): {}".format(', '.join(db_files)))
    if flag_b_AA_charges == 1:
        print("AA-charges were read for the b-molecule(s) from file(s): {}".format(', '.join(db_files)))
    if flag_a_UA_charges == 1:
        print("UA-charges were read for the a-molecule(s) from file(s): {}".format(', '.join(db_files)))
    if flag_b_UA_charges == 1:
        print("UA-charges were read for the b-molecule(s) from file(s): {}".format(', '.join(db_files)))

    ################################################################################
    # Perform charge rescaling to conserve the total electrostatics in the UA_fits #
    ################################################################################

    # Only rescale the UA charges if none of the molecules being fit have their UA_charges read from the FF dictionary.
    # (those being read are expected to be fixed, so charge rescaling needs to be avoided)
    if flag_a_UA_charges == 0 and flag_b_UA_charges == 0:

        # Calculate total intermolecular AA-electrostatic energy
        E_C_Tot_AA = 0.0
        for i in Names:
            E_C_Tot_AA += E_C_Tot(i,'AA')

        # Calculate total intermolecular UA-electrostatic energy (uncorrected)
        E_C_Tot_UA = 0.0
        for i in Names:
            E_C_Tot_UA += E_C_Tot(i,'UA')

        # calculate scale factor (the scale factor is the sqrt of the ratio because coulomb evaluations involve the product of charges)
        scale_factor = (E_C_Tot_AA/E_C_Tot_UA)**(0.5)

        # Print diagnostic and calculate scale factor 
        print("\n{:40s} {:< 12.6f} (kcal/mol)".format("All-atom EC_tot:",E_C_Tot_AA))
        print("{:40s} {:< 12.6f} (kcal/mol)".format("United-atom EC_tot:",E_C_Tot_UA))
        print("{:40s} {:< 12.6f}".format("United-atom EC rescaling factor:",scale_factor))


    # Calculate interaction energies with the coulombic contribution subtracted out
    for i in Names:
        if "dft" in job_type:
            Data[i]["e_dft_fit_AA"] = Data[i]["e_dft"] - E_C_Tot_drude(i,'AA')       # DFT fit potential (subtracted electrostatic component)
            Data[i]["e_dft_fit_AA_ori"] = Data[i]["e_dft"] - E_C_Tot(i,'AA')       # DFT fit potential (subtracted electrostatic component)
            Data[i]["e_dft_fit_UA"] = Data[i]["e_dft"] - E_C_Tot(i,'UA')       # DFT fit potential (subtracted electrostatic component)
        if "mp2" in job_type:
            Data[i]["e_mp2_fit_AA"] = Data[i]["e_mp2"] - E_C_Tot(i,'AA')       # MP2 fit potential (subtracted electrostatic component)
            Data[i]["e_mp2_fit_UA"] = Data[i]["e_mp2"] - E_C_Tot(i,'UA')       # MP2 fit potential (subtracted electrostatic component)

    return Data,job_type

# Description: scrape charges from the supplied db file. Any charges read from the file are used regardless of the
#              charge parsing option.
def grab_db_charges(db_files):
    
    FF_charges = {}
    for i in db_files:
        with open(i,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == 'charge':
                    FF_charges[fields[1]] = float(fields[2])

    return FF_charges

# Return the coordinates and atom types from an xyz file
# atom types are expected as a fifth column and molecule
# index is expected in the sixth column. These are atypical
# an xyz file, but automatically outputed by the vdw_gen.py 
# program.
def scrape_xyz(name,parse_charges=0):

    atom_count = 0
    with open(name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()

            # Grab the number of atoms in the geometry and initialize the parsing lists
            if lc == 0:
                elem = ["X"]*int(fields[0])
                geo = np.zeros([int(fields[0]),3]) 
                types = ["X"]*int(fields[0])
                a_ind = []
                b_ind = []
                
                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges =np.zeros([int(fields[0])])

            # Within the geometry block parse the geometry, charges, mol_ids, and atom types from the enhanced xyz file
            if lc > 1 and len(fields) >= 7:
                geo[atom_count,:] = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                if fields[4] == '0':
                    a_ind += [atom_count]
                else:
                    b_ind += [atom_count]
                elem[atom_count] = fields[0]

                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges[atom_count] = float(fields[5])
                types[atom_count] = fields[6]
                atom_count += 1    

    # If flag is enabled return the charges
    if parse_charges == 1:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               [ i for count_i,i in enumerate(types) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(types) if count_i in b_ind ],\
               charges[a_ind],charges[b_ind]

    # Don't return the charges list if the flag isn't enabled
    if parse_charges == 0:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               [ i for count_i,i in enumerate(types) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(types) if count_i in b_ind ]

# same as scrape_xyz but don't read type
# because some drude_xyz file doesn't contain type info
def scrape_xyz_drude(name,parse_charges=0):

    atom_count = 0
    with open(name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()

            # Grab the number of atoms in the geometry and initialize the parsing lists
            if lc == 0:
                elem = ["X"]*int(fields[0])
                geo = np.zeros([int(fields[0]),3]) 
                a_ind = []
                b_ind = []
                
                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges =np.zeros([int(fields[0])])

            # Within the geometry block parse the geometry, charges, mol_ids, and atom types from the enhanced xyz file
            if lc > 1 and len(fields) >= 6:
                geo[atom_count,:] = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                if fields[4] == '0':
                    a_ind += [atom_count]
                else:
                    b_ind += [atom_count]
                elem[atom_count] = fields[0]

                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges[atom_count] = float(fields[5])
                atom_count += 1    

    # If flag is enabled return the charges
    if parse_charges == 1:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               charges[a_ind],charges[b_ind]

    # Don't return the charges list if the flag isn't enabled
    if parse_charges == 0:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind]

# Description: scrape the charges for both molecules from the first CHELPG calculation
#              of the QC calculation. 
def scrape_charges(filename,atomtypes_a,atomtypes_b):

    # Read in the charges
    charges_flag  = 0
    break_flag    = 0
    charges_count = 0
    charges_a = np.zeros(len(atomtypes_a))
    charges_b = np.zeros(len(atomtypes_b))
    with open(filename,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()

            # Break after all charges have been read
            if break_flag == 1:
                break

            # Find the start of the first CHELPG calculation
            if len(fields) == 2 and fields[0] == "CHELPG" and fields[1] == "Charges":
                charges_flag = 1
                continue

            # Necessary because of a "------" filler line in the Orca CHELPG output
            if charges_flag == 1:
                charges_flag = 2
                continue

            # Start parsing charges in molecule b
            if charges_flag == 3:
                charges_b[charges_count-len(charges_a)] = float(fields[3])
                charges_count += 1

                # Once all charges have been read reset the count and flag.
                if charges_count == len(charges_a)+len(charges_b):
                    break_flag = 1

            # Start parsing charges in molecule a
            if charges_flag == 2:
                charges_a[charges_count] = float(fields[3])
                charges_count += 1
                if charges_count == len(charges_a):
                    charges_flag += 1
    

    # Average over like types of molecule a
    for i in sorted(set(atomtypes_a)):
        ind = []
        charge = 0.0

        # find like types and sum charges
        for count_j,j in enumerate(charges_a):
            if atomtypes_a[count_j] == i:
                ind += [count_j]
                charge += j

        # Assign averaged charges to all like types
        for j in ind:
            charges_a[j] = charge/float(len(ind))        

    # Average over like types of molecule b
    for i in sorted(set(atomtypes_b)):
        ind = []
        charge = 0.0

        # find like types and sum charges
        for count_j,j in enumerate(charges_b):
            if atomtypes_b[count_j] == i:
                ind += [count_j]
                charge += j

        # Assign averaged charges to all like types
        for j in ind:
            charges_b[j] = charge/float(len(ind))        

    return charges_a,charges_b

def find_jobtype(name,functional):

    job_type = []
    with open(name,'r') as f:
        for lines in f:
            fields = lines.split()
            if len(fields) >= 4 and fields[3] == functional and 'dft' not in job_type:
                job_type += ['dft']
            if len(fields) >= 4 and fields[3] == "RI-MP2" and 'mp2' not in job_type:
                job_type += ['mp2']
    return job_type

# This function expects data to be a globally defined dictionary. The function
# outputs the coulomb contribution to the interaction energy of a pair of molecules
def E_C_Tot(ind,eng_type):

    # Initialize the energy
    E_C = 0.0


    # Calculate the atoms_A atom_B separations
    r_dist = cdist(Data[str(ind)]['geo_a'],Data[str(ind)]['geo_b'])

    # Cumulatively add up the electrostatic energy
    for count_i,i in enumerate(r_dist):
        for count_j,j in enumerate(i):
            if eng_type == 'AA':
                E_C += ( Data[str(ind)]["charges_a"][count_i]*Data[str(ind)]["charges_b"][count_j] / j )*Ecoul_Const
            elif eng_type == 'UA':
                E_C += ( Data[str(ind)]["charges_a_UA"][count_i]*Data[str(ind)]["charges_b_UA"][count_j] / j )*Ecoul_Const

    return E_C

# This function expects data to be a globally defined dictionary. The function
# outputs the coulomb contribution to the interaction energy of a pair of molecules
def E_C_Tot_drude(ind,eng_type):

    # Initialize the energy
    E_C = 0.0


    # Calculate the atoms_A atom_B separations
    r_dist = cdist(Data[str(ind)]['geo_a_drude'],Data[str(ind)]['geo_b_drude'])

    # Cumulatively add up the electrostatic energy
    for count_i,i in enumerate(r_dist):
        for count_j,j in enumerate(i):
            if eng_type == 'AA':
                E_C += ( Data[str(ind)]["charges_a_drude"][count_i]*Data[str(ind)]["charges_b_drude"][count_j] / j )*Ecoul_Const


    return E_C

# Description: This function searches through the configurations being parsed for unique
#              molecules. All FF_parameters are linked examples of the molecule(s) used
#              for their parametrization. 
# Note:        identifying unique molecules is generally a graph isomorphism problem
#              (P difficult). The algorithm used here assumes that replica molecules have
#              same atomic ordering (which the vdw_gen.py parser should guarrantee) which 
#              makes identifying the molecule types much easier.
#             
# Inputs:      folder: folder for saving the molecule.db file
#              Data:   the dictionary holding the configuration data 
def parse_molecules(folder,Data):

    # Import min_pairs
    global min_pairs

    # All atom types are specified in terms of atomic number. This dictionary is used to convert them to the corresponding element label. 
    atom2element_dict = { 1:'H' ,  2:'He',
                          3:'Li',  4:'Be',                                                                                            5:'B' ,  6:'C' ,  7:'N' ,  8:'O' ,  9:'F' , 10:'Ne',
                         11:'Na', 12:'Mg',                                                                                           13:'Al', 14:'Si', 15:'P' , 16:'S' , 17:'Cl', 18:'Ar',
                         19:'K' , 20:'Ca', 21:'Sc', 22:'Ti', 23:'V' , 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr',
                         37:'Rb', 38:'Sr', 39:'Y' , 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I' , 54:'Xe',
                         55:'Cs', 56:'Ba', 57:'La', 72:'Hf', 73:'Ta', 74:'W' , 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn'}

    # Initialize lists to hold unique adj_mats, type_lists, element_lists, and geometries
    mol_adj_mats = []
    mol_types    = []
    mol_elements = []
    mol_geos     = []
    mol_charges  = []

    # Iterate over all configurations
    for i in list(Data.keys()):

        # Find elements for mol_a and mol_b
        elem_a = [ atom2element_dict[int(j.split('[')[1].split(']')[0])] for j in Data[i]['types_a'] ]
        elem_b = [ atom2element_dict[int(j.split('[')[1].split(']')[0])] for j in Data[i]['types_b'] ]

        # Generate adjacency matrix for the geo_a and geo_b
        adj_mat_a = Table_generator(elem_a,Data[i]['geo_a'])            
        adj_mat_b = Table_generator(elem_b,Data[i]['geo_b'])            

        # if adj_mat_a isn't in mol_adj_mats then it is a new molecule
        if True not in [ np.array_equal(adj_mat_a,j) for j in  mol_adj_mats ]:
            mol_adj_mats += [adj_mat_a]
            mol_types    += [Data[i]['types_a']]
            mol_elements += [elem_a]
            mol_geos     += [Data[i]['geo_a']]
            mol_charges  += [Data[i]['charges_a']]

        # if adj_mat_a is in mol_adj_mats but there isn't a match for type lists then it is a new molecule
        else:
            match = 0
            for count_j,j in enumerate(mol_adj_mats):
                if np.array_equal(adj_mat_a, j) and Data[i]['types_a'] == mol_types[count_j]: match = 1
            if match == 0:
                mol_adj_mats += [adj_mat_a]
                mol_types    += [Data[i]['types_a']]
                mol_elements += [elem_a]
                mol_geos     += [Data[i]['geo_a']]
                mol_charges  += [Data[i]['charges_a']]

        # if adj_mat_b isn't in mol_adj_mats then it is a new molecule
        if True not in [ np.array_equal(adj_mat_b,j) for j in  mol_adj_mats ]:
            mol_adj_mats += [adj_mat_b]
            mol_types    += [Data[i]['types_b']]
            mol_elements += [elem_b]
            mol_geos     += [Data[i]['geo_b']]
            mol_charges  += [Data[i]['charges_b']]
        
        # if adj_mat_b is in mol_adj_mats but there isn't a match for type lists then it is a new molecule
        else:
            match = 0
            for count_j,j in enumerate(mol_adj_mats):
                if np.array_equal(adj_mat_b, j) and Data[i]['types_b'] == mol_types[count_j]: match = 1
            if match == 0:
                mol_adj_mats += [adj_mat_b]
                mol_types    += [Data[i]['types_b']]
                mol_elements += [elem_b]
                mol_geos     += [Data[i]['geo_b']]
                mol_charges  += [Data[i]['charges_b']]

    # Write examples of each molecule in the sim to a molecule.db file. This is used when
    # adding the VDW parameters to the master database. Each parameter in the master is linked
    # to the molecule(s) used to generate it. In the case of VDW params, one example configuration
    # for each molecule type is used. 
    with open(folder+'/molecules.db','w') as f:
        for count_i,i in enumerate(mol_geos):
            f.write("\nmol {:6d} start\n".format(count_i+1))
            # Save an xyz for viewing the configuration in VMD
            f.write('{}\n\n'.format(len(i)))
            for count_j,j in enumerate(i):
                f.write('  {:<10s} {:< 20.6f} {:< 20.6f} {:< 20.6f} {:< 12.6f} {:<60s}\n'.\
                format(mol_elements[count_i][count_j],j[0],j[1],j[2],mol_charges[count_i][count_j],mol_types[count_i][count_j]))
            f.write("mol {:6d} end\n".format(count_i+1))

    # Assemble the list of atomtypes whose molecules are minimal structures
    min_list = []
    for count_i,i in enumerate(mol_types):
        for count_j,j in enumerate(i):
            if minimal_structure(j,mol_geos[count_i],mol_elements[count_i],adj_mat=mol_adj_mats[count_i],atomtypes=mol_types[count_i],gens=2) is True:
                min_list += [j]

    # Assemble all permutations of minimal types                
    min_list = set(min_list)
    min_pairs = []
    for i in min_list:
        for j in min_list:
            min_pairs += [(i,j),(i+'-UA',j+'-UA')]            
    min_pairs = set(min_pairs)

    # Print diagnostic and return
    print("{:70s} {}".format("Number of unique molecules in the fit configurations:",len(mol_types)))
    return

# Description:   Checks is the supplied geometry corresponds to the minimal structure of the molecule
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
            print("ERROR in minimal_structure: While trying to automatically assign atomtypes, the elements argument must have dimensions equal to geo. Exiting...")
            quit()

        # Generate the adjacency matrix
        # NOTE: the units are converted back angstroms
        adj_mat = Table_generator(elements,geo)

        # Generate the atomtypes
        atomtypes = id_types(elements,adj_mat,gens)
        
    # Check minimal conditions
    count = 0
    for count_i,i in enumerate(atomtypes):

        # If the current atomtype matches the atomtype being searched for then proceed with minimal geo check
        if i == atomtype:
            count += 1

            # Initialize lists for holding indices in the structure within "gens" bonds of the seed atom (count_i)
            keep_list = [count_i]
            new_list  = [count_i]
            
            # Carry ount a "gens" bond deep search
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

# This function drives the fitting algorithm for UA and AA pairwise interaction
# Calls subfunctions "iterative_fit" and "plot_convergence" to perform the fitting
# and generate convergence plots.
def Fit_pairs(Folder,Data,fit,weight,delta_xhi2_thresh,min_cycles,max_cycles,FF_db,QC_type="DFT",charges_opt=1,UA_vol_scale=1.0,UA_method='radius',r_min_scale=1.0,avoid_read_flag=0,mixing_rule="none",L2_s=0.0,L2_e=0.0,outlier_option=True):

    global VDW_dict,fit_type,Pair_min_dict

    # Get unique pairs and number of occurrences
    Fit_Pairs = []                                                   # Initialize list holding unique pair interactions
    Fit_Pairs_hist = {}                                              # Initialize list indexed to pairs to hold the occurence of each pair interaction (there can be multiple per pair of molecules).
    for i in list(Data.keys()):                                            # Iterate over all configurations
        for j in Data[i]["types_a"]:                                 # Iterate over the list of types in molecule A
            for k in Data[i]["types_b"]:                             # Iterate over the list of types in molecule B
                if j >= k:                                           # Higher ranked elements go first in the definition of the type
                    if (j,k) not in Fit_Pairs:                       # If the tuple isn't in the Pairs list then add it and add an element in Pairs_hist
                        Fit_Pairs += [(j,k)]
                        Fit_Pairs_hist[(j,k)] = 1
                    else:                                            # If the pair already exists Pairs then increment the appropriate 
                        Fit_Pairs_hist[(j,k)] += 1                   # element in Fit_Pairs_hist
                else:
                    if (k,j) not in Fit_Pairs:                       # If the tuple isn't in the Pairs list then add it and add an element in Pairs_hist
                        Fit_Pairs += [(k,j)]
                        Fit_Pairs_hist[(k,j)] = 1
                    else:                                            # If the pair already exists Pairs then increment the appropriate 
                        Fit_Pairs_hist[(k,j)] += 1                   # element in Fit_Pairs_hist  
    
    # Generate a list of UA pairs (pairs between atoms that do not contain hydrogen)
    UA_Pairs = []
    for i in Fit_Pairs:
        if return_UA_H(i[0]) != 1 and return_UA_H(i[1]) != 1:
            UA_Pairs += [i]

    # Read in fixed parameters from the fixed_params.db file
    FIXED_Data = {"vdws":{}}
    if os.path.isfile('fixed_params.db'):         
        FIXED_Data = get_FF_data(['fixed_params.db'])

    # If a FF_db file is supplied then the program avoids fitting any of the
    # pairs that are pairs are in the database
    DB_Pairs = []
    FF_Data = {}
    if FF_db != [] or len(list(FIXED_Data['vdws'].keys())) > 0:
        
        # Read in parameters from user-supplied database
        FF_Data = get_FF_data(FF_db,mixing_rule=mixing_rule)

        # Add FIXED_Data parameters while avoiding overwriting those read from the user-supplied database
        for i in [ j for j in list(FIXED_Data["vdws"].keys()) if j not in list(FF_Data["vdws"].keys()) ]:
            FF_Data["vdws"][i] = FIXED_Data["vdws"][i]

        if fit == 'lj':

            # Move AA-pairs that are already in the database to DB_Pairs
            if mixing_rule == "none":
                move_list  = [ count_i for count_i,i in enumerate(Fit_Pairs) if ( (i[0],i[1],'lj') in list(FF_Data["vdws"].keys()) or (i[1],i[0],'lj') in list(FF_Data["vdws"].keys()) ) ]
            elif mixing_rule in ['wh','lb']:
#                move_list  = [ count_i for count_i,i in enumerate(Fit_Pairs) if ( (i[0],i[1],'lj') in FF_Data["vdws"].keys() or (i[1],i[0],'lj') in FF_Data["vdws"].keys() ) and i[0] == i[1] ]
                move_list  = [ count_i for count_i,i in enumerate(Fit_Pairs) if ( (i[0],i[0],'lj') in list(FF_Data["vdws"].keys()) and (i[1],i[1],'lj') in list(FF_Data["vdws"].keys()) ) ]
            else:
                move_list  = []
            DB_Pairs       = [ i for count_i,i in enumerate(Fit_Pairs) if count_i in move_list ]
            Fit_Pairs      = [ i for count_i,i in enumerate(Fit_Pairs) if count_i not in move_list ]

            # Move UA-pairs that are already in the database to DB_Pairs
            if mixing_rule == "none":
                move_list  = [ count_i for count_i,i in enumerate(UA_Pairs) if ( (i[0]+'-UA',i[1]+'-UA','lj') in list(FF_Data["vdws"].keys()) or (i[1]+'-UA',i[0]+'-UA','lj') in list(FF_Data["vdws"].keys()) ) ]
            elif mixing_rule in ['wh','lb']:
#                move_list  = [ count_i for count_i,i in enumerate(UA_Pairs) if ( (i[0]+'-UA',i[1]+'-UA','lj') in FF_Data["vdws"].keys() or (i[1]+'-UA',i[0]+'-UA','lj') in FF_Data["vdws"].keys() ) and i[0] == i[1] ]
                move_list  = [ count_i for count_i,i in enumerate(UA_Pairs) if ( (i[0]+'-UA',i[0]+'-UA','lj') in list(FF_Data["vdws"].keys()) and (i[1]+'-UA',i[1]+'-UA','lj') in list(FF_Data["vdws"].keys()) ) ]
            else:
                move_list  = []
#            DB_Pairs       = [ i for count_i,i in enumerate(Fit_Pairs) if count_i in move_list ]
            UA_Pairs       = [ i for count_i,i in enumerate(UA_Pairs) if count_i not in move_list ]
            
        if fit == 'buck':

            # Move AA-pairs that are already in the database to DB_Pairs
            if mixing_rule == "none":
                move_list  = [ count_i for count_i,i in enumerate(Fit_Pairs) if ( (i[0],i[1],'buck') in list(FF_Data["vdws"].keys()) or (i[1],i[0],'buck') in list(FF_Data["vdws"].keys()) ) ]
            else:
                move_list  = []
            DB_Pairs       = [ i for count_i,i in enumerate(Fit_Pairs) if count_i in move_list ]
            Fit_Pairs      = [ i for count_i,i in enumerate(Fit_Pairs) if count_i not in move_list ]

            # Move UA-pairs that are already in the database to DB_Pairs
            if mixing_rule == "none":                
                move_list  = [ count_i for count_i,i in enumerate(UA_Pairs) if ( (i[0]+'-UA',i[1]+'-UA','buck') in list(FF_Data["vdws"].keys()) or (i[1]+'-UA',i[0]+'-UA','buck') in list(FF_Data["vdws"].keys()) ) ]
            else:
                move_list  = []
            DB_Pairs      += [ (i[0]+'-UA',i[1]+'-UA') for count_i,i in enumerate(UA_Pairs) if count_i in move_list ]
            UA_Pairs       = [ i for count_i,i in enumerate(UA_Pairs) if count_i not in move_list ]

    # Generate a dictionary of minimum separation for pairs (used during fit to ensure that 
    # the minimum in the fit potential does not occur in a region where there is no sampling)
    Pair_min_dict = {}
    Pair_min_dict_UA = {}
    for i in Fit_Pairs: Pair_min_dict[i] = 1000.0; Pair_min_dict[(i[1],i[0])] = 1000.0
    for i in UA_Pairs: Pair_min_dict[i] = 1000.0; Pair_min_dict[(i[1],i[0])] = 1000.0     # Can eventually be removed once separate Pair_min_dict_AA and Pair_min_dict_UA are coded
    for i in list(Data.keys()):
        for count_j,j in enumerate(Data[i]["r_dist"]):
            for count_k,k in enumerate(j):

                # Order the pair to conform to VDW_dict key definitions
                if Data[i]['types_a'][count_j] >= Data[i]['types_b'][count_k]:
                    pair_type = (Data[i]['types_a'][count_j],Data[i]['types_b'][count_k])
                else:
                    pair_type = (Data[i]['types_b'][count_k],Data[i]['types_a'][count_j])

                # Replace minimum value
                if pair_type in list(Pair_min_dict.keys()) and k < Pair_min_dict[pair_type]: Pair_min_dict[pair_type] = k; Pair_min_dict[(pair_type[1],pair_type[0])] = k

    # Scale Pair_min_dict entries by the r_min_scale factor (when set to 1.0 no change is applied)
    for i in list(Pair_min_dict.keys()):
        Pair_min_dict[i] *= r_min_scale
    
    # Print diagnostics
    print("\nUnique AA pairs being fit (occurences,min_sep,min_sigma):\n")
    for count_i,i in enumerate(Fit_Pairs):
        print("\tpair {:80}: {:<10d}  {:<10.4f} {:<10.4f}".format(str(i),Fit_Pairs_hist[i],Pair_min_dict[i],Pair_min_dict[i]*2.0**(-1.0/6.0)))
    print("\nUnique UA pairs being fit (occurences):\n")
    for i in UA_Pairs:
        print("\tpair {:80}: {:<10d}".format(str(i),Fit_Pairs_hist[i]))

    # Generate a separation vs number histogam for all pair types
    print("\n\tGenerating separation histograms...")
#    pair_hist_plot(Folder+'/figures',Fit_Pairs,Data,QC_type=QC_type)

    ###############################
    # Fit DFT all-atom parameters #
    ###############################
    if QC_type == "DFT":
        
        print("\n{}".format("*"*144))
        print("* {:^140s} *".format("Fit DFT All-Atom VDW Parameters (Self-Consistent Pair-Type Fits to All Configurations)"))
        print("{}".format("*"*144))

        # Initialize fit arrays
        x_vals = list(Data.keys())
        y_vals = [ Data[i]["e_dft_fit_AA"] for i in x_vals ]

        # Initialize VDW_dict based on DFT_FF and UFF parameters    
        VDW_dict = init_VDW(Data,Fit_Pairs,DB_Pairs,Pair_min_dict,fit,FF_db,{},UA=0,avoid_read_flag=avoid_read_flag,mixing_rule=mixing_rule)

        # Initialize UFF reference dictionary
        VDW_0 = init_VDW(Data,Fit_Pairs,[],Pair_min_dict,fit,[],{},UA=0,avoid_read_flag=1,mixing_rule=mixing_rule,verbose=False)        

        # Iteratively fit all parameters        
        fit_vals = global_fit(x_vals,y_vals,Fit_Pairs,weight,fit,delta_xhi2_thresh,min_cycles,max_cycles,mixing_rule=mixing_rule,L2_s=L2_s,L2_e=L2_e,VDW_0=VDW_0,outlier_option=outlier_option)

        # Calculate the VDW energy with the all-atom parameters
        E_VDW_before = 0.0
        if fit == 'lj':
            for i in x_vals:
                E_VDW_before += E_LJ_Tot(i)
        if fit == 'buck':
            for i in x_vals:
                E_VDW_before += E_BUCK_Tot(i)

        # Print diagnostic indicating which pairs probably have sampling errors based if their sigma values are within
        # 0.05% of their minimum separation threshold (i.e. the element in Pair_min_dict generated based on the configurational distribution)
        Warning_list = [ i for i in Fit_Pairs if VDW_dict[i][1]*2.0**(1.0/6.0) < Pair_min_dict[i]*1.05 ]
        if len(Warning_list) > 0:
            print("{}".format("*"*144))
            print("* {:^140s} *".format("The following pairs have minima in their LJ interaction within 5% of their closest sampled separation."))
            print("* {:^140s} *".format("This is usually indicative of insufficient sampling."))
            print("{}".format("*"*144))
            print("\n\t{:<80s} {:<20s} {:<20s} {:<20s}".format("Pair","sigma (ang)","r_min (ang)","min_sep (ang)"))
            for i in Warning_list:
                print("\t{:80s} {:<20.6f} {:<20.6f} {:<20.6f}".format(str(i),VDW_dict[i][1],VDW_dict[i][1]*2.0**(1.0/6.0),Pair_min_dict[i]))
            print("")

        # Write VDW parameters to database file (note VDW_dict is a global variable so the subfunction has access to the parameters)
        write_params(Folder+'/DFT-AA',fit,type='AA',write_charges=charges_opt) 
#        write_params(Folder+'/DFT-AA.db',fit,type='AA',write_charges=1) # Always writes parameters

        # Generate individual plots of the fit potential
        plot_fit_potentials(Folder+'/figures',fit,type='DFT-AA',r_min=0.1,r_max=10.0)

        # Generate convergence plot
        plot_convergence(Folder+'/figures',Data,fit_vals,plot_num=5,type='DFT-AA')

    ##################################
    # Fit DFT united-atom parameters #
    ##################################

    # depricateed

    ###############################
    # Fit MP2 all-atom parameters #
    ###############################

    # depricated

    ##################################
    # Fit MP2 united-atom parameters #
    ##################################

    # depricated

    return

# This function reads in an atomtype label, generates its adjacency matrix, 
# and determines if it is a UA-hydrogen based on the first atomtype and its connections
def return_UA_H(atype):

    # Initialize lists/array 
    brackets = [ i for i in atype if i == "[" or i == "]" ]
    atoms    = atype.replace('['," ").replace(']'," ").split()
    adj_mat = np.zeros([len([ i for i in brackets if i == "[" ]),len([ i for i in brackets if i == "[" ])])

    # The basis of adj_mat corresponds to the atoms in "atoms", equivalently the "[" objects in brackets
    for count_i,i in enumerate(adj_mat):        

        # Indices are counted based on the number of left_brackets encountered (i.e. left_count's values)
        # Connections occur when bracket_count is equal to 1 and the current object in brackets is a left-bracket.
        # The algorithm works by incrementing bracket_count by +/- 1 for every bracket encountered.
        bracket_count = -1
        left_count = -1

        # Loop over the brackets list, bracket_count only gets incremented for objects in brackets satisfying left_count >= count_i 
        for count_j,j in enumerate(brackets):
            if j == "[": left_count += 1
            if left_count < count_i: continue
            if j == "[": bracket_count += 1; 
            if j == "]": bracket_count -= 1
            if left_count == count_i: continue
            if bracket_count == 1 and j == "[":
                adj_mat[count_i,left_count] = 1
                adj_mat[left_count,count_i] = 1

    # If the current atom is a hydrogen (i.e. atomtype "1") and it is connected to at least on carbon (i.e. atomtype "6") then it is a UA-type
    if atoms[0] == "1" and len([ i for count_i,i in enumerate(adj_mat[0]) if i == 1 and atoms[count_i] == "6" ]) > 0:
        return 1
    # Else, it is not a UA-type
    else:
        return 0

# Description: Reads in the FF parameters
# NOTE: db_files
def get_FF_data(db_files,keep_types=[ "atom", "vdw", "bond", "angle", "torsion", "dihedral", "charge" ],mixing_rule=None):

    Data = {"atoms":{},"bonds":{},"angles":{},"dihedrals":{},"vdws":{},"charges":{}}
    
    # Read in Atom data
    for i in db_files:
        with open(i,'r') as f:
            for lines in f:
                fields = lines.split()

                if len(fields) == 0: continue
                if fields[0] == "atom" and "atom" in keep_types:
                    Data["atoms"][fields[1]] = fields[1:]
                if fields[0] == "vdw" and "vdw" in keep_types:
                    Data["vdws"][(fields[1],fields[2],fields[3])] = fields[1:]
                if fields[0] == "bond" and "bond" in keep_types:
                    Data["bonds"][(fields[1],fields[2],fields[3])] = fields[1:]
                if fields[0] == "angle" and "angle" in keep_types:
                    Data["angles"][(fields[1],fields[2],fields[3],fields[4])] = fields[1:]
                if fields[0] == "diehdral" or fields[0] == "torsion" and ( "diehdral" in keep_types or "torsion" in keep_types ):
                    Data["dihedrals"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = fields[1:]
                if fields[0] == "charge" and "charge" in keep_types:
                    Data["charges"][fields[1]] = fields[1:]
  
    # Generate the cross terms via the mixing rule if specified
    if mixing_rule == "lb":

        # AA-Types
        self_types = [ i[0] for i in Data["vdws"] if i[0] == i[1] and i[2] == 'lj' and "-UA" not in i[0] ]
        for count_i,i in enumerate(self_types):
            for count_j,j in enumerate(self_types):
                if count_j > count_i:
                    eps   = ( float(Data["vdws"][(i,i,'lj')][3]) * float(Data["vdws"][(j,j,'lj')][3]) )**(0.5)
                    sigma = ( float(Data["vdws"][(i,i,'lj')][4]) + float(Data["vdws"][(j,j,'lj')][4]) ) / 2.0
                    Data["vdws"][(i,j,'lj')] = [i,j,'lj', eps, sigma]
                    Data["vdws"][(j,i,'lj')] = [j,i,'lj', eps, sigma]

        # UA-Types
        self_types = [ i[0] for i in Data["vdws"] if i[0] == i[1] and i[2] == 'lj' and "-UA" in i[0] ]
        for count_i,i in enumerate(self_types):
            for count_j,j in enumerate(self_types):
                if count_j > count_i:
                    eps   = ( float(Data["vdws"][(i,i,'lj')][3]) * float(Data["vdws"][(j,j,'lj')][3]) )**(0.5)
                    sigma = ( float(Data["vdws"][(i,i,'lj')][4]) + float(Data["vdws"][(j,j,'lj')][4]) ) / 2.0
                    Data["vdws"][(i,j,'lj')] = [i,j,'lj', eps, sigma]
                    Data["vdws"][(j,i,'lj')] = [j,i,'lj', eps, sigma]

    if mixing_rule == "wh":

        # AA-Types
        self_types = [ i[0] for i in Data["vdws"] if i[0] == i[1] and i[2] == 'lj' and "-UA" not in i[0] ]
        for count_i,i in enumerate(self_types):
            for count_j,j in enumerate(self_types):
                if count_j > count_i:
                    sigma  = (( float(Data["vdws"][(i,i,'lj')][4])**(6.0) + float(Data["vdws"][(j,j,'lj')][4])**(6.0))/2.0)**(1.0/6.0)
                    eps    = (float(Data["vdws"][(i,i,'lj')][3])*float(Data["vdws"][(i,i,'lj')][4])**(6.0) * float(Data["vdws"][(j,j,'lj')][3])*float(Data["vdws"][(j,j,'lj')][4])**(6.0) )**(0.5) / sigma**(6.0)
                    Data["vdws"][(i,j,'lj')] = [i,j,'lj', eps, sigma]
                    Data["vdws"][(j,i,'lj')] = [j,i,'lj', eps, sigma]

        # UA-Types
        self_types = [ i[0] for i in Data["vdws"] if i[0] == i[1] and i[2] == 'lj' and "-UA" in i[0] ]
        for count_i,i in enumerate(self_types):
            for count_j,j in enumerate(self_types):
                if count_j > count_i:
                    sigma  = (( float(Data["vdws"][(i,i,'lj')][4])**(6.0) + float(Data["vdws"][(j,j,'lj')][4])**(6.0))/2.0)**(1.0/6.0)
                    eps    = (float(Data["vdws"][(i,i,'lj')][3])*float(Data["vdws"][(i,i,'lj')][4])**(6.0) * float(Data["vdws"][(j,j,'lj')][3])*float(Data["vdws"][(j,j,'lj')][4])**(6.0) )**(0.5) / sigma**(6.0)
                    Data["vdws"][(i,j,'lj')] = [i,j,'lj', eps, sigma]
                    Data["vdws"][(j,i,'lj')] = [j,i,'lj', eps, sigma]
      
    return Data

# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def init_VDW(Data,fit_pairs,read_pairs,Pair_min_dict,fit,FF_db,VDW_dict_prev,UA=0,avoid_read_flag=0,mixing_rule="lb",verbose=True):

    if verbose is True:
        if avoid_read_flag == 0 and mixing_rule in ["none","lb"]:
            print("\nReading VDW parameters or initializing based on UFF values and Lorentz-Berthelot mixing rules\n")
        elif avoid_read_flag == 0 and mixing_rule in ["wh"]:
            print("\nReading VDW parameters or initializing based on UFF values and Waldman-Hagler mixing rules\n")
        elif avoid_read_flag == 1 and mixing_rule in ["none","lb"]:
            print("\nInitializing VDW parameters based on UFF values and Lorentz-Berthelot mixing rules\n")
        elif avoid_read_flag == 1 and mixing_rule in ["wh"]:
            print("\nInitializing VDW parameters based on UFF values and Waldman-Hagler mixing rules\n")

        print("\t{:85s}  {:22s} {:20s}".format("pair_interaction","params","origin"))

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

    # Initialize initial_guess dictionary based on the initial_params.db file (if present). 
    # The default is to use these parameters in place of UFF if the file is present.     
    INIT_dict = {}
    if os.path.isfile('initial_params.db') and avoid_read_flag == 0:
        with open('initial_params.db','r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "vdw":
                    INIT_dict[(fields[1],fields[2])] = [ float(i) for i in fields[4:] ]
                    INIT_dict[(fields[2],fields[1])] = [ float(i) for i in fields[4:] ]

    # Collect types that will be converted into UA
    if UA==1:
        UA_types = {}
        UA_H_types = []
        for i in list(Data.keys()):

            # Iterate over the molecule a adj_mat and identify hydrogens attached to carbons. 
            for count_j,j in enumerate(Data[i]["adj_mat_a"]):
                if int(Data[i]["types_a"][count_j].split('[')[1].split(']')[0]) != 6: continue
                UA_H_types += [ Data[i]["types_a"][count_k] for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_a'][count_k].split('[')[1].split(']')[0]) == 1 ]

            # Iterate over the molecule b adj_mat and identify hydrogens attached to carbons.
            for count_j,j in enumerate(Data[i]["adj_mat_b"]):
                if int(Data[i]["types_b"][count_j].split('[')[1].split(']')[0]) != 6: continue
                UA_H_types += [ Data[i]["types_b"][count_k] for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_b'][count_k].split('[')[1].split(']')[0]) == 1 ]

    # Initialize VDW_dict with fit types, first guess based on element types and Lorentz-Berthelot mixing rules
    VDW_dict = {}
    for i in fit_pairs:

        # Grab the base atomic number from the atomtype
        type_1 = int(i[0].split('[')[1].split(']')[0])
        type_2 = int(i[1].split('[')[1].split(']')[0])
        
        # for UA fits, the eps of Hydrogen containing pairs are set to zero
        # so that they do not participate in the fitting. 
        if UA==1 and ( i[0] in UA_H_types or i[1] in UA_H_types ):
            eps    = 0.0
            sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0
            VDW_dict[i] = (eps,sigma)
            VDW_dict[(i[1],i[0])] = (eps,sigma)

            # Print diagnostic
            if verbose is True:
                print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"Zeroed"))

        # Non-hydrogen UA types read their values from the AA-fit
        elif UA==1:
            eps    = VDW_dict_prev[i][0]
            sigma  = VDW_dict_prev[i][1]
            VDW_dict[i] = (eps,sigma)
            VDW_dict[(i[1],i[0])] = (eps,sigma)

            # Print diagnostic
            if verbose is True:
                print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"AA-Fit"))

        # For all atom fits, initial guesses for eps and sigma are based on scaled UFF parameters with a 
        # check to make sure that the radii aren't too small for the fit set. 
        else:
            if i in INIT_dict:
                eps    = INIT_dict[i][0]
                sigma  = INIT_dict[i][1]
                if sigma < Pair_min_dict[i]*2.**(-1./6.): sigma = Pair_min_dict[i]*2.**(-1./6.)
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)
                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"initial_params.db"))
            else:
                if mixing_rule in ["lb",'none']:
                    eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) 
                    sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0
                elif mixing_rule == "wh":
                    sigma  = ((UFF_dict[type_1][1]**(6.0)+UFF_dict[type_2][1]**(6.0))/2.0)**(1.0/6.0)
                    eps    = (UFF_dict[type_1][0]*UFF_dict[type_1][1]**(6.0) * UFF_dict[type_2][0]*UFF_dict[type_2][1]**(6.0) )**(0.5) / sigma**(6.0)
                else: 
                    print("ERROR in init_VDW: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
                    quit()
                if sigma < Pair_min_dict[i]*2.**(-1./6.): sigma = Pair_min_dict[i]*2.**(-1./6.)
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)
        
                # Print diagnostic
                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"UFF"))

    # For diagnostic purposes, the fixed_params.db and user_supplied FF files are reparsed here in order to 
    # print the proper origin information for the user.
    # Read in fixed parameters from the fixed_params.db file
    FIXED_dict = {"vdws":{}}
    if os.path.isfile('fixed_params.db'):         
        FIXED_dict = get_FF_data(['fixed_params.db'])

    # If a FF_db file is supplied then the program avoids fitting any of the
    # pairs that are pairs are in the database
    READ_dict = {"vdws":{}}
    if FF_db != []:        
        READ_dict = get_FF_data(FF_db,mixing_rule=mixing_rule)

    # Initialize read_pairs with entries from the FF_db.
    for i in read_pairs:

        # Grab the base atomic number from the atomtype
        type_1 = int(i[0].split('[')[1].split(']')[0])
        type_2 = int(i[1].split('[')[1].split(']')[0])
        
        # for UA fits, the eps of Hydrogen containing pairs are set to zero
        # so that they do not participate in the fitting. 
        if UA==1 and ( i[0] in UA_H_types or i[1] in UA_H_types ):
            VDW_dict[i] = (0.0,1.0)
            VDW_dict[(i[1],i[0])] = (0.0,1.0)
            continue
        else:        

            # Mixing rule protocols
            if mixing_rule != "none": 

                # Use mixing rule protocol for read parameters
                if (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):

                    eps_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for fixed parameters
                elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):

                    eps_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for fixed/read parameters
                elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):

                    eps_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for read/fixed parameters
                elif (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):

                    eps_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][4])

                # Apply mixing rule
                if mixing_rule in ["lb"]:
                    eps    = (eps_1*eps_2)**(0.5)
                    sigma  = (sigma_1+sigma_2)/2.0
                elif mixing_rule == "wh":
                    sigma  =  ( ( sigma_1**(6.0) + sigma_2**(6.0) ) / 2.0 )**(1.0/6.0)
                    eps    =  ( eps_1*sigma_1**(6.0) * eps_2*sigma_2**(6.0) )**(0.5) / sigma**(6.0)    
                else: 
                    print("ERROR in init_VDW: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
                    quit()
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)

                # Print commands
                if verbose is True:
                    if (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"supplied_file"))
                    elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed_params.db"))
                    elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed/supplied mixing"))
                    elif (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed/supplied mixing"))



            # Add the pair type (along with the reverse combination) using the parameters found in the READ_dict
            elif (i[0],i[1],'lj') in list(READ_dict["vdws"].keys()) or (i[1],i[0],'lj') in list(READ_dict["vdws"].keys()):
                if (i[0],i[1],'lj') in list(READ_dict["vdws"].keys()):
                    VDW_dict[i] = (float(READ_dict["vdws"][(i[0],i[1],'lj')][3]),float(READ_dict["vdws"][(i[0],i[1],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(READ_dict["vdws"][(i[0],i[1],'lj')][3]),float(READ_dict["vdws"][(i[0],i[1],'lj')][4]))            

                else:
                    VDW_dict[i] = (float(READ_dict["vdws"][(i[1],i[0],'lj')][3]),float(READ_dict["vdws"][(i[1],i[0],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(READ_dict["vdws"][(i[1],i[0],'lj')][3]),float(READ_dict["vdws"][(i[1],i[0],'lj')][4]))

                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"supplied_file"))

            # Add the pair type (along with the reverse combination) using the parameters found in the FIXED_dict
            elif (i[0],i[1],'lj') in list(FIXED_dict["vdws"].keys()) or (i[1],i[0],'lj') in list(FIXED_dict["vdws"].keys()):

                if (i[0],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                    VDW_dict[i] = (float(FIXED_dict["vdws"][(i[0],i[1],'lj')][3]),float(FIXED_dict["vdws"][(i[0],i[1],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(FIXED_dict["vdws"][(i[0],i[1],'lj')][3]),float(FIXED_dict["vdws"][(i[0],i[1],'lj')][4]))            

                else:
                    VDW_dict[i] = (float(FIXED_dict["vdws"][(i[1],i[0],'lj')][3]),float(FIXED_dict["vdws"][(i[1],i[0],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(FIXED_dict["vdws"][(i[1],i[0],'lj')][3]),float(FIXED_dict["vdws"][(i[1],i[0],'lj')][4]))

                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(str(i),", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed_params.db"))

            else:
                print("ERROR in init_VDW: Information for pair type {} was expected in the the fixed_params.db".format((i[0],i[1],'lj')))
                print("                   and/or the supplied FF files, but it could not be found. Please check that the parameters exist. Exiting...")
                quit()


    # Perform conversions to buckingham style. Conversions follow the approach of Rappe et al. 
    # with a "form factor" of 12.0 which reproduces the long-range LJ behavior
    if fit == 'buck':
        for i in list(VDW_dict.keys()):
            A = VDW_dict[i][0]*0.5*exp(12.0)
            B = VDW_dict[i][1]/12.0
            C = VDW_dict[i][0] * 2.0 * VDW_dict[i][1]**(6.0)
            VDW_dict[i] = (A,B,C)

    return VDW_dict

# Description: This function drives the pair-wise fit of the VDW parameters. It accepts
#              various arguments that modify the convergence criteria and optimization
#              algorithm. 
# 
# Inputs:      x_vals            : a list holding the keys in Data for configurations being fit
#              y_vals            : a list holding the interaction energies for each configuration being fit 
#              Pairs             : a list of the unique pair-types being fit.
#              weight            : a scalar that controls the weight of the update at each step.
#              delta_xhi2_thresh : a scalar that controls the xhi^2 convergence threshold for the fit. 
#              min_cycles        : an integer that controls the minimum number of fit cycles to perform.
#              max_cycles        : an integer that controls the maximum number of fit cycles to perform. 
# 
# Returns:     fit_vals          : a list of arrays, each holding the FF based interaction energies in each cycle of the fit.
#
def global_fit(x_vals,y_vals,Pairs,weight,fit,delta_xhi2_thresh=0.00001,min_cycles=100,max_cycles=1000,mixing_rule='none',L2_s=0.0,L2_e=0.0,VDW_0=None,outlier_option=True):

    # In python you need to expand the scope of global variables at every level
    global VDW_dict,fit_type,Pair_min_dict,Data

    # BOTH OF THE SELF TERMS NEED TO BE DEFINED FOR MIXING RULES TO BE SENSIBLY APPLIED
    if mixing_rule in ['wh','lb']:
        # Pairs = [ (k,k) for k in set([ j for i in Pairs for j in i ]) ] # OLD
        Pairs = [ i for i in Pairs if i[0] == i[1] ] # NEW
    elif mixing_rule != "none":
        print("ERROR in global_fit: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
        quit()

    # Intialize cycle count and array of fit values
    fit_type = 0 # dummy value so that all pair-potentials are taken from VDW_dict (see E_LJ_Fit for more details)
    cycle_count = 0

    if fit == 'lj':
        fit_vals = [ E_LJ_Fit(x_vals,0.0,0.0) ] * (max_cycles+5)
    elif fit == 'buck':
        fit_vals = [ E_BUCK_Fit(x_vals,0.0,0.0,0.0) ] * (max_cycles+5)
    else:
        print("ERROR in global_fit: {} is not an implemented option for the fit variable. Exiting...".format(fit))
        quit()
    xhi2_previous = np.mean((y_vals-fit_vals[cycle_count])**2) 

    # If no pairs are being fit (as happens if only UA types are being fit and AA are being read from file) just return the fit_vals
    if len(Pairs) == 0: print(" "); return fit_vals[:cycle_count+1]

    # Check the rank of the fit data by assembling the lstsq fit matrix
    if fit in ['lj','buck'] :
        A = np.zeros([len(x_vals),len(Pairs)*2]) # For the lj fit there are 2 linearly parameters per pair type (A,B). The rank of the LJ matrix is used in place of the buck 
                                              # A matrix owing to the non-linearity of the buckingham potential
    # Build up A matrix
    for count_i,i in enumerate(x_vals):
        for count_j,j in enumerate(Pairs):

            if fit in ['lj','buck']:
                # Pairs only contains I J or J I combinations whereas Data[i]["pairs"] has both combinations
                # To be consistent the contributions from both I J and J I configurations need to be added to 
                # the elements of A
                tmp_6 = 0.0
                tmp_12 = 0.0

                if (j[0],j[1]) in list(Data[i]["pairs"].keys()):
                    tmp_6  += sum(Data[i]["pairs"][(j[0],j[1])]**(-6.0))
                    tmp_12 += sum(Data[i]["pairs"][(j[0],j[1])]**(-12.0))
                if j[0] != j[1] and (j[1],j[0]) in list(Data[i]["pairs"].keys()):
                    tmp_6  += sum(Data[i]["pairs"][(j[1],j[0])]**(-6.0))
                    tmp_12 += sum(Data[i]["pairs"][(j[1],j[0])]**(-12.0))

                A[count_i,(count_j*2)+0] += tmp_12
                A[count_i,(count_j*2)+1] -= tmp_6

    # Calculate the rank of the fit matrix
    rank = matrix_rank(A)

    # Print diagnostics
    print("\nInitializing L-BFGS-B fit of pair-wise interactions:\n")
    print("\t{:<30s} {}".format("delta_xhi2_thresh:",delta_xhi2_thresh))
    print("\t{:<30s} {}".format("mixing rule:",mixing_rule))
    print("\t{:<30s} {}".format("Fit rank:",rank))
    print("\t{:<30s} {}".format("Fit parameters:",len(A[0])))
    print("\t{:<30s} {}".format("L2 sigma weight:",L2_s))
    print("\t{:<30s} {}".format("L2 epsilon weight:",L2_e))
    print("\t{:<30s} {:<12.6g}".format("Fit Tolerance (kcal/mol):",delta_xhi2_thresh))

    # If Reference VDW dictionary isn't supplied, then the VDW_dict is copied and used for regularization
    if VDW_0 is None:
        VDW_0 = deepcopy(VDW_dict)

    if outlier_option:
        Mi,outlier = detect_outlier(y_vals)   
        print("Number of outliers removed: {}".format(len([out for out in outlier if out])))
        fit_func = lambda x: global_fit_LJ(*x,ind=x_vals,outlier=outlier,E_config=y_vals,VDW_dict=VDW_dict,Data=Data,Fit_Pairs=Pairs,\
                                      Pair_min_dict=Pair_min_dict,penalty=0.0,mixing_rule=mixing_rule,L2_sigma=L2_s,L2_eps=L2_e,VDW_0=VDW_0)

    else:
        # Initialize anonymous fit function
        fit_func = lambda x: global_fit_LJ(*x,ind=x_vals,E_config=y_vals,VDW_dict=VDW_dict,Data=Data,Fit_Pairs=Pairs,\
                                         Pair_min_dict=Pair_min_dict,penalty=0.0,mixing_rule=mixing_rule,L2_sigma=L2_s,L2_eps=L2_e,VDW_0=VDW_0)

    # Initialize guesses and parameter bounds
    initial_guess = []
    initial_guess_sqrt = []
    bounds = []
    bounds_sqrt = []
    for i in Pairs:
        initial_guess += [ j for j in VDW_dict[i] ] # UFF GUESS
        initial_guess_sqrt += [ j**(1/2) for j in VDW_dict[i] ]
#        initial_guess += [ 0.002,Pair_min_dict[i] ] # MIN GUESS
        bounds += [ (0.001,30.0),(Pair_min_dict[i]*2.**(-1./6.),10.0) ]

    # Perform the fit
    t0 = time.time()
    if fit.lower() == "lj":
        results = minimize(fit_func,initial_guess,bounds=bounds,method='L-BFGS-B',options={'maxiter':1000000,'maxfun':1000000,'ftol':delta_xhi2_thresh})
        new_params = results.x

    # Print diagnostics
    if results.success == True:
        print("\t{:<30s} {:<12s}".format("Fit termination condition:","Converged"))
    else:
        print("\t{:<30s} {:<12s}".format("Fit termination condition:","Unconverged"))
    print("\t{:<30s} {:<12.6f}".format("Initial xhi^2 (kcal/mol):",global_fit_LJ(*initial_guess,ind=x_vals,E_config=y_vals,VDW_dict=VDW_dict,Data=Data,Fit_Pairs=Pairs,\
                                                                                 Pair_min_dict=Pair_min_dict,penalty=0.0,mixing_rule=mixing_rule,L2_sigma=0.0,L2_eps=0.0)))
    print("\t{:<30s} {:<12.6f}".format("Final xhi^2 (kcal/mol):",global_fit_LJ(*new_params,ind=x_vals,E_config=y_vals,VDW_dict=VDW_dict,Data=Data,Fit_Pairs=Pairs,\
                                                                               Pair_min_dict=Pair_min_dict,penalty=0.0,mixing_rule=mixing_rule,L2_sigma=0.0,L2_eps=0.0)))
    print("\t{:<30s} {:<12.6f}".format("Time for completion (s):",time.time()-t0))

    # Update the VDW_dict    
    print("\nParameter update summary:\n")
    counter = 0
    for i in Pairs:      
        print("\tpair {:80}: {} ->     {}".format(str(i),", ".join([ "{:<15.4f}".format(j) for j in VDW_dict[i]]),", ".join([ "{:<15.4f}".format(j) for j in new_params[counter:counter+len(VDW_dict[i])] ])))
        VDW_dict[i] = tuple([ j for j in new_params[counter:counter+len(VDW_dict[i])] ])
        VDW_dict[(i[1],i[0])] = VDW_dict[i]
        counter += len(VDW_dict[i])

    # Update cross terms if mixing rules were used
    if mixing_rule == "wh":        
        set_of_types = set([ j for i in Pairs for j in i ])
        updated_pairs = [ i for i in list(VDW_dict.keys()) if i[0] > i[1] and ( i[0] in set_of_types or i[1] in set_of_types ) ]
        for i in updated_pairs:
            old = deepcopy(VDW_dict[i])
            sigma = (( VDW_dict[(i[0],i[0])][1]**(6.0) + VDW_dict[(i[1],i[1])][1]**(6.0) ) / 2.0 )**(1.0/6.0)
            VDW_dict[i] = ( ( VDW_dict[(i[0],i[0])][0]*VDW_dict[(i[0],i[0])][1]**(6.0) * VDW_dict[(i[1],i[1])][0]*VDW_dict[(i[1],i[1])][1]**(6.0) )**(0.5) / sigma**(6.0)  , sigma )
            VDW_dict[(i[1],i[0])] = VDW_dict[i]
            print("\tpair {:80}: {} ->     {}".format(str(i),", ".join([ "{:<15.4f}".format(k) for k in old ]),", ".join([ "{:<15.4f}".format(j) for j in VDW_dict[i] ])))
    if mixing_rule == "lb":
        set_of_types = set([ j for i in Pairs for j in i ])
        updated_pairs = [ i for i in list(VDW_dict.keys()) if i[0] > i[1] and ( i[0] in set_of_types or i[1] in set_of_types ) ]
        for i in updated_pairs:
            old = deepcopy(VDW_dict[i])
            VDW_dict[i] = ( ( VDW_dict[(i[0],i[0])][0]*VDW_dict[(i[1],i[1])][0] )**(0.5) , (VDW_dict[(i[0],i[0])][1]+VDW_dict[(i[1],i[1])][1]) / 2.0)
            VDW_dict[(i[1],i[0])] = VDW_dict[i]
            print("\tpair {:80}: {} ->     {}".format(str(i),", ".join([ "{:<15.4f}".format(k) for k in old ]),", ".join([ "{:<15.4f}".format(j) for j in VDW_dict[i] ])))
                    
    print(" ")

    # Calculate fit values and xhi2 with final fit parameters
    if fit.lower() == "lj":
        fit_vals[1] = E_LJ_Fit(x_vals,0.0,0.0)
    return fit_vals[:2]  # only return the result of using final fit params and UFF 

# Returns the LJ energy for each configuration (as an array) using the parameters currently defined in VDW_dict
# and the supplied parameters for the type being fit (fit_type is globally defined outside of the function). 
def E_LJ_Fit(ind,eps_fit,sigma_fit):

    # In python you need to expand the scope of global variables at every level
    global VDW_dict,fit_type,Pair_min_dict

    # If the attempted value for the fit type violates the sanity checks on the fit range then the energy is set to an astronomical value
    # NOTE: the restriction on sigma_min is to stop the minimum of the fit potential from lieing outside of the sampled configuration space
    if fit_type in list(Pair_min_dict.keys()) and ( eps_fit < 0.001 or eps_fit > 30.0 or sigma_fit > 10.0 or sigma_fit < Pair_min_dict[fit_type]*2.**(-1./6.) ): return np.ones(len(ind))*1.E20        

    # Initialize E_LJ for the current configuration 
    E_LJ = np.zeros(len(ind))

    # Initialize eps and sigma arrays for the vectorized calculations
    eps_array = np.array([ VDW_dict[i][0] for i in Data[ind[0]]["pair_type_vector"] ])
    sigma_array = np.array([ VDW_dict[i][1] for i in Data[ind[0]]["pair_type_vector"] ])

    for i in [ count_j for count_j,j in enumerate(Data[ind[0]]["pair_type_vector"]) if ( j == fit_type or (j[1],j[0]) == fit_type ) ]:        
        eps_array[i]   = eps_fit
        sigma_array[i] = sigma_fit

    # Iterate over the configuration index
    for count_d,d in enumerate(ind):

        # Iterate over list of pair-seps for each pair-type (Data[d]["pairs"] is a dictionary with keys for each pair type linked to arrays of the instances of pair separations for that type)
        E_LJ[count_d] = sum( ( 4.0*eps_array*sigma_array**(12.0)*Data[d]["pair_vector-12"] - 4.0*eps_array*sigma_array**(6.0)*Data[d]["pair_vector-6"] ) )

    return E_LJ

# Returns the BUCK energy for each configuration (as an array) using the parameters currently defined in VDW_dict
# and the supplied parameters for the type being fit (fit_type is globally defined outside of the function). 
def E_BUCK_Fit(ind,A,B,C):

    # In python you need to expand the scope of global variables at every level
    global VDW_dict,fit_type,Pair_min_dict

    # If the attempted value for the fit type violates the sanity checks on the fit range then the energy is set to an astronomical value
    # NOTE: the restriction on sigma_min is to stop the minimum of the fit potential from lieing outside of the sampled configuration space
    if fit_type in list(Pair_min_dict.keys()) and ( ( A < 0.001 or B < 0.001 or C < 0.001 ) or True not in [ (A*exp(-i/B) - C*i**(-6.0) ) > 10.0 for i in np.arange(0.1,5.1,0.01) ] ):
        return np.ones(len(ind))*1.E20        

    # Initialize E_BUCK for the current configuration
    E_BUCK = np.zeros(len(ind))

    # Iterate over the configuration index
    for count_d,d in enumerate(ind):

        # Iterate over list of pair-seps for each pair-type (Data[d]["pairs"] is a dictionary with keys for each pair type linked to arrays of the instances of pair separations for that type)
        for i in list(Data[d]["pairs"].keys()):
            if i == fit_type or (i[1],i[0]) == fit_type:
                E_BUCK[count_d] += sum( ( A*exp(-Data[d]["pairs"][i]/B) - C*Data[d]["pairs"][i]**(-6.0) ) )
            else:
                E_BUCK[count_d] += sum( ( VDW_dict[i][0]*exp(-Data[d]["pairs"][i]/VDW_dict[i][1]) - VDW_dict[i][2]*Data[d]["pairs"][i]**(-6.0) ) )

    return E_BUCK

def detect_outlier(E_config):
    # modified z-score: Mi=0.6745*(xi-x_m)/MAD 
    # x_m: median; 
    #MAD: mean absolute deviation = median(|xi-x_m|)
    # recommanded threshold: remove Mi>3.5 
    threshold = 6
    
    MAD = stats.median_absolute_deviation(E_config)
    x_m = np.median(E_config)
    # Mi holds the corresponding modified z-score
    # outliter is a bunch of true/false indicating whether that point is outlier
    Mi = np.zeros(len(E_config))
    outlier = np.zeros(len(E_config))
    for count_i,i in enumerate(E_config):
        #Mi[count_i] = abs(0.6745*(i-x_m)/MAD)
        Mi[count_i] = 0.6745*(i-x_m)/MAD
        if Mi[count_i] < -threshold: outlier[count_i] = True
        #if Mi[count_i] > threshold : outlier[count_i] = True
        else: outlier[count_i] = False
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    plot_handle, = ax.plot(E_config,Mi,marker='.',markersize=20,color='r',markeredgewidth=0.0,alpha=0.3,linestyle='None')
    ax.set_xlabel('E_config(kcal/mol)')
    ax.set_ylabel('modified z-score')
    ax.tick_params(axis='both', which='major',labelsize=12,pad=10,direction='out',width=2,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=12,pad=10,direction='out',width=2,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]
    # Save the figure
    # XXX: the plot path need to be changed
    plot_name = 'z-score'
    plt.savefig(plot_name, dpi=300)
    plt.close(fig)

    return Mi,outlier

# fit_type is a globally defined tuple
def global_fit_LJ(e_0=0, s_0=0, e_1=0,  s_1=0,  e_2=0,  s_2=0,  e_3=0,  s_3=0,  e_4=0,  s_4=0,  e_5=0,  s_5=0,  e_6=0,  s_6=0,  e_7=0,  s_7=0,  e_8=0,  s_8=0,  e_9=0,  s_9=0,  e_10=0, s_10=0,\
                                e_11=0, s_11=0, e_12=0, s_12=0, e_13=0, s_13=0, e_14=0, s_14=0, e_15=0, s_15=0, e_16=0, s_16=0, e_17=0, s_17=0, e_18=0, s_18=0, e_19=0, s_19=0, e_20=0, s_20=0,\
                                e_21=0, s_21=0, e_22=0, s_22=0, e_23=0, s_23=0, e_24=0, s_24=0, e_25=0, s_25=0, e_26=0, s_26=0, e_27=0, s_27=0, e_28=0, s_28=0, e_29=0, s_29=0, e_30=0, s_30=0,\
                                e_31=0, s_31=0, e_32=0, s_32=0, e_33=0, s_33=0, e_34=0, s_34=0, e_35=0, s_35=0, e_36=0, s_36=0, e_37=0, s_37=0, e_38=0, s_38=0, e_39=0, s_39=0, e_40=0, s_40=0,\
                                e_41=0, s_41=0, e_42=0, s_42=0, e_43=0, s_43=0, e_44=0, s_44=0, e_45=0, s_45=0, e_46=0, s_46=0, e_47=0, s_47=0, e_48=0, s_48=0, e_49=0, s_49=0, e_50=0, s_50=0,\
                                e_51=0, s_51=0, e_52=0, s_52=0, e_53=0, s_53=0, e_54=0, s_54=0, e_55=0, s_55=0, e_56=0, s_56=0, e_57=0, s_57=0, e_58=0, s_58=0, e_59=0, s_59=0, e_60=0, s_60=0,\
                                e_61=0, s_61=0, e_62=0, s_62=0, e_63=0, s_63=0, e_64=0, s_64=0, e_65=0, s_65=0, e_66=0, s_66=0, e_67=0, s_67=0, e_68=0, s_68=0, e_69=0, s_69=0, e_70=0, s_70=0,\
                                e_71=0, s_71=0, e_72=0, s_72=0, e_73=0, s_73=0, e_74=0, s_74=0, e_75=0, s_75=0, e_76=0, s_76=0, e_77=0, s_77=0, e_78=0, s_78=0, e_79=0, s_79=0, e_80=0, s_80=0,\
                                e_81=0, s_81=0, e_82=0, s_82=0, e_83=0, s_83=0, e_84=0, s_84=0, e_85=0, s_85=0, e_86=0, s_86=0, e_87=0, s_87=0, e_88=0, s_88=0, e_89=0, s_89=0, e_90=0, s_90=0,\
                                e_91=0, s_91=0, e_92=0, s_92=0, e_93=0, s_93=0, e_94=0, s_94=0, e_95=0, s_95=0, e_96=0, s_96=0, e_97=0, s_97=0, e_98=0, s_98=0, e_99=0, s_99=0, e_100=0,s_100=0,\
                                e_101=0,s_101=0,e_102=0,s_102=0,e_103=0,s_103=0,e_104=0,s_104=0,e_105=0,s_105=0,e_106=0,s_106=0,e_107=0,s_107=0,e_108=0,s_108=0,e_109=0,s_109=0,e_110=0,s_110=0,\
                                e_111=0,s_111=0,e_112=0,s_112=0,e_113=0,s_113=0,e_114=0,s_114=0,e_115=0,s_115=0,e_116=0,s_116=0,e_117=0,s_117=0,e_118=0,s_118=0,e_119=0,s_119=0,e_120=0,s_120=0,\
                                e_121=0,s_121=0,e_122=0,s_122=0,e_123=0,s_123=0,e_124=0,s_124=0,e_125=0,s_125=0,e_126=0,s_126=0,e_127=0,s_127=0,e_128=0,s_128=0,e_129=0,s_129=0,e_130=0,s_130=0,\
                                e_131=0,s_131=0,e_132=0,s_132=0,e_133=0,s_133=0,e_134=0,s_134=0,e_135=0,s_135=0,e_136=0,s_136=0,e_137=0,s_137=0,e_138=0,s_138=0,e_139=0,s_139=0,e_140=0,s_140=0,\
                                e_141=0,s_141=0,e_142=0,s_142=0,e_143=0,s_143=0,e_144=0,s_144=0,e_145=0,s_145=0,e_146=0,s_146=0,e_147=0,s_147=0,e_148=0,s_148=0,e_149=0,s_149=0,e_150=0,s_150=0,\
                                e_151=0,s_151=0,e_152=0,s_152=0,e_153=0,s_153=0,e_154=0,s_154=0,e_155=0,s_155=0,e_156=0,s_156=0,e_157=0,s_157=0,e_158=0,s_158=0,e_159=0,s_159=0,e_160=0,s_160=0,\
                                e_161=0,s_161=0,e_162=0,s_162=0,e_163=0,s_163=0,e_164=0,s_164=0,e_165=0,s_165=0,e_166=0,s_166=0,e_167=0,s_167=0,e_168=0,s_168=0,e_169=0,s_169=0,e_170=0,s_170=0,\
                                e_171=0,s_171=0,e_172=0,s_172=0,e_173=0,s_173=0,e_174=0,s_174=0,e_175=0,s_175=0,e_176=0,s_176=0,e_177=0,s_177=0,e_178=0,s_178=0,e_179=0,s_179=0,e_180=0,s_180=0,\
                                e_181=0,s_181=0,e_182=0,s_182=0,e_183=0,s_183=0,e_184=0,s_184=0,e_185=0,s_185=0,e_186=0,s_186=0,e_187=0,s_187=0,e_188=0,s_188=0,e_189=0,s_189=0,e_190=0,s_190=0,\
                                e_191=0,s_191=0,e_192=0,s_192=0,e_193=0,s_193=0,e_194=0,s_194=0,e_195=0,s_195=0,e_196=0,s_196=0,e_197=0,s_197=0,e_198=0,s_198=0,e_199=0,s_199=0,e_200=0,s_200=0,\
                                ind=[],outlier=[],E_config=[],VDW_dict=[],Data=[],Fit_Pairs=[],Pair_min_dict=[],penalty=100.0,mixing_rule="none",L2_sigma=0.0,L2_eps=0.0,VDW_0=None):

    # Initialize local variable dictionary (used for determining what has been defined)
    local_vars = locals()

    # Find defined e and s variables. *vars holds the variable names *vals holds the variable values
    e_vals = np.array([ local_vars[j] for j in natural_sort([ i for i in local_vars if 'e_' in i and "__" not in i and local_vars[i] != 0 ]) ])
    s_vals = np.array([ local_vars[j] for j in natural_sort([ i for i in local_vars if 's_' in i and "__" not in i and local_vars[i] != 0 ]) ])

    # Consistency check
    if len(e_vals) != len(Fit_Pairs):
        print("ERROR in global_fit_LJ: the function expects the len(Fit_Pairs) and number of non-zero e_* variable to be equal. Exiting...")
        quit()

    # Initialize eps and sigma arrays for the vectorized calculations
    eps_array = np.array([ VDW_dict[i][0] for i in Data[ind[0]]["pair_type_vector"] ])
    sigma_array = np.array([ VDW_dict[i][1] for i in Data[ind[0]]["pair_type_vector"] ])

    # Update the elements being fit with the corresponding e_vals and s_vals
    for count_c,c in enumerate(Fit_Pairs):        
        for i in [ count_j for count_j,j in enumerate(Data[ind[0]]["pair_type_vector"]) if ( j == c or (j[1],j[0]) == c ) ]:        
            eps_array[i]   = e_vals[count_c]
            sigma_array[i] = s_vals[count_c]

    # Apply mixing rules
    # If/else conditions handle cases where one or both of the involved types is in the fit set 
    if mixing_rule == "wh":
        fit_types = { i[0]:count_i for count_i,i in enumerate(Fit_Pairs) }
        for count_i,i in [ (count_j,j) for count_j,j in enumerate(Data[ind[0]]["pair_type_vector"]) if ( j[0] in fit_types or j[1] in fit_types )]:            
            if i[0] in fit_types:
                if i[1] in fit_types:
                    sigma_array[count_i] = (( s_vals[fit_types[i[0]]]**(6.0) + s_vals[fit_types[i[1]]]**(6.0) ) / 2.0 )**(1.0/6.0)
                    eps_array[count_i]   = (e_vals[fit_types[i[0]]]*s_vals[fit_types[i[0]]]**(6.0) * e_vals[fit_types[i[1]]]*s_vals[fit_types[i[1]]]**(6.0) )**(0.5) / sigma_array[count_i]**(6.0)
                else:
                    sigma_array[count_i] = (( s_vals[fit_types[i[0]]]**(6.0) + VDW_dict[(i[1],i[1])][1]**(6.0) ) / 2.0 )**(1.0/6.0)
                    eps_array[count_i]   = (e_vals[fit_types[i[0]]]*s_vals[fit_types[i[0]]]**(6.0) * VDW_dict[(i[1],i[1])][0]*VDW_dict[(i[1],i[1])][1]**(6.0) )**(0.5) / sigma_array[count_i]**(6.0)
            elif i[1] in fit_types:
                sigma_array[count_i] = (( s_vals[fit_types[i[1]]]**(6.0) + VDW_dict[(i[0],i[0])][1]**(6.0) ) / 2.0 )**(1.0/6.0)
                eps_array[count_i]   = (e_vals[fit_types[i[1]]]*s_vals[fit_types[i[1]]]**(6.0) * VDW_dict[(i[0],i[0])][0]*VDW_dict[(i[0],i[0])][1]**(6.0) )**(0.5) / sigma_array[count_i]**(6.0)

    elif mixing_rule == "lb":
        fit_types = { i[0]:count_i for count_i,i in enumerate(Fit_Pairs) }

        for count_i,i in [ (count_j,j) for count_j,j in enumerate(Data[ind[0]]["pair_type_vector"]) if ( j[0] in fit_types or j[1] in fit_types )]:            
            if i[0] in fit_types:
                if i[1] in fit_types:
                    sigma_array[count_i] = (s_vals[fit_types[i[0]]] + s_vals[fit_types[i[1]]] ) / 2.0 
                    eps_array[count_i]   = (e_vals[fit_types[i[0]]] * e_vals[fit_types[i[1]]] )**(0.5)
                else:
                    sigma_array[count_i] = (s_vals[fit_types[i[0]]] + VDW_dict[(i[1],i[1])][1] ) / 2.0 
                    eps_array[count_i]   = (e_vals[fit_types[i[0]]] * VDW_dict[(i[1],i[1])][0] )**(0.5)
            elif i[1] in fit_types:
                sigma_array[count_i] = (VDW_dict[(i[0],i[0])][1] + VDW_dict[(i[1],i[1])][1] ) / 2.0 
                eps_array[count_i]   = (VDW_dict[(i[0],i[0])][0] * VDW_dict[(i[1],i[1])][0] )**(0.5)

    # Initialize L2 regularization sums
    if VDW_0 is not None:

        # only loop over the values in the sigma and eps arrays (they are the only ones that affect the fit!)
        L2_e = np.mean([ (e_vals[count_i]-VDW_0[i][0])**(2.0) for count_i,i in enumerate(Fit_Pairs) if i[0] == i[1] ]) 
        L2_s = np.mean([ (s_vals[count_i]-VDW_0[i][1])**(2.0) for count_i,i in enumerate(Fit_Pairs) if i[0] == i[1] ]) 

    else:
        L2_e = np.mean(eps_array**(2.0))
        L2_s = np.mean(sigma_array**(2.0))
        
    # Initialize E_LJ for the current configuration
    xhi2 = 0

    # Iterate over the configuration index
    data_count = 0
    for count_d,d in enumerate(ind):
       
        if len(outlier) != 0: 
            if outlier[count_d]: continue
        data_count += 1
        # Iterate over list of pair-seps for each pair-type (Data[d]["pairs"] is a dictionary with keys for each pair type linked to arrays of the instances of pair separations for that type)
        xhi2 += ( sum( ( 4.0*eps_array*sigma_array**(12.0)*Data[d]["pair_vector-12"] - 4.0*eps_array*sigma_array**(6.0)*Data[d]["pair_vector-6"] ) ) - E_config[count_d] )**(2.0)

    return xhi2/float(data_count) + L2_sigma*L2_s + L2_eps*L2_e

# This function expects data to be a globally defined dictionary. The function
# outputs the LJ contribution to the interaction energy of a pair of molecules
def E_LJ_Tot(ind):

    global Data,VDW_dict

    # Initialize the energy
    E_LJ = 0
    
    # Iterate over list of pair-seps for each pair-type (Data[d]["pairs"] is a dictionary with keys for each pair type linked to arrays of the instances of pair separations for that type)
    for i in list(Data[ind]["pairs"].keys()):
        E_LJ += sum( ( 4.0*VDW_dict[i][0]*VDW_dict[i][1]**(12.0)*Data[ind]["pairs"][i]**(-12.0) - 4.0*VDW_dict[i][0]*VDW_dict[i][1]**(6.0)*Data[ind]["pairs"][i]**(-6.0) ) )

    return E_LJ

# This function expects data to be a globally defined dictionary. The function
# outputs the Buckingham contribution to the interaction energy of a pair of molecules
def E_BUCK_Tot(ind):

    global Data,VDW_dict

    # Initialize the energy
    E_BUCK = 0
    
    # Iterate over list of pair-seps for each pair-type (Data[d]["pairs"] is a dictionary with keys for each pair type linked to arrays of the instances of pair separations for that type)
    for i in list(Data[ind]["pairs"].keys()):
        E_BUCK += sum( ( VDW_dict[i][0]*exp(-Data[ind]["pairs"][i]/VDW_dict[i][1]) - VDW_dict[i][2]*Data[ind]["pairs"][i]**(-6.0) ) )

    return E_BUCK

# Description: A simple wrapper for the file writing commands for the vdw parameters.
#              The function takes arguments for the filename, fit type (e.g. lj), and
#              united-atom (UA) vs all-atom (AA) style.
def write_params(name,fit,eng_type='AA',write_charges=0):

    # In python you need to expand the scope of global variables at every level
    global VDW_dict,Data,min_pairs

    # Collect minimum types
    min_types = []
    for i in min_pairs:
        for j in i:
            min_types += [j,j+'-UA']
    min_types = set(min_types)
    
    # First clean up VDW_dict by removing reverse types with the wrong priority ordering
    # NOTE: parameters are saved using type_a > type_b ordering
    # NOTE: Hydrogen containing pairs are removed when the fit_type is 'UA'
    pop_list = []
    keys = list(VDW_dict.keys())
    for i in keys:
        if i[0] > i[1]:
            if (i[1],i[0]) in list(VDW_dict.keys()): VDW_dict.pop((i[1],i[0]))
        if eng_type=="UA" and ( return_UA_H(i[0]) == 1 or return_UA_H(i[1]) == 1 ):
            if i in list(VDW_dict.keys()): VDW_dict.pop(i)
            if (i[1],i[0]) in list(VDW_dict.keys()): VDW_dict.pop((i[1],i[0]))

    # Save params
    with open(name+'.db','w') as f:

        # Write VDW definitions
        if fit == 'lj':
            f.write("\n# VDW definitions\n#\n{:10s} {:60s} {:60s} {:11s} {:20s} {:20s}\n".format("#","Atom_type","Atom_type","Potential","eps (kcal/mol)","sigma (ang)"))        
        if fit == 'buck':
            f.write("\n# VDW definitions\n#\n{:10s} {:60s} {:60s} {:11s} {:20s} {:20s} {:20s}\n".format("#","Atom_type","Atom_type","Potential","A (kcal/mol)","B (ang)","C (ang kcal/mol)"))        
        if eng_type == 'AA':
            for i in list(VDW_dict.keys()):
                f.write("{:10s} {:60s} {:60s} {:10s} {}\n".format("vdw",str(i[0]),str(i[1]),fit," ".join([ "{:< 20.6f}".format(j) for j in VDW_dict[i]])))        
        elif eng_type == 'UA':
            for i in list(VDW_dict.keys()):
                f.write("{:10s} {:60s} {:60s} {:10s} {}\n".format("vdw",str(i[0])+'-UA',str(i[1])+'-UA',fit," ".join([ "{:< 20.6f}".format(j) for j in VDW_dict[i]])))        

        # Write charges is the flag is enabled
        if write_charges == 1:

            # Collect AA-charges and UA-charges
            atom_types_AA = []
            atom_types_UA = []
            charges_AA    = []
            charges_UA    = []
            for i in list(Data.keys()):
                for count_j,j in enumerate(Data[i]["types_a"]):
                    if j not in atom_types_AA:
                        atom_types_AA += [j]
                        charges_AA    += [Data[i]["charges_a"][count_j]]
                        if return_UA_H(j) != 1:
                            atom_types_UA += [j+'-UA']
                            charges_UA    += [Data[i]["charges_a_UA"][count_j]]
                for count_j,j in enumerate(Data[i]["types_b"]):
                    if j not in atom_types_AA:
                        atom_types_AA += [j]
                        charges_AA    += [Data[i]["charges_b"][count_j]]
                        if return_UA_H(j) != 1:
                            atom_types_UA += [j+'-UA']
                            charges_UA    += [Data[i]["charges_b_UA"][count_j]]
                    
            # Write charges
            f.write("\n# Charge definitions\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","Charge"))
            if eng_type == 'AA':
                for count_i,i in enumerate(atom_types_AA):
                    f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(i),charges_AA[count_i]))
            if eng_type == 'UA':
                for count_i,i in enumerate(atom_types_UA):
                    f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(i),charges_UA[count_i]))

    # Save minimum pair info
    with open(name+'-min.db','w') as f:

        # Write VDW definitions
        if fit == 'lj':
            f.write("\n# VDW definitions\n#\n{:10s} {:60s} {:60s} {:11s} {:20s} {:20s}\n".format("#","Atom_type","Atom_type","Potential","eps (kcal/mol)","sigma (ang)"))        
        if fit == 'buck':
            f.write("\n# VDW definitions\n#\n{:10s} {:60s} {:60s} {:11s} {:20s} {:20s} {:20s}\n".format("#","Atom_type","Atom_type","Potential","A (kcal/mol)","B (ang)","C (ang kcal/mol)"))        
        if eng_type == 'AA':
            for i in [ j for j in list(VDW_dict.keys()) if j in min_pairs ]:
                f.write("{:10s} {:60s} {:60s} {:10s} {}\n".format("vdw",str(i[0]),str(i[1]),fit," ".join([ "{:< 20.6f}".format(j) for j in VDW_dict[i]])))        
        elif eng_type == 'UA':
            for i in [ j for j in list(VDW_dict.keys()) if j in min_pairs ]:
                f.write("{:10s} {:60s} {:60s} {:10s} {}\n".format("vdw",str(i[0])+'-UA',str(i[1])+'-UA',fit," ".join([ "{:< 20.6f}".format(j) for j in VDW_dict[i]])))        

        # Write charges is the flag is enabled
        if write_charges == 1:

            # Collect AA-charges and UA-charges
            atom_types_AA = []
            atom_types_UA = []
            charges_AA    = []
            charges_UA    = []
            for i in list(Data.keys()):
                for count_j,j in enumerate(Data[i]["types_a"]):
                    if j not in atom_types_AA:
                        atom_types_AA += [j]
                        charges_AA    += [Data[i]["charges_a"][count_j]]
                        if return_UA_H(j) != 1:
                            atom_types_UA += [j+'-UA']
                            charges_UA    += [Data[i]["charges_a_UA"][count_j]]
                for count_j,j in enumerate(Data[i]["types_b"]):
                    if j not in atom_types_AA:
                        atom_types_AA += [j]
                        charges_AA    += [Data[i]["charges_b"][count_j]]
                        if return_UA_H(j) != 1:
                            atom_types_UA += [j+'-UA']
                            charges_UA    += [Data[i]["charges_b_UA"][count_j]]
                    
            # Write charges
            f.write("\n# Charge definitions\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","Charge"))
            if eng_type == 'AA':
                for count_i,i in enumerate(atom_types_AA):
                    if i not in min_types:
                        continue
                    else:
                        f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(i),charges_AA[count_i]))
            if eng_type == 'UA':
                for count_i,i in enumerate(atom_types_UA):
                    if i not in min_types:
                        continue
                    else:
                        f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(i),charges_UA[count_i]))

    return

def plot_fit_potentials(folder,fit='lj',eng_type='DFT-AA',r_min=0.1,r_max=10.0):

    # In python you need to expand the scope of global variables at every level
    global VDW_dict

    # Initialize the rij range
    seps = np.arange(r_min,r_max,0.01)

    # Plot all pair potentials
    # NOTE: write_params call cleans up VDW_dict.keys to remove redundant pair_types 
    for i in list(VDW_dict.keys()):

        if fit == 'lj':
            y_vals = 4*VDW_dict[i][0]*( (VDW_dict[i][1]/seps)**12 - (VDW_dict[i][1]/seps)**6 )
        elif fit == 'buck':
            y_vals = VDW_dict[i][0]*exp(-seps/VDW_dict[i][1]) - VDW_dict[i][2]*seps**(-6.0) 

        # Initialize figure and generate plot
        fig = plt.figure(figsize=(6,4))
        ax = plt.subplot(111)
        ax.plot(seps,y_vals,color=(0.05,0.35,0.75),linestyle='-',linewidth=2.0)
        
        # Set limits based on the largest range
        y_min = -5.0
        if fit == 'lj':
            while (1):
                if y_min > min(y_vals):
                    y_min -= 5
                else:
                    break
        ax.set_xlim([0,r_max])
        ax.set_ylim([y_min,5])

        # Set Labels and Save Name
        ax.set_ylabel("$\mathrm{E_{LJ} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        ax.set_xlabel("$r\mathrm{_{ij} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')

        # Set save name
        if eng_type == 'DFT-AA':
            plot_name = i[0]+'-'+i[1]+'_DFT-AA.png'
        elif eng_type == 'DFT-UA':
            plot_name = i[0]+'-'+i[1]+'_DFT-UA.png'
        elif eng_type == 'MP2-AA':
            plot_name = i[0]+'-'+i[1]+'_MP2-AA.png'
        elif eng_type == 'MP2-UA':
            plot_name = i[0]+'-'+i[1]+'_MP2-UA.png'

        # Format ticks and plot box
        ax.tick_params(axis='both', which='major',labelsize=12,pad=10,direction='out',width=2,length=4)
        ax.tick_params(axis='both', which='minor',labelsize=12,pad=10,direction='out',width=2,length=3)
        [j.set_linewidth(2) for j in ax.spines.values()]

        # Save the figure
        plt.savefig(folder+'/'+plot_name, dpi=300, bbox_inches='tight')
        plt.close(fig)

    return

# Description: This function is a simple wrapper for the matplotlib plotting commands
#              to generate the scatterplot of the correlation between the QC interaction
#              energies and the FF interaction energies. The function also saves the 
#              convergence data to file.
def plot_convergence(folder,Data,fit_vals,plot_num=5,eng_type='DFT-AA'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []

    # Generate subset of plot values ("plot_num" argument determines the number of convergence cycles
    # to explicitly plot, although all are saved to the .txt file; e.g. plot_num plots the first and last
    # cycles plus the 1/4 1/2 and 3/4 converged cycles).
    if len(fit_vals) <= plot_num:
        ind = list(range(0,len(fit_vals)))
    else:
        ind = [ int(round((len(fit_vals)-1)/float(plot_num-1)*(i))) for i in range(plot_num) ]

    #
    # Plot VDW_fit correlation plot
    #
    x_vals = [ i for i in list(Data.keys()) ]

    # Set the y_vals to the total dft energy
    if eng_type == 'DFT-AA':
        y_vals2 = [ Data[i]["e_dft_fit_AA_ori"] for i in list(Data.keys()) ]    
        y_vals = [ Data[i]["e_dft_fit_AA"] for i in list(Data.keys()) ]    
    if eng_type == 'DFT-UA':
        y_vals = [ Data[i]["e_dft_fit_UA"] for i in list(Data.keys()) ]    
    if eng_type == 'MP2-AA':
        y_vals = [ Data[i]["e_mp2_fit_AA"] for i in list(Data.keys()) ]    
    if eng_type == 'MP2-UA':
        y_vals = [ Data[i]["e_mp2_fit_UA"] for i in list(Data.keys()) ]    

    # Plot the scatter data
    for count_i,i in enumerate(ind):   
        legend = ['initial_guess','Final fitting']
        if count_i == 0: plot_handle, = ax.plot(y_vals2,fit_vals[i],marker='.',markersize=20,color=color_list[count_i],markeredgewidth=0.0,alpha=0.3,linestyle='None',label=legend[count_i])
        if count_i == 1: plot_handle, = ax.plot(y_vals,fit_vals[i],marker='.',markersize=20,color=color_list[count_i],markeredgewidth=0.0,alpha=0.3,linestyle='None',label=legend[count_i])
        #plot_handle, = ax.plot(y_vals,fit_vals[i],marker='.',markersize=20,color=color_list[count_i],markeredgewidth=0.0,alpha=0.3,linestyle='None',label=legend[count_i])

    # Find min and max based on the final fit_vals
    y_min_fit = 1000000.0
    y_max_fit =-1000000.0
    for i in fit_vals:
        if min(i) < y_min_fit: y_min_fit = np.floor(min(i))
        if max(i) > y_min_fit: y_max_fit = np.ceil(max(i))        
            
    # Set limits based on the largest range
    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    # If y_min or y_max are impractically larger than the limits of the final fit values, then the y limits are reassigned
    if abs(y_min-y_min_fit) > 10: y_min = y_min_fit
    if abs(y_max-y_max_fit) > 10: y_max = y_max_fit

    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")


    ####### I found a way to implement .format method in mathrm, please see plot_convergence() function in Parser/plot_vdw_convergence

    # Set Labels and Save Name
    ax.set_ylabel("$\mathrm{E_{FF,VDW} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
    if eng_type == 'DFT-AA':
        ax.set_xlabel("$\mathrm{E_{DFT,VDW} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'DFT-AA_vdw_convergence.png'
        data_name = 'DFT-AA_vdw_convergence.txt'
    elif eng_type == 'DFT-UA':
        ax.set_xlabel("$\mathrm{E_{DFT,VDW} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'DFT-UA_vdw_convergence.png'
        data_name = 'DFT-UA_vdw_convergence.txt'
    elif eng_type == 'MP2-AA':
        ax.set_xlabel("$\mathrm{E_{MP2,VDW} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'MP2-AA_vdw_convergence.png'
        data_name = 'MP2-AA_vdw_convergence.txt'
    elif eng_type == 'MP2-UA':
        ax.set_xlabel("$\mathrm{E_{MP2,VDW} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'MP2-UA_vdw_convergence.png'
        data_name = 'MP2-UA_vdw_convergence.txt'

    # Generate Legend
    ax.legend(loc='best',frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=12,pad=10,direction='out',width=2,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=12,pad=10,direction='out',width=2,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(folder+'/'+plot_name, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    # Save the convergence data to file
    with open(folder+'/'+data_name,'w') as f:
        if eng_type in ['DFT-AA','DFT-UA']:
            #f.write(' {:<15s} {}\n'.format("E_INT_VDW_DFT",' '.join([ "{:<15s}".format("E_INT_VDW_FF_"+str(count_i)) for count_i,i in enumerate(fit_vals) ])))
            f.write(' {:<15s} {:<15s} {}\n'.format("E_INT_VDW_DFT","E_VDW_DFT_ori",' '.join([ "{:<15s}".format("E_INT_VDW_FF_"+str(count_i)) for count_i,i in enumerate(fit_vals) ])))
        elif eng_type in ['MP2-AA','MP2-UA']:
            f.write(' {:<15s} {}\n'.format("E_INT_VDW_MP2",' '.join([ "{:<15s}".format("E_INT_VDW_FF_"+str(count_i)) for count_i,i in enumerate(fit_vals) ])))
        for count_i,i in enumerate(y_vals):
            #f.write('{:< 15.8f} {}\n'.format(i,' '.join([ "{:< 15.8f}".format(j[count_i]) for j in fit_vals ])))
            f.write('{:< 15.8f} {:< 15.8f} {}\n'.format(i,y_vals2[count_i],' '.join([ "{:< 15.8f}".format(j[count_i]) for j in fit_vals ])))

    #
    # Plot Total_fit correlation plot (Total interaction energy = VDW + Coul contribution
    #
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []

    # Initialize list of indices
    x_vals = [ i for i in list(Data.keys()) ]

    # Set the y_vals to the total dft energy
    if eng_type == 'DFT-AA' or eng_type == 'DFT-UA':
        y_vals = [ Data[i]["e_dft"] for i in list(Data.keys()) ]    

    # Set the y_vals to the total mp2 energy
    if eng_type == 'MP2-AA' or eng_type == 'MP2-UA':
        y_vals = [ Data[i]["e_mp2"] for i in list(Data.keys()) ]    
    
    # Create all-atom/united-atom coulombic array
    if eng_type == 'DFT-AA' or eng_type == 'MP2-AA':
        coul_array = np.array([ E_C_Tot_drude(i,'AA') for i in x_vals ])
        coul_array_ori = np.array([ E_C_Tot(i,'AA') for i in x_vals ])
    if eng_type == 'DFT-UA' or eng_type == 'MP2-UA':
        coul_array = np.array([ E_C_Tot(i,'UA') for i in x_vals ])

    # Add coulombic contribution to the fit values
    for i in range(len(fit_vals)):
        #fit_vals[i] += coul_array 
        if i == 0: fit_vals[i] += coul_array_ori
        else:fit_vals[i] += coul_array

    # Plot the scatter data
    for count_i,i in enumerate(ind):   
        plot_handle, = ax.plot(y_vals,fit_vals[i],marker='.',markersize=20,color=color_list[count_i],markeredgewidth=0.0,alpha=0.3,linestyle='None',label='cycle: {}'.format(i))

    # Find min and max based on the final fit_vals
    y_min_fit = 1000000.0
    y_max_fit =-1000000.0
    for i in fit_vals:
        if min(i) < y_min_fit: y_min_fit = np.floor(min(i))
        if max(i) > y_min_fit: y_max_fit = np.ceil(max(i))        
            
    # Set limits based on the largest range
    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    # If y_min or y_max are impractically larger than the limits of the final fit values, then the y limits are reassigned
    if abs(y_min-y_min_fit) > 10: y_min = y_min_fit
    if abs(y_max-y_max_fit) > 10: y_max = y_max_fit

    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")

    # Set Labels and Save Name
    ax.set_ylabel("$\mathrm{E_{FF,TOT} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
    if eng_type == 'DFT-AA':
        ax.set_xlabel("$\mathrm{E_{DFT,TOT} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'DFT-AA_tot_convergence.png'
        data_name = 'DFT-AA_tot_convergence.txt'
    elif eng_type == 'DFT-UA':
        ax.set_xlabel("$\mathrm{E_{DFT,TOT} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'DFT-UA_tot_convergence.png'
        data_name = 'DFT-UA_tot_convergence.txt'
    elif eng_type == 'MP2-AA':
        ax.set_xlabel("$\mathrm{E_{MP2,TOT} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'MP2-AA_tot_convergence.png'
        data_name = 'MP2-AA_tot_convergence.txt'
    elif eng_type == 'MP2-UA':
        ax.set_xlabel("$\mathrm{E_{MP2,TOT} \, (kcal \cdot mol^{-1})}$",fontsize=32,labelpad=10,fontweight='bold')
        plot_name = 'MP2-UA_tot_convergence.png'
        data_name = 'MP2-UA_tot_convergence.txt'

    # Generate Legend
    ax.legend(loc='best',frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=12,pad=10,direction='out',width=2,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=12,pad=10,direction='out',width=2,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(folder+'/'+plot_name, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    # Save the convergence data to file
    with open(folder+'/'+data_name,'w') as f:
        if eng_type in ['DFT-AA','DFT-UA']:
            f.write(' {:<15s} {}\n'.format("E_INT_TOT_DFT",' '.join([ "{:<15s}".format("E_INT_TOT_FF_"+str(count_i)) for count_i,i in enumerate(fit_vals) ])))
        elif eng_type in ['MP2-AA','MP2-UA']:
            f.write(' {:<15s} {}\n'.format("E_INT_TOT_MP2",' '.join([ "{:<15s}".format("E_INT_TOT_FF_"+str(count_i)) for count_i,i in enumerate(fit_vals) ])))
        for count_i,i in enumerate(y_vals):
            f.write('{:< 15.8f} {}\n'.format(i,' '.join([ "{:< 15.8f}".format(j[count_i]) for j in fit_vals ])))

    return

# Create logger to save stdout to logfile
class Logger(object):

    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/extract_vdw.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

if __name__ == "__main__":
   main(sys.argv[1:])
