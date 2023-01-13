#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,getpass,json,fnmatch

# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
from monitor_jobs import *
from parse_all import *
import frag_gen_inter,frag_gen,xyz_to_orca,paramgen,merge_FF,extract_intramolecular_params,gen_md_for_sampling,gen_jobs_for_charges 
# Old version TAFFI
import frag_gen_old,xyz_to_orca_old,paramgen_old,extract_intramolecular_params_old
import restart_scans,dihedral_restart
import codecs,json


def main(argv):
    user = getpass.getuser()
    run_dir = os.getcwd()

    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-o',dest='output',default='taffi_intra.log',
                        help = 'log name for stdout')

    args=parser.parse_args()    
    sys.stdout = Logger(run_dir+'/'+args.output)
    
    # parse configuration dictionary (c)
    print("parsing configurations...")
    c = parse_configuration(args.config)
    # Make config file absolute
    args.config = os.path.abspath(args.config)
    if (c['module_string'] == None):
      c['module_string'] = ''
    

    print("Working on final intramolecular parameters...")  

    # If local parameters are available then use them for the FF argument
    local_FF = run_dir+'/Params_for_Batch.db'
    if os.path.isfile(local_FF) :
         c['ff'] = local_FF

##############################
# Run the frag optimizations #
##############################

    # Generate Intermolecular Modes Fragments
    if os.path.isdir('Intramolecular_Modes') is False:

         # Assemble shell submit string for intramolecular fragment generation
         # On some systems running simple python scripts on the head node isn't allowed, so it is safer to just always run on compute nodes
         jobids = []
         sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
         sublist.append('python3 {}/FF_functions/frag_gen_old.py \'*.xyz\' -FF \'{}\' -gens {} -q {} -o Intramolecular_Modes'.format(c["taffi_path"],c["ff"],c['gens'],c['charge']))
         sublist= sublist + ('frag_gen_intra  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
         output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
         output = str(output,'utf-8')
         jobids += [ output.split("\n")[-2].split()[-1]]

         # Wait until jobs complete
         monitor_jobs(jobids,user)

         # QC geometry optimize starting fragments      
         os.chdir('Intramolecular_Modes')
         subprocess.Popen('chmod 777 mass_geoopt.sh'.split()) 

         # Check if there are any intramolecular modes in need of parameterization
         xyzfiles = [f for f in os.listdir('.') if os.path.isfile(f) and f.split('.')[-1] == 'xyz']
         if xyzfiles == []:
            print("No intramolecular jobs in need of execution. Exiting...")
            quit()

         # Generate the geometry optimization files
         for i in xyzfiles:
            name = i.split('.')[0]
            if(c['functional'] == 'wB97X-D3'):
               command_list = ('{} -p {} -r 0.0 --no_freq -q {} -f {} -b {} --no_D3 -o {}_geoopt.in'.format(i,c['param_geoopt_procs'],c['charge'],c['functional'],c['basis'],name)).split()
            else:
               command_list = ('{} -p {} -r 0.0 --no_freq -q {} -f {} -b {}  -o {}_geoopt.in'.format(i,c['param_geoopt_procs'],c['charge'],c['functional'],c['basis'],name)).split()
            xyz_to_orca_old.main(command_list)
            os.makedirs(name+'_geoopt')
            shutil.move(name+'_geoopt.in',name+'_geoopt')
         
         # Submit geometry optimizations
         jobids = []
         substring = "{}/Automation_Scripts/orca_submit.py -p {} -d geoopt -t {} -ppn {} -q {} -sched {} -size {} -o optimizations -path_to_exe {}  --silent"
         substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_ppn"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_size"],c["orca_exe"])
         output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
         output = str(output,'utf-8')
         jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
         monitor_jobs(jobids,user)

         # Remove job submission files
         delete_wildcard('optimization')

         # Remove *.gbw files (take up a LOT of memory)
         delete_gbw()

         # cd back to main running directory
         os.chdir('..')

    # Check if geoopt jobs are incomplete/need resubmission
    os.chdir('Intramolecular_Modes')
    xyzfiles = [f for f in os.listdir('.') if os.path.isfile(f) and f.split('.')[-1] == 'xyz']
    for j in range(0,10):
         resub_flag = 0
         # Check for incompletes/resubmissions
         for i in xyzfiles:
            name = i.split('.')[0]
            if os.path.isdir(name+'_geoopt') is False:
               print("ERROR: {}_geoopt folder wasn't found in {}...".format(name,i))
               quit()
            outfile = '{}_geoopt/{}_geoopt.out'.format(name,name)
            if os.path.isfile(outfile) is False:
               print("WARNING: {} wasn't found in {}...".format(outfile,i)) 
               resub_flag = 1
               break
            flag = find_string(outfile,b'THE OPTIMIZATION HAS CONVERGED')
            if flag is False:
               print("WARNING: {}_geoopt didn't run to completion".format(name))
               resub_flag = 1
               break
         if(resub_flag):
            print("resub of geoopt has been triggered")
            jobids = []
            substring = "{}/Automation_Scripts/orca_submit.py -d geoopt -p {} -t {} -ppn {} -q {} -sched {} -size {} -o optimizations -path_to_exe {}  --silent --resubmit"
            substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_ppn"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
         
            # Wait until jobs complete
            monitor_jobs(jobids,user)

         # If nothing is in need of resubmission then break
         else:
            delete_gbw()
            break
    os.chdir('..')   
    print("COMPLETE: geoopt jobs")
      
######################
# Run the mode scans #
######################

    # Check if bond/angle mode scans are needed (flag is set to 1 if so)
    ba_flag = 0 
    os.chdir('Intramolecular_Modes')
    # Check if the bonds folders are present for all fragments
    ba_folders =  [f for f in os.listdir('.') if os.path.isfile(f) is False and f.split('_')[1] == 'bonds'] 
    if ba_folders == []: ba_flag = 1
    os.chdir('..')

    # Check if dihedral mode scans are needed (flag is set to 1 if so)
    d_flag = 0
    os.chdir('Intramolecular_Modes')
    # Check if the dihedrals folders are present for all fragments
    d_folders =  [f for f in os.listdir('.') if os.path.isfile(f) is False and f.split('_')[-1] == 'dihedrals'] 
    if d_folders == []: d_flag =1
    os.chdir('..')

    print("COMPLETE: check on missing modes")

    # Generate bond/angle mode scan jobs if needed
    jobids = []
    if ba_flag == 1:
         if(c['functional'] == 'wB97X-D3'):
            command_list = ('Intramolecular_Modes -p {} -q {} -gens {} -theory dft -f {} -b {} --no_D3 -modes'.format(c['param_ba_procs'],c['charge'],c['gens'],c['functional'],c['basis'])).split()
         else:
            command_list = ('Intramolecular_Modes -p {} -q {} -gens {} -theory dft -f {} -b {} -modes'.format(c['param_ba_procs'],c['charge'],c['gens'],c['functional'],c['basis'])).split()
         command_list.append("bonds angles")
         paramgen_old.main(command_list)
         substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d bonds_angles -p {} -t {}  -q {} -sched {} -ppn {} -size {} -path_to_exe {} --silent -o bonds_angles"
         substring = substring.format(c["taffi_path"],c["param_ba_procs"],c["param_ba_wt"],c["param_ba_q"],c["param_ba_sched"],c["param_ba_ppn"],c["param_ba_size"],c["orca_exe"])
         output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
         output = str(output,'utf-8')
         jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
         os.chdir('..')
    
    # Generate dihedral mode scan jobs if needed
    if d_flag == 1:
         if(c['functional'] == 'wB97X-D3'):
            command_list = ('Intramolecular_Modes -p {} -q {} -gens {} -theory dft -f {} -b {} --no_D3 -d_step 10.0 --scan -modes'.format(c['param_d_procs'],c['charge'],c['gens'],c['functional'],c['basis'])).split()
         else:
            command_list = ('Intramolecular_Modes -p {} -q {} -gens {} -theory dft -f {} -b {} -d_step 10.0 --scan -modes'.format(c['param_d_procs'],c['charge'],c['gens'],c['functional'],c['basis'])).split()
         command_list.append("dihedrals")
         paramgen_old.main(command_list)
         substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d dihedrals -p {} -t {}  -q {} -sched {} -ppn {} -size {} -path_to_exe {} --silent -o dihedrals"
         substring = substring.format(c["taffi_path"],c["param_d_procs"],c["param_d_wt"],c["param_d_q"],c["param_d_sched"],c["param_d_ppn"],c["param_d_size"],c["orca_exe"])
         output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
         output = str(output,'utf-8')
         jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
         os.chdir('..')

    # Wait until jobs complete
    monitor_jobs(jobids,user)

    print("COMPLETE: QC jobs for missing modes")

    # Remove old test and REDO folders (sometimes the result of an interuption mid-cycle)
    os.chdir('Intramolecular_Modes')
    if(os.path.isdir('test')):
         shutil.rmtree('test')
    for item in os.listdir('.'):
         if os.path.isdir(item):
            if (item.find('REDO') != -1) or (item.find('RESTART') != -1):
               shutil.rmtree(item)
    os.chdir('..')

    # Check for bond/angle resubmissions
    for i in range(0,10):
         # If final_params data is already in place then avoid the resubmissions
         if os.path.isdir('Intramolecular_Modes/final_params') is True and os.path.isfile('Intramolecular_Modes/final_params/Intramolecular_Modes-DFT.db') is True:
            break
         # Check for needed resubmissions using  the extract_intramolecular_params.py script running in --find_restarts mode.
         sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
         sublist.append('{}{}/FF_functions/extract_intramolecular_params_old.py -f Intramolecular_Modes -gens {} --find_restarts -lammps_exe {} -FF \'{}\' -o test -modes \'bonds angles harm_dihedrals\''.format(c["module_string"],c["taffi_path"],c['gens'],c['lammps_exe'],c['ff']))
         sublist= sublist + ('restart  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
         output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf-8').communicate()[0]
         #output = str(output,'utf-8')
         jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
         monitor_jobs(jobids,user)
         flag = find_string("Intramolecular_Modes/test/extract_intramolecular.log",b"Generating the input files for a reoptimization based on the lowest energy")
         shutil.rmtree('Intramolecular_Modes/test')
       
         # Handle flexible dihedral resubmission
         substring = "python3 {}/FF_functions/dihedral_restart.py Intramolecular_Modes --verbose".format(c["taffi_path"])
         output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf-8').communicate()[0]
         #output = str(output,'utf-8')
         if (output.split("\n")[0]) != '':
            D_flag = len(output.split("\n"))
         else:
            D_flag = 0

         if D_flag != 0:
            print("Some flexible dihedral scans did not complete. Resuming...")
            # Submit restarted dihedral scans
            D_ids_sub = []
            os.chdir("Intramolecular_Modes")
            substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d RESTART -p {} -t {}  -q {} -sched {} -ppn {} -size {} -path_to_exe {} --silent -o dihedrals_restart"
            substring = substring.format(c["taffi_path"],c["param_d_procs"],c["param_d_wt"],c["param_d_q"],c["param_d_sched"],c["param_d_ppn"],c["param_d_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf-8').communicate()[0]
            #output = str(output,'utf-8')
            D_ids_sub += [ m.split()[-1] for m in output.split("\n")[:-1]]
            os.chdir('..')

         if flag is True:
            print("Some bond/angle jobs failed to converge. Resubmitting...")
            jobids = []

            # Run geometry optimizations of the unconverged modes
            os.chdir('Intramolecular_Modes')
            substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d REDO -p {} -t {}  -q {} -sched {} -ppn {} -size {} -o ba_geo_resub -path_to_exe {} --silent"
            substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_ppn"],c["param_geoopt_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
            os.chdir('..') # Back to main folder
      
            # Wait until job complete 
            monitor_jobs(jobids,user)

            # Generate and run new mode scans
            if(c['functional'] == 'wB97X-D3'):
               restart_scans.main(("Intra -f {} -b {} -gens {} --no_D3".format(c['basis'],c['functional'],c['gens'])).split())
            else:
               restart_scans.main(("Intra -f {} -b {} -gens {}".format(c['basis'],c['functional'],c['gens'])).split())
            os.chdir('Intramolecular_Modes')
            substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in  -p {} -t {}  -q {} -sched {} -ppn {} -size {} -o ba_scan_resub -path_to_exe {} --silent"
            substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_ppn"],c["param_geoopt_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
            os.chdir('..') # Back to main folder
            monitor_jobs(jobids,user)

            # Remove REDO folder for a fresh start on the next cycle in case something went wrong
            for item in os.listdir('.'):
               if os.path.isdir(item):
                  if(item.find('REDO') != -1):
                     shutil.rmtree(item)

         # Run the dihedral post processing if jobs were submitted.
         if D_flag != 0:  

            # Wait until jobs complete
            monitor_jobs(D_ids_sub,user)

            # Run post processing
            dihedral_restart.main("Intramolecular_Modes --post".split()) 

            # Remove RESTART folder for a fresh start on the next cycle in case something went wrong
            for item in os.listdir('.'):
               if os.path.isdir(item):
                  if(item.find('RESTART') != -1):
                     shutil.rmtree(item)

         # Break out of the restart loop if no resubmissions occured
         if D_flag == 0 and flag == False: break

    print("COMPLETE: QC data for intramolecular mode fit")

    # If there was an aborted attempt to fit then remove the old folder
    if (os.path.isdir('Intramolecular_Modes/final_params')) and (os.path.isfile('Intramolecular_Modes/final_params/Intramolecular_Modes-DFT.db') is False):
         print("Removing old final_params folder...")
         shutil.rmtree('Intramolecular_Modes/final_params')

    # Fit final params
    if os.path.isdir('Intramolecular_Modes/final_params') is False:

         # remove *.gbw files (take up a LOT of memory)
         delete_gbw()

         # Assemble shell submit string
         jobids = []
         sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
         sublist.append('{}{}/FF_functions/extract_intramolecular_params_old.py -f Intramolecular_Modes -o final_params -FF \'{}\' -min_cycles 5 -gens {} -max_cycles 10 -N_sc 3 -lammps_exe {} --mixing_rule wh'.format(c["module_string"],c["taffi_path"],c["ff"],c['gens'],c["lammps_exe"]))
         sublist= sublist + ('final_params  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
         output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
         output = str(output,'utf-8')
         jobids += [ output.split("\n")[-2].split()[-1]]

         # Wait until jobs complete
         monitor_jobs(jobids,user)

    # Check the results of the fit
    if (os.path.isdir('Intramolecular_Modes/final_params') is False) or (os.path.isfile('Intramolecular_Modes/final_params/Intramolecular_Modes-DFT.db') is False):
         print("An error occured during the intramolecular mode fit.")
         quit()

    print("COMPLETE: intramolecular mode fit")
    quit()

def delete_gbw():

    fileList = glob.glob('./*/*.gbw')

    for filePath in fileList:
       try:
          os.remove(filePath)
       except:
          print("Error while deleting file : ", filePath)

def delete_wildcard(pattern):
    fileList = glob.glob('./'+pattern+'.*')
      
    for filePath in fileList:
       try:
          os.remove(filePath)
       except:
          print("Error while deleting file : ", filePath)
         

def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s


def find_string(textfile,finding):
   # finding has to a byte string, ex:b'blabla'
   import mmap

   with open(textfile, 'rb', 0) as file, \
     mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
      if s.find(finding) != -1:
         return True
      else: 
         return False 

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder, "a",buffering = 1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
