#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,getpass,json

# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Parsers')
from monitor_jobs import *
from parse_all import *
import frag_gen_inter,frag_gen,xyz_to_orca,paramgen,merge_FF,extract_intramolecular_params,gen_md_for_sampling,gen_jobs_for_charges 
import plot_vdw_convergence,check_density,rescale_data,gen_jobs_for_vdw
import codecs,json
from file_parsers import xyz_parse



def main(argv):
    
    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-FF', dest='FF_db', default='',
                        help = 'This variable holds the filename(s) of the .db files holding the necessary force-field parameters (default: FF.db)')

    args=parser.parse_args()    

    # start logger
    sys.stdout = Logger(os.getcwd()+'/vdw_loop.log')

    print("{}\n\nPROGRAM CALL: python {}\n".format('-'*150,' '.join([i for i in sys.argv])))
            
    # parse configuration dictionary (c)
    c = parse_configuration(args.config)
    c['user'] = getpass.getuser()
    c['run_dir'] = os.getcwd()
    if(args.FF_db != ''):
      c['ff'] = (args.FF_db)

    # Make config file absolute
    args.config = os.path.abspath(args.config)
    # XXX this should be set by parse_configuration and made consistent across all scripts. 7/20/20 BMS
    if (c['module_string'] == None):
      c['module_string'] = ''
    
    # Hard-coded run parameters.
    # XXX these should be set by parse_configuration as defaults. 7/20/20 BMS
    d_N_thresh="0.01"   # number density threshold for rerunning MD at NVT
    d_N_fix="0.05"      # jobs whose number density falls below d_N_thresh are resubmitted with a number density of d_N_fix
    L2_s="0.1"
    L2_e="0.5"
    max_cycles=20
    max_attempts=10     # Maximum number of md resubmissions to attempt.

    # Determine the starting cycle if jobs are already present    
    if os.path.isdir('vdw'):

        # Find the number of cycle folders
        counter=0
        cyclefolders = []
        for item in os.listdir(c['run_dir']+'/vdw'):
            if os.path.isdir(os.path.join(c['run_dir']+'/vdw', item)):
                if(item.find('cycle') != -1):
                    cyclefolders.append(item)
        counter = len(cyclefolders)

        # Determine the starting cycle
        cycle=1
        rem_flag=0

        # If no cycle folders are present, start on cycle 1 and delete old folders
        if (counter == 0 ):

            cycle = 1

            # remove first md folder if incomplete
            if(os.path.isdir('vdw/md-1')):
                if not (os.path.isfile('vdw/md-1/md-1.end.data')):
                    shutil.rmtree('vdw/md-1')  

            # remove second md folder if incomplete
            if(os.path.isdir('vdw/md-2')):
                if not (os.path.isfile('vdw/md-2/md-2.end.data')):
                    shutil.rmtree('vdw/md-2')  

            # remove third md folder if incomplete
            if(os.path.isdir('vdw/md-3')):
                if not (os.path.isfile('vdw/md-3/md-3.end.data')):
                    shutil.rmtree('vdw/md-3')  

            # remove configs folder if present
            if(os.path.isdir('vdw/configs')):
                shutil.rmtree('vdw/configs')

        # Else, start after the first complete cycle
        else:
            
            for i in range(1,max_cycles+1):

                # Check the status of existing cycles
                # If incomplete, restart here and delete the rest
                md = i+2
                if (os.path.isfile('vdw/cycle-'+str(i)+'/DFT-AA.db')) is False or (os.path.isfile('vdw/md-'+str(md)+'/md-'+str(md)+'.end.data')) is False or (os.path.isdir('vdw/configs/'+str(md)+'-0')) is False: 

                    # Find the number of cycle folders
                    mdfolders = []
                    for item in os.listdir(c['run_dir']+'/vdw'):
                        if os.path.isdir(os.path.join(c['run_dir']+'/vdw', item)):
                            if(item.find('md') != -1):
                                mdfolders.append(item)

                    # loop over md folders and remove any that are greater than or equal to the current md
                    for j in mdfolders:
                        num = j.split('-')[-1]
                        if(int(num) >= md):
                            print("removing md-{}... (rm -r {})".format(num,j))
                            shutil.rmtree(c['run_dir']+'/vdw/'+j)

                            # Remove MD configs for this cycle if present
                            configfolders = []
                            for item in os.listdir(c['run_dir']+'/vdw/configs'):
                                if os.path.isdir(os.path.join(c['run_dir']+'/vdw/configs', item)):
                                    if(item.find(num+'-') != -1):
                                        configfolders.append(item)

                            for k in configfolders:
                                shutil.rmtree(c['run_dir']+'/vdw/configs/'+k)

                    # loop over cycle folders and remove any that are greater than or equal to the current cycle
                    for j in cyclefolders:
                       num = j.split('-')[-1]
                       if (int(num) >= cycle):
                          print("removing cycle-{}... (rm -r {})".format(num,j))
                          shutil.rmtree(c['run_dir']+'/vdw/'+j)

                    # Break out of checking folders
                    break

                else:
                    cycle = i 

                # If complete up to this point, increment starting cycle
                cycle += 1


    # If no job data is present then start with 1
    else:
        cycle=1

    # Outer loop controls cycles
    inchi = os.getcwd().split('/')[-1] # XXX Kluge. We need better file/directory managment for this script. Should probably accept an inchi key/list of keys as argument.
    fix_d=0
    while (cycle <= max_cycles): 

        # Step 1 - Run LAMMPS jobs to sample configurations
        # while loop handles resubmission conditions if the number density of the md simulation drops below d_N_thresh
        cond=0
        attempts=1
        while(cond == 0):

            # Check if the maximum number of resubmission attempts has been reached
            if(attempts > max_attempts):
                print( "ERROR: md failed. Maximum number of resubmission attempts reached.")
                quit()

            jobids = []
            if(cycle == 1):

                # Collate the preliminary intramolecular modes 
                if (os.path.isfile('vdw/intermediate_params.db')) is False:
                    print("ERROR: {}/vdw/intermediate_params.db doesn't exist. Should have been created by driver.py.".format(os.getcwd()))
                    quit()

                # STEP 1: RUN MD (cycle == 1)
                elements,geo = xyz_parse('{}.xyz'.format(inchi))
                N_atoms = len(elements)
                N_mol = int(1000/N_atoms) + 1

                # Make runfolder
                # XXX shouldn't be needed anymore because the driver generates this folder. 7/20/20
                if (os.path.isdir('vdw')) is False:
                    os.makedirs('vdw')

                # For resubmissions use a lower density
                if (attempts >= 1):
                    ts=1.0
                    d_N=0.08
                    t_A=1E5
                else:
                    ts=1.0
                    d_N=0.09
                    t_A=1E5

                # Make runs (if/then are in place in case of resubmission)
                if (os.path.isdir('vdw/md-1')) is False:
                    sublist = ("{}.xyz -N {} -o vdw/md-1 -d_N {} -sigma_scale 0.6 -charge_scale 0.0 -T_A 10 -T 298 -t_A {} -t 1E5 --UFF -q 0 -ts {} -gens {} --molecule -FF".format(inchi,N_mol,d_N,t_A,ts,c['gens'])).split()
                    sublist.append("vdw/intermediate_params.db {}".format(c['ff']))
                    gen_md_for_sampling.main(sublist)
                if (os.path.isdir('vdw/md-2')) is False:
                    sublist = ("{}.xyz -N {} -o vdw/md-2 -d_N {} -sigma_scale 0.8 -charge_scale 0.0 -T_A 10 -T 298 -t_A {} -t 1E5 --UFF -q 0 -ts {} -gens {} --molecule -FF".format(inchi,N_mol,d_N,t_A,ts,c['gens'])).split()
                    sublist.append("vdw/intermediate_params.db {}".format(c['ff']))
                    gen_md_for_sampling.main(sublist)
                if (os.path.isdir('vdw/md-3')) is False:
                    sublist = ("{}.xyz -N {} -o vdw/md-3 -d_N {} -sigma_scale 1.0 -charge_scale 0.0 -T_A 10 -T 298 -t_A {} -t 1E5 --UFF -q 0 -ts {} -gens {} --molecule -FF".format(inchi,N_mol,d_N,t_A,ts,c['gens'])).split()
                    sublist.append("vdw/intermediate_params.db {}".format(c['ff']))
                    gen_md_for_sampling.main(sublist)

                # Submit the jobs
                os.chdir('vdw')
                mdfolders = []
                for item in os.listdir(c['run_dir']+'/vdw'):
                    if os.path.isdir(os.path.join(c['run_dir']+'/vdw', item)):
                        if(item.find('md') != -1):
                            mdfolders.append(item)
                for j in mdfolders:
                    os.chdir(j)
                    name = j
                    print("Submitting {}...".format(name))
                    substring = "python3 {}/Automation_Scripts/lammps_submit.py -f {}.in.init -p {} -t {} -q {} -o {} -sched {} -ppn {} -size {} -npp {} -path_to_exe {} --silent"
                    substring = substring.format(c["taffi_path"],name,c['vdw_md_procs'],c['vdw_md_wt'],c['vdw_md_q'],name,c['vdw_md_sched'],c['vdw_md_ppn'],c['vdw_md_size'],c['vdw_md_npp'],c['lammps_exe'])
                    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                    output = str(output,'utf-8')
                    jobids += [m.split()[-1] for m in output.split("\n")[:-1]]
                    os.chdir('..')
                os.chdir('..')

            # STEP 1: RUN MD (cycle>1)
            else:

                md = cycle +2
                elements,geo = xyz_parse('{}.xyz'.format(inchi))
                N_atoms = len(elements)
                N_mol = int(1000/N_atoms) + 1
                jobids = []
                # If the constant number density flag has been triggered, run NVT
                if (fix_d == 1):
                   gen_md_for_sampling.main(("vdw -N {} -sigma_scale 0.6 -o vdw -T 298 -t_A 1E5 -t 1E5 -T_A 100 -d_N {} -gens {} --molecule".format(N_mol,d_N_fix,c['gens'])).split())

                # Else, run NPT
                else:
                   gen_md_for_sampling.main(("vdw -N {} -sigma_scale 0.6 -o vdw -T 298 -t_A 1E5 -t 1E5 -T_A 100 -gens {} --molecule".format(N_mol,c['gens'])).split())

                # Submit the MD simulation to the cluster (ids_sub is used to keep track of jobs awaiting completion)
                # If/else loops bracketing the directory switches have been added throughout the script for graceful failure.  
                os.chdir(c['run_dir'])
                if (os.path.isdir('vdw/md-'+str(md))):
                   os.chdir('vdw/md-'+str(md))
                else:
                   print("ERROR in vdw_loop: failed to create vdw/md-{}. Exiting...".format(md))
                   quit()
                # Check if input file was actually generated
                if os.path.isfile('md-'+str(md)+'.in.init') is False:
                   print("ERROR in vdw_loop: failed to create vdw/md-{}.in.init.  Exiting...".format(md))    
                   quit()

                name='md-'+str(md)
                substring = "python3 {}/Automation_Scripts/lammps_submit.py -f {}.in.init -p {} -t {} -q {} -o {} -sched {} -ppn {} -size {} -npp {} -path_to_exe {} --silent"
                substring = substring.format(c["taffi_path"],name,c['vdw_md_procs'],c['vdw_md_wt'],c['vdw_md_q'],name,c['vdw_md_sched'],c['vdw_md_ppn'],c['vdw_md_size'],c['vdw_md_npp'],c['lammps_exe'])
                output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                output = str(output,'utf-8')
                jobids += [m.split()[-1] for m in output.split("\n")[:-1]]
                os.chdir(c['run_dir'])

            # Wait until jobs complete
            monitor_jobs(jobids,c['user'])

            # First cycle, check that the job completed
            if (cycle == 1):
                cond=1
                folders=["md-1","md-2", "md-3"]
                for j in folders:
                    if os.path.isfile('vdw/'+j+'/'+j+'.end.data') is False:
                        print("WARNING: vdw/{} didn't complete. Attempting resubmission...".format(j))
                        cond=0
                        shutil.rmtree('vdw/'+j)

            # Other cycles, check that the system hasn't evaporated and that the job completed
            else:
                # Check the number density in the latest MD simulation
                # Since in the first cycle the number density is held constant, the break condition is guarranteed to be met
                md = cycle +2
                d_N=check_density.main(["vdw/md-{}/md-{}.end.data".format(md,md)])
                if (d_N < 0.01):
                    # Print diagnostic
                    print("Switching over to NVT simulations at a constant number density of {} atoms/A^3...".format(d_N_fix)) 

                    # Toggle NVT flag, remove the most recent run, backup the penultimate *end.data file, 
                    # and rescale the box in the penultimate *end.data file for use in the resubmitted job
                    fix_d=1
                    shutil.rmtree('vdw/md-'+str(md))
                    md_old = cycle +1
                    source = "vdw/md-{}/md-{}.end.data".format(md_old,md_old)
                    destination = "vdw/md-{}/md-{}.end.data.unscaled".format(md_old,md_old)
                    dest = shutil.copyfile(source, destination) 
                    rescale_data.main(("vdw/md-{}/md-{}.end.data -d_N {}".format(md_old,md_old,d_N_fix)).split())

                elif os.path.isfile('vdw/md-'+str(md)+'/thermo.avg') is False or os.path.isfile('vdw/md-'+str(md)+'/equil.lammpstrj') is False:
                    print( "WARNING: md-{} didn't complete. Attempting resubmission...".format(md))
                    cond=0
                    shutil.rmtree('vdw/md-'+str(md))
                else:
                    cond=1

            # Increment attempts
            attempts += 1
      	
        # STEP 2: RUN ORCA JOBS
        jobids=[]

        # 20 cycles
        sublist = "vdw -r_min_scale 1.5 -every 10 -N 100 -p {} {} -QC_types dft -f {} --remove_cross -FF".format(c['vdw_qc_procs'],c['nod3'],c['functional']).split()
        sublist.append(c['ff'])
        gen_jobs_for_vdw.main(sublist)

        """
        if(cycle == 1):
            sublist = "vdw -r_min_scale 1.5 -every 10 -N 100 -p {} {} -QC_types dft -f {} --remove_cross -FF".format(c['vdw_qc_procs'],c['nod3'],c['functional']).split()
            sublist.append(c['ff'])
            gen_jobs_for_vdw.main(sublist)
        else:
            sublist = "vdw -r_min_scale 1.5 -every 10 -N 1700 -p {} {} -QC_types dft -f {} --remove_cross -FF".format(c['vdw_qc_procs'],c['nod3'],c['functional']).split()
            sublist.append(c['ff'])
            gen_jobs_for_vdw.main(sublist)
        """

        # Submit QC jobs
        os.chdir('vdw')
        substring = "python3 {}/Automation_Scripts/orca_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -path_to_exe {}  --silent"
        substring = substring.format(c["taffi_path"],c["vdw_qc_procs"],c["vdw_qc_wt"],c["vdw_qc_ppn"],c["vdw_qc_q"],c["vdw_qc_sched"],c["vdw_qc_size"],c["orca_exe"])
        output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
        output = str(output,'utf-8')
        jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
        os.chdir(c['run_dir'])

        # Wait until jobs complete
        monitor_jobs(jobids,c['user'])

        # remove *.gbw files (take up a LOT of memory)
        os.chdir('vdw')
        delete_gbw()

        # Remove submission file
        delete_wildcard('orca_mass')
        os.chdir(c['run_dir'])

        # STEP 3: FIT VDW PARAMETERS
        jobids = []
        sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
        sublist.append('{}/FF_functions/extract_vdw.py -f vdw -FF_DFT \'intermediate_params.db {}\' -o cycle-{} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} -gens {}'.format(c['taffi_path'],c['ff'],cycle,L2_s,L2_e,c['functional'],c['gens']))
        sublist= sublist + ('vdw_parse  -p 1 -t {} -q {} -ppn {}'.format(c["param_fit_wt"],c["param_fit_q"],c["param_fit_ppn"])).split()
        output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
        output = str(output,'utf-8')
        jobids += [ output.split("\n")[-2].split()[-1]]

        # Wait until jobs complete
        monitor_jobs(jobids,c['user'])

        # Update convergence condition
        cycle += 1
        if ( cycle <= max_cycles) :
            with open("vdw/status.txt",'w') as f:
                f.write("incomplete")
        else:
            with open("vdw/status.txt",'w') as f:
                f.write("complete")

    # Plot parameter convergence
    # XXX should put all of these inside of inchi/vdw not in inchi/. 7/21/20
    if os.path.isfile('vdw_cycle.eps.pdf') is False or os.path.isfile('vdw_cycle.sigmas.pdf') is False:
        plot_vdw_convergence.main(['vdw'])

    # Merge modes and charges into intermediate param dictionary
    # XXX should put all of these inside of inchi/vdw not in inchi/. 7/21/20
    merge_FF.main(("final_vdw.db vdw/cycle-{}/DFT-AA.db -R".format(max_cycles)).split())
    merge_FF.main(("final_vdw.db vdw/cycle-{}/DFT-UA.db -R".format(max_cycles)).split())
    
    print("---------VDW LOOP ENDS-----------\n")
    quit()

# XXX Add documentation for all help functions. 7/21/20
def delete_gbw():

    fileList = glob.glob('./*.gbw')

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

def find_string(textfile,finding):
   # finding has to a byte string, ex:b'blabla'
   import mmap

   with open(textfile, 'rb', 0) as file, \
     mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
      if s.find(finding) != -1:
         return True
      else: 
         return False 
      
   

def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

def read_json(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    for key in jload:
      jload[key]['geo'] = np.array(jload[key]['geo'])
    
    return jload

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder, "a",buffering = 1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
        self.flush()
    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
