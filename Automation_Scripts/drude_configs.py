#!/bin/env python                                                                                                                                                              

import sys,os,argparse,subprocess,shutil,time,glob,getpass,json,fnmatch,codecs,random

# For plotting (Agg called needed for cluster image generation)
from pylab import *
import numpy as np
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
from monitor_jobs import *
from file_parsers import xyz_parse
from parse_all import parse_configuration
import generate_connolly


def main(argv):
    

    global user,run_dir,c
    user = getpass.getuser()
    current_dir = os.getcwd()

    

    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-folder', dest='folder', default='',
                        help = 'inchikey folder name')
    parser.add_argument('-o',dest='output',default='polar.log',
                        help = 'log name for stdout')
    parser.add_argument('-tag',dest='tag',default='',
                        help = '_tag name default: None')
    parser.add_argument('-nonpolar_dir',dest='nonpolar_dir',default='',
                        help = 'nonpolar directory, this is the directory you execute normal taffi in')
    parser.add_argument('--write_bash', dest='write_bash', default=False, action='store_const', const=True,
                        help = 'When this flag is supplied, will only write out bash with generation info') 

    args=parser.parse_args()    

    if args.nonpolar_dir == '':
      print('ERROE: nonpolar directory not specified, Exiting....')
      quit()

    if args.write_bash:
       # Read all.json (contains model compound and dependency information for the batch 
       status = read_alljson('./all.json')
       # Find inchi keys for this generation
       mc_keys = [ _ for _ in status["mc"].keys() ]
       # Find inchi keys for this generation
       #mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]
       gens = sorted(status["gens"].keys())
       with open('drude_gens.sh','w') as f:
         f.write("#!/bin/bash\n")
         for g in gens:
            keys = [ _ for _ in status["mc"].keys() if int(g) == status["mc"][_]["gen"][0] ] 
            f.write("gen{}=(".format(g))
            for i in keys:
               f.write('\'{}\' '.format(i))
            f.write(')\n')
         f.write("\nfor i in \"${gen1[@]}\"; do \n")
         f.write("   screen -X -S $i quit \n")
         f.write("   screen -dmS $i \n")
         f.write("   screen -S $i -X stuff \"ana3 ; ~/taffi_beta/Automation_Scripts/drude_configs.py -nonpolar_dir {} -c {} -folder $i \\n\" \n".format(args.nonpolar_dir,args.config))
         f.write("done \n")
       subprocess.call("chmod u+x drude_gens.sh",shell=True)
       quit()
    
    sys.stdout = Logger(current_dir+'/'+args.folder+'/'+args.output)
    
    # parse configuration dictionary (c)
    print("parsing configurations...")
    c = parse_configuration(args.config)
    # Make config file absolute
    args.config = os.path.abspath(args.config)
    if (c['module_string'] == None):
      c['module_string'] = ''


    fit_Drude_custom(current_dir,args.folder,args,args.tag)

    return


# Normal Drude fitting, fit all heavy atoms
def fit_Drude_custom(driver_dir,folder,args,tag):

    global run_dir

    run_dir = '{}/{}'.format(driver_dir,folder)
    os.chdir(run_dir) 
 
    nonpolar = args.nonpolar_dir 
    c_config ='{}/{}/charges/CHELPG_calcs/configs'.format(nonpolar,folder)
    nconf = sorted([int(f) for f in os.listdir(c_config) if os.path.isdir(os.path.join(c_config, f))])[-1]+1
    if nconf > 100: nconf = 100
   
    # placing perturbation charges 
    jobids = []
    print("Start placing perturbation charges...")

    # try five times incase something went wrong
    for retry in range(0,5):
       for count_i,i in enumerate(list(range(0,nconf))):
         fitfolder = 'fitcharge_{}'.format(i)
         if os.path.isdir(run_dir+'/'+fitfolder) is False:
            os.makedirs(run_dir+'/'+fitfolder)
         xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
         element,geo = xyz_parse(xyz)
         jobids += place_perturbation_submit(fitfolder,xyz,dplus=retry*100)
       if jobids == []: break
       monitor_jobs(jobids,user)
    clean_pcsubmit(nconf)

    # do QM calculation of polarizability and dipole for the origianl xyz for later FOM to compare to
    print("Starting QM polar calculation...")
    for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
      element,geo = xyz_parse(xyz)
      polar_orca(fitfolder,geo,element)
  
    os.chdir(run_dir) 
    substring = "{}/Automation_Scripts/orca_submit.py -f polar.in -p {}  -t {} -ppn {} -q {} -sched {} -size {} -o optimizations -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],c['charges_qc_procs'],c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_size"],c["orca_exe"])
    orcaids = execute_orca_submit(substring)
    
    print("Starting CHELPG potential calculation...")
    for retry in range(0,5):
       jobids = []
       for q in [-0.5,0.5]:
          for count_i,i in enumerate(list(range(0,nconf))):
            fitfolder = 'fitcharge_{}'.format(i)
            xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
            element,geo = xyz_parse(xyz)
            total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
            jobids += calc_CHELPG(fitfolder,xyz,total_pc,pc_layer1,pc_layer2,pc_layer3,geo,element,q,'_{}'.format(q))
       monitor_jobs(jobids,user)
       clean_chelpg(nconf)

    print("Starting reevaluation of potential at Connolly grid...")
    for retry in range(0,5):
       jobids = []
       for q in [-0.5,0.5]:
          for count_i,i in enumerate(list(range(0,nconf))):
            fitfolder = 'fitcharge_{}'.format(i)
            xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
            total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
            connolly_grid = get_connolly_grid(fitfolder,xyz)
            jobids += calc_connolly_pot(fitfolder,total_pc,connolly_grid,c,'_{}'.format(q))
       monitor_jobs(jobids,user)
       clean_connolly(nconf)
   
    print("Starting rewrite Connolly potential format...")
    for q in [-0.5,0.5]:
       for count_i,i in enumerate(list(range(0,nconf))):
         fitfolder = 'fitcharge_{}'.format(i)
         xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
         total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
         connolly_grid = get_connolly_grid(fitfolder,xyz)
         rewrite_pot_wrap(fitfolder,total_pc,connolly_grid,'_{}'.format(q))
   
    print("Starting writing charmm data...") 
    # Combine data to data_ dir
    for q in [-0.5,0.5]:
       for count_i,i in enumerate(list(range(0,nconf))):
         fitfolder = 'fitcharge_{}'.format(i)
         xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
         element,geo = xyz_parse(xyz)
         total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
         write_charmm_data(geo,fitfolder,fitfolder,list(range(0,total_pc)),'_{}'.format(q))
       write_datadir('_{}'.format(q))

    # tar configs files to save space
    for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      tar_file(fitfolder)

    monitor_jobs(orcaids,user)
    clean_polar(nconf)
    

    # XXX: need to check if this cleans the files real good
    delete_wildcard('','optimization')
    print("Drude QM calculation succeed!")
    quit()
   
    return

def place_perturbation_submit(calc_folder,xyz,dplus=0):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    jobids = []
    if os.path.isdir('pointcharge_layer1') is False:
       d = 700+int(dplus)
       commands = 'cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 13 -c 1 -extra 0.1 -o layer1 {} -ofolder pointcharge_layer1 \n'.format(run_dir,calc_folder,c['taffi_path'],d,xyz)
       write_shell(c["param_fit_q"],1,'pc_layer1',c["param_fit_wt"],'pc_layer1',commands)
       jobid = submit_job_nomonitor('pc_layer1.sh')
       jobids += jobid
    else:
       if len(os.listdir('pointcharge_layer1')) == 0:
         shutil.rmtree('pointcharge_layer1')
    if os.path.isdir('pointcharge_layer2') is False:
       d = 500 +int(dplus)
       commands = 'cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 22 -c 1 -extra 0.02 -o layer2 {} -ofolder pointcharge_layer2 \n'.format(run_dir,calc_folder,c['taffi_path'],d,xyz)
       write_shell(c["param_fit_q"],1,'pc_layer2',c["param_fit_wt"],'pc_layer2',commands)
       jobid = submit_job_nomonitor('pc_layer2.sh')
       jobids += jobid
    else:
       if len(os.listdir('pointcharge_layer2')) == 0:
         shutil.rmtree('pointcharge_layer2')
    if os.path.isdir('pointcharge_layer3') is False:
       d = 500 +int(dplus)
       commands = 'cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 40 -c 1 -extra 0.003 -o layer3 {} -ofolder pointcharge_layer3 \n'.format(run_dir,calc_folder,c['taffi_path'],d,xyz)) 
       write_shell(c["param_fit_q"],1,'pc_layer3',c["param_fit_wt"],'pc_layer3',commands)
       jobid = submit_job_nomonitor('pc_layer3.sh')
       jobids += jobid
    else:
       if len(os.listdir('pointcharge_layer3')) == 0:
         shutil.rmtree('pointcharge_layer3')

    print("Submitting place_charge jobs for {}...".format(calc_folder))
   
    os.chdir(ori_dir)

    return jobids

def execute_orca_submit(substring):
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    orcaids = [ m.split()[-1] for m in output.split("\n")[:-1]]

    return orcaids

def submit_job_nomonitor(shname):
    command  = 'sbatch {}'.format(shname).split()
    output = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf-8').communicate()[0]
    jobid = [ output.split("\n")[-2].split()[-1]]
   
    return jobid 

def write_shell(queue,cpu,shname,walltime,jobname,commands):

   with open('{}.sh'.format(shname),'w') as f:
      f.write('#!/bin/bash\n')
      f.write('#\n')
      f.write('#SBATCH --job-name {}\n'.format(jobname))
      f.write('#SBATCH -o sh.out\n')
      f.write('#SBATCH -e sh.err\n')
      f.write('#SBATCH -A {}\n'.format(queue))
      f.write('#SBATCH -N 1\n')
      f.write('#SBATCH -n {}\n'.format(cpu)) 
      f.write('#SBATCH -t {}:00:00\n'.format(walltime)
      """
      if queue == 'bsavoie': 
         f.write('#SBATCH -t 14-00:00:00\n')
         #f.write('#SBATCH -t 1-00:00:00\n')
      elif queue == 'standby':
         f.write('#SBATCH -t 4:00:00\n')
      elif queue == 'debug':
         f.write('#SBATCH -t 00:30:00\n')
      else:
         print('ERROR: unrecognized queue: {}'.format(queue))
         quit()
      """
      f.write('\n')
      f.write(commands)
      f.write('\n')

   return

def polar_orca(calc_dir,Geometry,Elements):
   os.chdir('{}/{}'.format(run_dir,calc_dir))
   random_factor = 0.0
   with open("polar.in",'w') as f:
      f.write("# calculation of polarizability\n")
      f.write("! wB97X-D3 def2-TZVP TIGHTSCF CHELPG PMODEL PAL4\n")

      f.write("\n%base \"polar\"\n\n")

      f.write("%elprop Polar 1\n")
      f.write("        dipole true\n")
      f.write("        end\n")
      f.write("* xyz {} {}\n".format(0,1))
      for count_i,i in enumerate(Geometry):
         f.write("  {:3s}".format(Elements[count_i]))
         for j in i:
             f.write(" {:< 16.8f}".format((random.random()*2.0-1.0)*random_factor+j))
         f.write("\n")
      f.write("*\n")
   return

def get_pc(calc_folder):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    if os.path.isdir('pointcharge_layer1'):
       if len(os.listdir('pointcharge_layer1')) == 0:
         shutil.rmtree('pointcharge_layer1')
         print("ERROR: charge for layer 1 not placed for {}, try increasing density".format(calc_folder))
       else:
         pc_layer1 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer1') if os.path.isfile(os.path.join('pointcharge_layer1', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 1 not placed for {}".format(calc_folder))
    if os.path.isdir('pointcharge_layer2'):
       if len(os.listdir('pointcharge_layer2')) == 0:
         shutil.rmtree('pointcharge_layer2')
         print("ERROR: charge for layer 2 not placed for {}, try increasing density".format(calc_folder))
       else:
         pc_layer2 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer2') if os.path.isfile(os.path.join('pointcharge_layer2', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 2 not placed for {}".format(calc_folder))
    if os.path.isdir('pointcharge_layer3'):
       if len(os.listdir('pointcharge_layer3')) == 0:
         shutil.rmtree('pointcharge_layer3')
         print("ERROR: charge for layer 3 not placed for {}, try increasing density".format(calc_folder))
       else:
         pc_layer3 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer3') if os.path.isfile(os.path.join('pointcharge_layer3', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 3 not placed for {}".format(calc_folder))
   
    total_pc = pc_layer1 + pc_layer2 + pc_layer3

    print("Successfully place {} perturbations charges for {}, layer1: {}, layer2: {}, layer3: {}".format(total_pc,calc_folder,pc_layer1,pc_layer2,pc_layer3))
   
    os.chdir(ori_dir)

    return total_pc,pc_layer1,pc_layer2,pc_layer3

# Calculating CHELPG potential with pertubation charge
def calc_CHELPG(calc_folder,xyz,total_pc,pc_layer1,pc_layer2,pc_layer3,geo,element,Qmag,tag):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)
    config_folder = 'configs{}'.format(tag)

    if os.path.isdir(config_folder) is False:
         os.makedirs(config_folder)
    os.chdir("{}/{}/{}".format(run_dir,calc_folder,config_folder))
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else:
            scfpfile='{}/{}_charges.scfp'.format(i,i)
            if  os.path.isfile(scfpfile) is True: 
               continue
         os.chdir(str(i))
         #shutil.copyfile(xyz,'{}.xyz'.format(i))
         if i<pc_layer1:
            shutil.copyfile(run_dir+"/{}/pointcharge_layer1/pointcharges.{}.pc".format(calc_folder,i),'pointcharges.pc')
         elif i< pc_layer1+pc_layer2:
            shutil.copyfile(run_dir+"/{}/pointcharge_layer2/pointcharges.{}.pc".format(calc_folder,i-pc_layer1),'pointcharges.pc')
         else:
            shutil.copyfile(run_dir+"/{}/pointcharge_layer3/pointcharges.{}.pc".format(calc_folder,i-(pc_layer1+pc_layer2)),'pointcharges.pc')
   
         # change charge magnitude
         with open('pointcharges.pc', 'r') as f:
            # read a list of lines into data
            lines = f.readlines()
         fields = lines[1].split()
         fields[0] = '{}'.format(str(Qmag)) 
         lines[1] = " ".join(fields)
         with open('pointcharges.pc', 'w') as f:
           f.writelines( lines )
           
         # Write orca input file   
         name = "{}_charges".format(i)
         Write_charge_in(name,geo,element)
         scfpfile='{}_charges.scfp'.format(i)
         os.chdir("{}/{}/{}".format(run_dir,calc_folder,config_folder))
         

    # Create directory for unperturbed charge potential
    if os.path.isdir("unperturbed") is False: os.makedirs("unperturbed") 
    scfpfile='unperturbed/unperturbed_charges.scfp'
    if  os.path.isfile(scfpfile) is False: 
       os.chdir("unperturbed")
       #shutil.copyfile(xyz,'unperturbed.xyz')
       Write_charge_in('unperturbed_charges',geo,element)
       os.chdir("{}/{}/{}".format(run_dir,calc_folder,config_folder))

    # Submit orca job
    substring = "python3 {}/Automation_Scripts/orca_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],c["charges_qc_procs"],c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_size"],c["orca_exe"])
    jobids = execute_orca_submit(substring)

    print("COMPLETE: CHELPG potential submission for {}".format(calc_folder))
    os.chdir(ori_dir)

    return jobids

def Write_charge_in(name,geo,element):
   with open('charge.in','w') as f:
      # XXX: hard coding, need to sort this out once rewrite orca_submit.py 07/26/2022
      procs = 8
      D3_option = "D3BJ "
      if(c['functional'] == 'wB97X-D3'):
         D3_option=""
      f.write('#Run DFT-D3 single point of the dimer\n')
      f.write('! {} {} TIGHTSCF {} CHELPG PAL{} KeepDens\n\n'.format(c['functional'],c['basis'],D3_option,procs))
      f.write('%base "{}"\n\n'.format(name))
      if name != 'unperturbed_charges':
         f.write('% pointcharges "pointcharges.pc"\n\n')
      # Increase default max scf iterations
      f.write('\n\n%scf\nMaxIter 1000\nend\n')

      # Add print statements for different charge definitions
      f.write('%output\n')
      f.write('  {}\n  {}\n  {}\n  {}\n  {}\nend\n'.format('Print[ P_Mayer ] 1','Print[ P_NatPop ] 1','Print[ P_Hirshfeld ] 1', 'Print[ P_Mulliken ] 1','Print [ P_Loewdin ] 1'))
      # Write Coordinate header
      f.write('\n* xyz 0 1\n')
      # Write Geometry 
      for count_j,j in enumerate(geo):
          f.write('  {:<20s}   {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(element[count_j],j[0],j[1],j[2]))
      f.write('*\n')
   return

def get_connolly_grid(calc_folder,xyz):

    ori_dir = os.getcwd()
    # Genererating connolly surface grid
    os.chdir(run_dir+'/'+calc_folder)
    if os.path.isfile('map.xyz') is False:
         sublist = ('{} -d 1.41'.format(xyz)).split()
         connolly_grid = generate_connolly.main(sublist)
    else:
         with open('map.xyz','r') as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0:
                  connolly_grid = int(fields[0]) 
                  break
    os.chdir(ori_dir)

    return connolly_grid

# Revaluating potential at connolly grid
def calc_connolly_pot(calc_folder,total_pc,connolly_grid,c,tag):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag)

    if os.path.isdir(config_folder) is False:
         os.makedirs(config_folder)
    os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    for i in range(0,total_pc):
         # make folders and check for completion
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile(str(i)+'/tmp.output') is True: 
               # write an out file so that bundle_submit will leave it alone
               with open(str(i)+'/charge.out','w') as f:
                  f.write("complete")
               continue
         os.chdir(str(i))
         
         # copy necessary file from configs folder
         #prepare(run_dir,calc_folder,i,tag)
         if os.path.isfile('tmp.output') is False:
            with open('charge.in','w') as f:
               f.write("#!/bin/bash\n")
               f.write("{} {}".format((c['module_string'].split()[-2]).replace('\\n','\n'),(c['module_string'].split()[-1]).replace('\\n','\n')))
               f.write("/depot/bsavoie/apps/orca_4_1_2/orca_vpot {0:}/{1:}/{2:}/{3:}/{3:}_charges.gbw {0:}/{1:}/{2:}/{3:}/{3:}_charges.scfp {0:}/{1:}/map.xyz tmp.output\n".format(run_dir,calc_folder,config_ori_folder,i))
         os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
   
    # Submit orca job
    substring = "python3 {}/Automation_Scripts/bundle_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],1,c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_ppn"],c["orca_exe"])
    jobids = execute_orca_submit(substring)
       
    if os.path.isdir('unperturbed') is False: os.makedirs('unperturbed')
    os.chdir('unperturbed')
    i='unperturbed'
    #prepare(run_dir,calc_folder,'unperturbed',tag)
    if os.path.isfile('tmp.output') is False:
       with open('charge.sh','w') as f:
            f.write("#!/bin/bash\n")
            f.write("{} {}".format((c['module_string'].split()[-2]).replace('\\n','\n'),(c['module_string'].split()[-1]).replace('\\n','\n')))
            f.write("/depot/bsavoie/apps/orca_4_1_2/orca_vpot {0:}/{1:}/{2:}/{3:}/{3:}_charges.gbw {0:}/{1:}/{2:}/{3:}/{3:}_charges.scfp {0:}/{1:}/map.xyz tmp.output\n".format(run_dir,calc_folder,config_ori_folder,i))
       shell_exe('sh charge.sh')
    oldpot = '{}/{}/{}/{}/{}_charges.vpot'.format(run_dir,calc_folder,config_ori_folder,i,i)
    newpot = '{}_charges.vpot'.format(i)
    rewrite_potential(oldpot,newpot,connolly_grid)
    print("COMPLETE: revaluation of potential at connolly grid for {}".format(calc_folder))
    os.chdir(run_dir+'/'+calc_folder)
    os.chdir(ori_dir)

    return jobids

def shell_exe(command):
    command = command.split()
    process = subprocess.Popen(command)
    status = process.poll()
    while status == None:
       status = process.poll()
       time.sleep(3)
    return

# Transform the output potential from orca_vpot to chelpg output potential
# needed the chelpg output for mol's xyz 
#    tmp.output(output from orca_vpot) format: 
#    gridx gridy gridz pot
#    gridx gridy gridz pot
#    ....
#    *_charges.vpot(original from chelpg) format:
#    natoms total_grid_pt
#         xyz of mol
#         ...
#         ...
#    pot  gridx gridy gridz
#    pot  gridx gridy gridz
#         ....
def rewrite_potential(oldpot,newpot,connolly_grid):
   with open(oldpot,'r') as f:
      with open(newpot,'w') as g:
        atom_count = 0
        pot_count = 0
        for lc,lines in enumerate(f):
            fields = lines.split()

            # Initialize lists/arrays
            if lc == 0: 
                n_atom    = int(fields[0])
                cosmo_grid    = int(fields[1])
                g.write("{} {}\n".format(n_atom,connolly_grid))
                flag = 1
                continue

            # Parse the molecule
            if flag == 1:
                if len(fields) == 0: continue
                g.write(lines)
                atom_count += 1
                if atom_count == n_atom:
                    flag = 2
                    break

   with open('tmp.output','r') as f:
        with open(newpot,'a') as g:
           for lc,lines in enumerate(f):
               fields = lines.split()
            # Initialize lists/arrays
               if lc == 0: 
                   if connolly_grid != int(fields[0]):
                        print("ERROR: QM output connolly grid number: {} doesn't match map.xyz:{} ".format(fields[0],connolly_grid))
                        quit()
                   flag = 1
                   continue
               if flag == 1 :
                  g.write("{:< 10.7e}     {:< 10.7e}  {:< 10.7e}  {:< 10.7e}\n".format(float(fields[3]),float(fields[0]),float(fields[1]),float(fields[2])))
   return

def rewrite_pot_wrap(calc_folder,total_pc,connolly_grid,tag):
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag)
    os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile(str(i)+'_charges.vpot') is True: continue
         os.chdir(str(i))

         oldpot = '{}/{}/{}/{}/{}_charges.vpot'.format(run_dir,calc_folder,config_ori_folder,i,i)
         newpot = '{}_charges.vpot'.format(i)
         rewrite_potential(oldpot,newpot,connolly_grid)
         os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    return


# calc_folder: folder that has config_connolly
# fit_folder: folder where the fitting is taking place, the file are written in fit_folder/data/ folder
# start_pc, end_pc : control the number of perturbation charges
def write_charmm_data(geo,fitfolder,calc_folder,samples,tag):
  
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag) 
    skip_grid = 0
    
    data_dir = '{}/{}/data{}'.format(run_dir,fitfolder,tag)
    # Rewrite everything to CHARMM compatible format
    if os.path.isdir(data_dir) is False:
       os.makedirs(data_dir)

    Ang2Bohr = 1.889725989 
    # ESP unit in CHARMM: e/A
    # ESP unit in ORCA: e/Bohr
    # length unit: CHARMM: Angstrom, ORCA: Bohr
    fitpot = "{}/fitcharge.pot".format(data_dir)
    with open(fitpot,'w') as g:
       for i in samples:
          filename = "{}/{}/{}/{}/{}_charges.vpot".format(run_dir,calc_folder,config_folder,i,i) 
          with open(filename,'r') as f:
              atom_count = 0
              pot_count = 0
              for lc,lines in enumerate(f):
                  fields = lines.split()

                  # Initialize lists/arrays
                  if lc == 0: 
                      n_atom    = int(fields[0])
                      n_grid    = int(fields[1])
                      flag = 1
                      continue

                  # Parse the grid and potential
                  #if flag == 2:
                  if flag == 2 and lc > (n_atom+skip_grid):
                      g.write("{:< 10.7e}\n".format(float(fields[0])*Ang2Bohr))
                      pot_count += 1
                      if pot_count == n_grid:
                          break

                  # Parse the molecule
                  if flag == 1:
                      if len(fields) == 0: continue
                      atom_count += 1
                      if atom_count == n_atom:
                          flag = 2
                          continue
      
    fitpot0 = "{}/fitcharge.pot0".format(data_dir)
    with open(fitpot0,'w') as g:
       i = 'unperturbed'
       filename = "{}/{}/{}/{}/{}_charges.vpot".format(run_dir,calc_folder,config_folder,i,i) 
       with open(filename,'r') as f:
           atom_count = 0
           pot_count = 0
           for lc,lines in enumerate(f):
               fields = lines.split()

               # Initialize lists/arrays
               if lc == 0: 
                   n_atom    = int(fields[0])
                   n_grid    = int(fields[1])
                   flag = 1
                   continue

               # Parse the grid and potential
               #if flag == 2:
               if flag == 2 and lc > (n_atom+skip_grid):
                   g.write("{:< 10.7e} {:< 10.7e} {:< 10.7e} {:< 10.7e}\n".format(float(fields[1])/Ang2Bohr,float(fields[2])/Ang2Bohr,float(fields[3])/Ang2Bohr,float(fields[0])*Ang2Bohr))
                   pot_count += 1
                   if pot_count == n_grid:
                       break

               # Parse the molecule
               if flag == 1:
                   if len(fields) == 0: continue
                   atom_count += 1
                   if atom_count == n_atom:
                       flag = 2
                       continue
    qpos = "{}/fitcharge.qpos".format(data_dir)
    with open(qpos,'w') as g:
       for i in samples:
          for count_j,j in enumerate(geo):
              g.write('{:< 15.6f} {:< 15.6f} {:< 15.6f}\n'.format(j[0],j[1],j[2]))
          filename = "{}/{}/{}/{}/pointcharges.pc".format(run_dir,calc_folder,config_ori_folder,i) 
          with open(filename,'r') as f:
              for lc,lines in enumerate(f):
                  fields = lines.split()

                  # Initialize lists/arrays
                  if lc == 1: 
                      g.write('{:< 15.6f} {:< 15.6f} {:< 15.6f}\n'.format(float(fields[1]),float(fields[2]),float(fields[3])))

    
    return

def write_datadir(tag):
    with open(run_dir+'/datadir{}.def'.format(tag),'w') as f:
       f.write("* CHARMM data directory assignment\n"+\
               "*\n"+\
               "faster on\n"+\
               "set pnode =\n"+\
               "if ?numnode .gt. 1 set pnode = node 0\n"+\
               "set 0 data{}/     ! input data directory\n".format(tag)+\
               "set 9 scratch/  ! scratch directory\n"+\
               "set testcheck stream @0/test.str\n"+\
               "set qcheck stream @0/qtest.str\n"+
               "set testfail 0\n"+\
               "return\n")
    return

def tar_file(folder):
    os.chdir(run_dir+'/'+folder)
    ori_dir = os.getcwd()

    folder_names = ['configs_-0.5','configs_0.5','configs_connolly_-0.5','configs_connolly_0.5']
    for i in folder_names:
      subprocess.call('tar zcvf {}.tar.gz {}'.format(i,i),shell=True)
      subprocess.call('rm -r  {}'.format(i),shell=True)
    subprocess.call('tar zcvf pointcharge.tar.gz pointcharge_layer*'.format(i,i),shell=True)
    subprocess.call('rm -r  pointcharge_layer*',shell=True)
     
    os.chdir(ori_dir)

    return

def clean_pcsubmit(nconf):
   for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      delete_wildcard(fitfolder,'pc_layer')

   return

def clean_chelpg(nconf):
   folder = ['fitcharge_{}'.format(i) for i in range(0,nconf)]
   folder.append('unperturbed')
   for count_i,i in enumerate(folder):
      for q in [-0.5,0.5]:
         configfolder = 'configs_{}'.format(q)
         delete_wildcard('{}/{}'.format(i,configfolder),'pc_chelpg',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'prop',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'property.txt',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.in',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.out',addlayer=True)

   return

def clean_connolly(nconf):
   folder = ['fitcharge_{}'.format(i) for i in range(0,nconf)]
   folder.append('unperturbed')
   for count_i,i in enumerate(folder):
      for q in [-0.5,0.5]:
         configfolder = 'configs_connolly_{}'.format(q)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.in',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.out',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'K.tmp',addlayer=True)

   return

def clean_polar(nconf):
   for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      delete_wildcard('{}'.format(fitfolder),'pc_chelpg')
      delete_wildcard('{}'.format(fitfolder),'prop')
      delete_wildcard('{}'.format(fitfolder),'property.txt')
      delete_wildcard('{}'.format(fitfolder),'gbw')
      delete_wildcard('{}'.format(fitfolder),'vpot')
   return

def delete_wildcard(folder,pattern,addlayer=False):
    if addlayer:
      fileList = glob.glob('{}/{}/*/*{}*'.format(run_dir,folder,pattern))
    else:
      fileList = glob.glob('{}/{}/*{}*'.format(run_dir,folder,pattern))
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

def read_alljson(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    return jload

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
