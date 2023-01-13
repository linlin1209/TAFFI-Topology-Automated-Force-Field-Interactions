#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,glob,getpass,json
import numpy as np

sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
from parse_all import parse_configuration
from file_parsers import xyz_parse
from adjacency import Table_generator
from id_types import id_types
import place_charge,db_to_charmm,merge_FF_drude,relax_drude_loop,log_to_db


def main(argv):
    

    global user,c,intraonly,chargeonly,VDWonly,driver_dir,nonpolar_dir,args
    user = getpass.getuser()
    driver_dir = os.getcwd()


    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-alpha_scale', dest='alpha_scale', default=1,
                        help = 'options for charmm to rescale atomic polarizability when fitting, default: 1')
    parser.add_argument('-outFF', dest='outFF', default='fit_CHARMM.db',
                        help = 'output db filename, default: fit_CHARMM.db')
    parser.add_argument('-nonpolar_FF', dest='nonpolar_FF', default='',
                        help = 'nonpolar_FF db filename')
    parser.add_argument('-nonpolar_dir', dest='nonpolar_dir', default='',
                        help = 'nonpolar folder name')
    parser.add_argument('-o',dest='output',default='polar.log',
                        help = 'log name for stdout')
    parser.add_argument('-tag',dest='tag',default='',
                        help = '_tag name default: None')

    args=parser.parse_args()    
    sys.stdout = Logger(driver_dir+'/'+args.output)
    
    # parse configuration dictionary (c)
    print("parsing configurations...")
    c = parse_configuration(args.config)
    # Make config file absolute
    args.config = os.path.abspath(args.config)
    if (c['module_string'] == None):
      c['module_string'] = ''

    if args.nonpolar_FF == '':
      print('ERROR: nonpolar_FF file not specified...,Exiting....')

    if args.nonpolar_dir == '':
      print('ERROR: nonpolar directory not specified...,Exiting....')
    else:
      nonpolar_dir = args.nonpolar_dir

    gather_single()
    #fit_Drude(args,c,args.tag)

    return

# Normal Drude fitting, fit all heavy atoms
def fit_Drude(c,tag):

    global intraonly,chargeonly,VDWonly 
    if tag != '': tag = '_{}'.format(tag) 

    seperate_db(args.nonpolar_FF)
    intraonly = '{}/intra_only.db'.format(driver_dir)  
    chargeonly = '{}/charge_only.db'.format(driver_dir)
    VDWonly = '{}/vdw_only.db'.format(driver_dir)
    if os.path.isfile(intraonly) is False:
      print("ERROR: {} not found. Exiting...".format(intraonly))
      quit()
    if os.path.isfile(chargeonly) is False:
      print("ERROR: {} not found. Exiting...".format(chargeonly))
      quit()
    if os.path.isfile(VDWonly) is False:
      print("ERROR: {} not found. Exiting...".format(VDWonly))
      quit()
    FinalFF = '{}/final_drude{}.db'.format(driver_dir,tag)
    shutil.copyfile(intraonly,FinalFF)

    # FF dependency:
    # input for:
    # 1. fitcharge.inp: intra_only,VDW_only,charge_only(for initial guess),lower gens' Drude
    # 2. VDW refit: extract_vdw: lower gens' local Drude, Drude charge just fitted
    #               initial_param.db (initial guess): original cycle-20 and lower gens' local Drude (cross term is needed for initial guess)
    # FF info:
    # 1. in each inchikey:
    #    inchikey/drude_Q.db <- Drude charge
    #    inchikey/vdw/vdw_drude/DFT-AA.db <- Drude reparametrize VDW (has cross term)
    # 2. in driver folder:
    #    genN_drude.db : Drude charge, Drude vdw for that generation (has cross term)
    #    total_drude: complete FF (no cross term)

    # Read all.json (contains model compound and dependency information for the batch 
    status = read_alljson('./all.json')

    gens = sorted(status["gens"].keys())
    gens = [ int(i) for i in gens]

    for g in gens: 
      shutil.copyfile(intraonly,'{}/gen{}_drude.db'.format(driver_dir,g))
      option = {'master':'{}/gen{}_drude.db'.format(driver_dir,g),'no_ct':False}

      mc_keys = [ _ for _ in status["mc"].keys() if int(g) == status["mc"][_]["gen"][0] ] 
  
      jobids = [] 
      #################################################################
      #  Fit charge for all configurations at same generation at once #
      #################################################################
      for inchi in mc_keys:
         run_dir = '{}/{}'.format(driver_dir,inchi)
         c_config = '{}/{}/charges/CHELPG_calcs/configs'.format(nonpolar_dir,inchi)
         nconf = sorted([int(f) for f in os.listdir(c_config) if os.path.isdir(os.path.join(c_config, f))])[-1]+1
         if nconf > 100: nconf = 100

         for q in [-0.5,0.5]:
            for i in range(0,nconf):
               fitfolder = 'fitcharge_{}'.format(i)
               c_xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
               total_pc = get_pc(c_xyz,'{}/{}/data_{}'.format(run_dir,fitfolder,q))
               write_fitinp(run_dir,fitfolder,total_pc,q,g,'_{}'.format(q),c_xyz)
               fit_exe(run_dir+'/'+fitfolder,'_{}'.format(q),'_{}'.format(int(q*10)))

      os.chdir(driver_dir) 
      for q in [-5,5]:
         substring = "python3 {}/Automation_Scripts/bundle_submit.py  -f fit_{}.in -p {} -t {} -ppn {} -q {} -sched {} -size {} -o fit{} -path_to_exe {} --silent "
         substring = substring.format(c["taffi_path"],q,1,c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_ppn"],q,c["orca_exe"])
         tmp_job = execute_orca_submit(substring)
         jobids += tmp_job
      monitor_jobs(jobids,user)

      #######################################################
      # calculate average charge over all configurations    #
      #######################################################

      for inchi in mc_keys:
         run_dir = '{}/{}'.format(driver_dir,inchi)
         c_config = '{}/{}/charges/CHELPG_calcs/configs'.format(nonpolar_dir,inchi)
         nconf = sorted([int(f) for f in os.listdir(c_config) if os.path.isdir(os.path.join(c_config, f))])[-1]+1
         if nconf > 100: nconf = 100
         xyz = '{}/{}/{}.xyz'.format(c_config,0,0)
         element,geo = xyz_parse(xyz)
         adj_mat = Table_generator(element,geo)
         atom_types = sorted(list(set(id_types(element,adj_mat,int(c['gens'])))),reverse=True)
         filename = '{}/chi2_Q.txt'.format(run_dir)
         f = open_chi(filename,atom_types)
         charge_dict = {}
         polar_dict = {}
         for q in [-0.5,0.5]:
            for count_i,i in enumerate(list(range(0,nconf))):
              fitfolder = 'fitcharge_{}'.format(i)
              title = 'fit_{}_{}'.format(q,i)
              c_xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
              total_pc = get_pc(c_xyz,'{}/{}/data_{}'.format(run_dir,fitfolder,q))
              os.chdir(run_dir+'/'+fitfolder)
              charge_dict,polar_dict = write_chi(title,q,xyz,total_pc,atom_types,f,charge_dict,polar_dict,c['gens'])
         f.close()
         dbname = '{}/drude_Q.db'.format(run_dir)
         for key in charge_dict:
            charge_dict[key] /= (nconf*2)
            polar_dict[key] /= (nconf*2)
         write_db(dbname,charge_dict,polar_dict)
         VDW_re(run_dir,args,tag,inchi,g,xyz)
         option['new_params'] = '{}/drude_Q.db'.format(run_dir)
         merge_FF_drude.fun(option)
         option['new_params'] = '{}/{}/vdw/vdw_drude{}/DFT-AA.db'.format(nonpolar_dir,inchi,tag)
         merge_FF_drude.fun(option)

    ### Final Merge
    option = {'master':FinalFF,'no_ct':True}
    for g in gens: 
      option['new_params'] = '{}/gen{}_drude.db'.format(driver_dir,g)
      merge_FF_drude.fun(option)
   
    return

def gather_single():
      inchi = 'OKKJLVBELUTLKV-UHFFFAOYSA-N'
      run_dir = '/scratch/bell/lin1209/Drude_benchmark_AMOEBA/OKKJLVBELUTLKV-UHFFFAOYSA-N/O_only'
      nonpolar_dir = '/scratch/bell/lin1209/Drude_benchmark_AMOEBA/'
      c_config = '{}/{}/charges/CHELPG_calcs/configs'.format(nonpolar_dir,inchi)
      nconf = sorted([int(f) for f in os.listdir(c_config) if os.path.isdir(os.path.join(c_config, f))])[-1]+1
      if nconf > 100: nconf = 100
      xyz = '{}/{}/{}.xyz'.format(c_config,0,0)
      #xyz = '{}/{}/geo_opt.xyz'.format(driver_dir,inchi)
      element,geo = xyz_parse(xyz)
      adj_mat = Table_generator(element,geo)
      atom_types = sorted(list(set(id_types(element,adj_mat,int(c['gens'])))),reverse=True)
      filename = '{}/chi2_Q.txt'.format(run_dir)
      f = open_chi(filename,atom_types)
      charge_dict = {}
      polar_dict = {}
      for q in [-0.5,0.5]:
         for count_i,i in enumerate(list(range(0,nconf))):
           fitfolder = 'fitcharge_{}'.format(i)
           title = 'fit_{}_{}'.format(q,i)
           c_xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
           total_pc = get_pc(c_xyz,'{}/{}/data_{}'.format(run_dir,fitfolder,q))
           #total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(run_dir,fitfolder)
           os.chdir(run_dir+'/'+fitfolder)
           #print(run_dir+'/'+fitfolder)
           charge_dict,polar_dict = write_chi(title,q,xyz,total_pc,atom_types,f,charge_dict,polar_dict,c['gens'])
      f.close()
      dbname = '{}/drude_Q.db'.format(run_dir)
      for key in charge_dict:
         charge_dict[key] /= (nconf*2)
         polar_dict[key] /= (nconf*2)
      write_db(dbname,charge_dict,polar_dict)

      return

# Seperate a db file into:
# 1. intra_only
# 2. charge_only
# 3. vdw_only 
# db files
def seperate_db(filename):
    g = open('intra_only.db','w')
    h = open('charge_only.db','w')
    l = open('vdw_only.db','w')
    with open(filename,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) > 0 and fields[0] == 'vdw':
            l.write(lines) 
         elif len(fields) >0 and (fields[0] == 'bond' or fields[0] == 'angle' or fields[0] == 'torsion'):
            g.write(lines)
         elif len(fields) >0 and fields[0] == 'charge':
            h.write(lines)
         else:
            g.write(lines)
            h.write(lines)
            l.write(lines)
    g.close()
    h.close()
    l.close()

    return

def read_alljson(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    return jload

def get_pc(xyz,datafolder):
    
    # get geometry length
    atom_count = 0
    with open(xyz,'r') as f:
      for lc,line in enumerate(f):
         if lc != 0 or lc != 1: atom_count += 1
    dataname = '{}/fitcharge.qpos'.format(datafolder)
    line_count = 0
    with open(dataname,'r') as f:
      for lc,line in enumerate(f):
         line_count += 1

    total_pc = line_count // (atom_count+1)
   
    return total_pc

# Wrtie CHARMM fitcharge script
def write_fitinp(run_dir,fitfolder,npert,mag_c,current_gen,tag,c_xyz):

    ori_dir = os.getcwd()
    os.chdir(run_dir)

    # need intra_only, charge_only, VDW_only, lower gens' local Drude
    #genN.db has charges, VDW drude parameters for that generation (cross term included)

    # XXX: hard coding for mixing_rule
    if os.path.isfile(fitfolder+'/fitcharge{}.inp'.format(tag)) is False:
      sublist = ('python3 {}/FF_functions/gen_md_for_CHARMM.py {}  -o {}  -gens {} -npert {} -ascale {} -mixing_rule wh --fit -mag_c {} -inpname fitcharge{} -FF '.format(c["taffi_path"],c_xyz,fitfolder,c['gens'],npert,args.alpha_scale,mag_c,tag)).split()
      FFstring = '{} {} {} '.format(intraonly,chargeonly,VDWonly)
      # append lower gens' local Drude FF
      for i in range(1,current_gen):
         FFstring += '{}/gen{}_drude.db '.format(driver_dir,i)
      sublist.append(FFstring)
      shell_exe(sublist,list_string=True,idle_time=0.1)
    os.chdir(ori_dir)

    return

def shell_exe(command,list_string=False,idle_time=3):
    # if the input command is not a list but a complete string, split it
    if not list_string:
      command = command.split()
    process = subprocess.Popen(command)
    status = process.poll()
    while status == None:
       status = process.poll()
       time.sleep(idle_time)
    return

def place_perturbation(run_dir):

    # these folders should already be created during "drude_config.py", but in case they're tarred due to file amount issue, it'd call place_charge 

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/fitcharge')

    # Place perturantion charges
    if os.path.isdir('pointcharge_layer1'):
       pc_layer1 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer1') if os.path.isfile(os.path.join('pointcharge_layer1', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       sublist = ('-d 500 -radii_factor 13 -c 1 -extra 0.1 -o layer1 fitcharge.xyz -ofolder pointcharge_layer1').split()
       pc_layer1 = place_charge.main(sublist)
    if os.path.isdir('pointcharge_layer2'):
       pc_layer2 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer2') if os.path.isfile(os.path.join('pointcharge_layer2', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       sublist = ('-d 500 -radii_factor 22 -c 1 -extra 0.02 -o layer2 fitcharge.xyz -ofolder pointcharge_layer2').split()
       pc_layer2 = place_charge.main(sublist)
    if os.path.isdir('pointcharge_layer3'):
       pc_layer3 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer3') if os.path.isfile(os.path.join('pointcharge_layer3', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       sublist = ('-d 500 -radii_factor 40 -c 1 -extra 0.003 -o layer3 fitcharge.xyz -ofolder pointcharge_layer3').split()
       pc_layer3 = place_charge.main(sublist)
   
    total_pc = pc_layer1 + pc_layer2 + pc_layer3

    print("Successfully place {} perturbations charges, layer1: {}, layer2: {}, layer3: {}".format(total_pc,pc_layer1,pc_layer2,pc_layer3))
   
    os.chdir(ori_dir)

    return total_pc,pc_layer1,pc_layer2,pc_layer3

def fit_exe(fit_dir,tag,intag):
    with open('{}/fit{}.in'.format(fit_dir,intag),'w') as f:
         f.write("#!/bin/bash\n")
         f.write("charmm < fitcharge{}.inp > fit{}.out\n".format(tag,intag))
    return

def execute_orca_submit(substring,sublist=False):
    if sublist:
      output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    else:
      output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    orcaids = [ m.split()[-1] for m in output.split("\n")[:-1]]

    return orcaids

def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

def open_chi(filename,atom_types):
    f = open(filename,'w')
    f.write("{:<20s} {:<20s} {:<20s} {:<20s}".format("title","chi2","molecular_polar_MM","molecular_polar_QM"))
    for i in atom_types:
      if i.split('[')[1] != '1':
        f.write("{:<20s} {:<20s}".format(i,'polar'))
      else:
        f.write("{:<20s} ".format(i))
    f.write("\n")
    return f

def write_chi(title,q,xyz,total_pc,atomtypes,f,charge_dict,polar_dict,gens):
    command = ('-logname fit_{}.out -o {}.db -xyz {} -gens {}'.format(int(q*10),title,xyz,gens)).split()
    tmpchi,tmp_molpolar = log_to_db.main(command)
    QM_dipole,QM_molpolar = read_polarout("polar.out")  
    tmp = read_db('{}.db'.format(title))
    f.write("{:<20s} {:<20.5f} {:<20.5f} {:<20.5f} ".format('{}'.format(title),tmpchi/total_pc,tmp_molpolar,QM_molpolar))
    for j in atomtypes:
         if j.split('[')[1] != '1':
            f.write("{:<20.5f} {:<20.5f} ".format(tmp[0][j],tmp[1][j])) 
         else:
            f.write("{:<20.5f} ".format(tmp[0][j])) 
    f.write("\n")
    if charge_dict == {}:
         charge_dict = tmp[0]
         polar_dict = tmp[1]
    else:
         for key in tmp[0]:
            charge_dict[key] += tmp[0][key]
            polar_dict[key] += tmp[1][key]
         
    return charge_dict,polar_dict

def write_db(filename,charge_dict,polar_dict):
    with open(filename,'w') as f:

        f.write("\n# Charge definitions\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","Charge"))
        for key in charge_dict:
            f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(key),charge_dict[key]))
        f.write("\n# Atomic Polarizability\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","polarizability(A^3)"))
        for key in polar_dict:
            f.write("{:10s} {:40s} {:< 20.6f}\n".format("polar",str(key),polar_dict[key]))
    return

def VDW_re(run_dir,args,tag,inchikey,current_gen,xyz):

    ori_dir = os.getcwd()
    print("Start vdw reparametrization...")
    # move to non_polar folder to relax Drude particles
    vdwfolder = nonpolar_dir + '/' + inchikey + '/vdw' 
    vdwdrude_dir = 'vdw_drude{}'.format(tag) 

    # if completed then skip
    if os.path.isfile('{}/{}/DFT-AA.db'.format(vdwfolder,vdwdrude_dir)) is True:
      print('{}/{}/DFT-AA.db found, skipping....'.format(vdwfolder,vdwdrude_dir))
      print("COMPLETE: reparametrize vdw")
      return
      

    if os.path.isdir(vdwfolder) is False:
      print("ERROR: nonpolar vdw folder {} doesn't exist. Exiting...".format(vdw))
    os.chdir(vdwfolder)
    print("move to {} to reparametrize...".format(vdwfolder))
    # converting the new db file to charmm format FF which is needed when relaxing via CHARMM
    FFname = 'FF{}'.format(tag)
    if os.path.isdir(vdwfolder+'/'+FFname) is False:
      # XXX hard coding mixing_rule
      sublist = ('{}  -o {}  -gens {} -mixing_rule wh --polar -FF '.format(xyz,FFname,c['gens'])).split()
      # append lower gens' local Drude FF
      FFstring = '{} {} {} {} {}/drude_Q.db'.format(intraonly,chargeonly,VDWonly,' '.join(['{}/gen{}_drude.db'.format(driver_dir,i) for i in range(1,current_gen)]),run_dir)
      sublist.append(FFstring)
      # XXX need to figure out whether do an outside call or a function call is more efficient
      db_to_charmm.main(sublist)

    ## XXX Added if configs folder has been compressed
    subprocess.call('tar zxvf configs.tar.gz',shell=True)

    os.chdir('{}/configs'.format(vdwfolder))
    relax_drude_loop.fun(args.tag)

    # initial_param.db serves as extract_vdw's initial guess
    # normally this would be the FF from the cycle before the very last cycle
    # we want to use the fully parametrized FF as our initial guess
    os.chdir(vdwfolder)
    shutil.copyfile('{}/cycle-20/DFT-AA.db'.format(vdwfolder),'{}/initial_params.db'.format(vdwfolder))

    option = {'master':'{}/initial_params.db'.format(vdwfolder),'no_ct':False,'replace':True}
    for i in range(1,current_gen):
      option['new_params'] = '{}/gen{}_drude.db'.format(driver_dir,i)
      merge_FF_drude.fun(option)
    
    # mv to inchikey folder to execute extract_vdw
    os.chdir('{}/{}'.format(nonpolar_dir,inchikey))

    L2_s = "0.1" 
    L2_e = "0.5"
    if os.path.isdir(vdwfolder+'/'+vdwdrude_dir) is False:
       # input for existing data base is genN.db (Drude vdw only)
       # initial_parms serves as initial guess, a mix of taffi VDW and Drude VDW
       FF_path = ' '.join(['{}/gen{}_drude.db'.format(driver_dir,i) for i in range(1,current_gen)]) 
       charge_FF = '{}/drude_Q.db'.format(run_dir) 
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       if tag == '':
         sublist.append('{}/FF_functions/extract_vdw_drude.py -f vdw -FF_DFT \'{} {}\' -o {} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} '.format(c['taffi_path'],FF_path,charge_FF,vdwdrude_dir,L2_s,L2_e,c['functional']))
       else:
         sublist.append('{}/FF_functions/extract_vdw_drude.py -f vdw -FF_DFT \'{} {}\' -o {} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} -tag {}'.format(c['taffi_path'],FF_path,charge_FF,vdwdrude_dir,L2_s,L2_e,c['functional'],tag))
         
       sublist= sublist + ('vdw_drude_parse  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       jobids = execute_orca_submit(sublist,sublist=True)
       # Wait until jobs complete
       monitor_jobs(jobids,user)

    if os.path.isfile(vdwfolder+'/'+vdwdrude_dir+'/DFT-AA.db') is False:
       print("ERROR: vdw/{}/DFT-AA.db not found, something went wrong w/ vdw_parse_drude".format(vdwdrude_dir))
       quit()

    else:
       ## XXX Added if configs folder has been compressed
       os.chdir(vdwfolder)
       subprocess.call('mkdir configs_drude',shell=True)
       subprocess.call('mv configs/*/*_drude.xyz configs_drude',shell=True)
       subprocess.call('tar zcvf configs_drude.tar.gz configs_drude',shell=True)
       subprocess.call('rm -r configs_drude',shell=True)
       subprocess.call('rm -r configs',shell=True)
       
       print("COMPLETE: reparametrize vdw")

    os.chdir(ori_dir)

    return

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
