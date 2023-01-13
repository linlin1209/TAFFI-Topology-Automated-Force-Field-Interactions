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
import log_to_db,merge_FF_drude,db_to_charmm,relax_drude_loop
from monitor_jobs import *
from parse_all import *
from file_parsers import *
import place_charge,generate_connolly,generate_connolly_new,gen_md_for_CHARMM_new
import codecs,json
# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *


def main(argv):
    

    global user,c,intraonly,chargeonly,VDWonly,driver_dir,nonpolar_dir
    user = getpass.getuser()
    driver_dir = os.getcwd()
    nonpolar_dir = '/home/lin1209/fhalstead/nitrogen_cycle20_N100/parameterization110320_cycle20_config100/'

    status = read_alljson('./all.json')

    gens = sorted(status["gens"].keys())
    gens = [ int(i) for i in gens]
     
    h = open('chi2_sum_tot.txt','w')
    h.write("{:<30s} {:<20s} {:<10s} {:<10s} {:<10s}\n".format('inchikey','smiles','ori','initial','final'))

    for g in gens: 
      mc_keys = [ _ for _ in status["mc"].keys() if int(g) == status["mc"][_]["gen"][0] ] 
      for mc in mc_keys:
         flag = 0
         if os.path.isfile('{}/{}/vdw/vdw_drude/extract_vdw.log'.format(nonpolar_dir,mc)) is False:
            print("{} vdw drude not found".format(mc))
            continue
         initial = 0
         final = 0
         ori = 0
         count = 0
         with open('{}/{}/vdw/vdw_drude/figures/DFT-AA_tot_convergence.txt'.format(nonpolar_dir,mc)) as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0: continue
               initial += (float(fields[1]) - float(fields[0]))**2
               final += (float(fields[2]) - float(fields[0]))**2
               count += 1
            initial /= count
            final /= count
         count = 0
         with open('{}/{}/vdw/cycle-20/figures/DFT-AA_tot_convergence.txt'.format(nonpolar_dir,mc)) as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0: continue
               ori += (float(fields[2]) - float(fields[0]))**2
               count += 1
            ori /= count
         smiles = get_smiles('{}/{}/{}.xyz'.format(nonpolar_dir,mc,mc))
         h.write("{:<30s} {:<20s} {:<10.6f} {:<10.6f} {:<10.6f}\n".format(mc,smiles,ori,initial,final))

    return

def get_smiles(xyz):

   p = subprocess.Popen("obabel -ixyz {} -osmiles".format(xyz).split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
   smiles, err = p.communicate()
   smiles = str(smiles, 'utf-8').split('\t')[0] ##might be a problem if there are characters utf-8 can't handle

   return smiles

# Normal Drude fitting, fit all heavy atoms
def fit_Drude(args,tag):

    global intraonly,chargeonly,VDWonly 
    if tag != '': tag = '_{}'.format(tag) 

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
    #               initial_param.db (initial guess): original cycle-2 and lower gens' local Drude (cross term is needed for initial guess)
    # FF info:
    # 1. in each inchikey:
    #    inchikey/drude_Q.db <- Drude charge
    #    inchikey/vdw/vdw_drude/DFT-AA.db <- Drude reparametrize VDW (has cross term)
    # 2. in driver folder:
    #    genN.db : Drude charge, Drude vdw for that generation (has cross term)
    #    total_drude: complete FF (no cross term)

    # Read all.json (contains model compound and dependency information for the batch 
    status = read_alljson('./all.json')

    gens = sorted(status["gens"].keys())
    gens = [ int(i) for i in gens]

    for g in gens: 
      shutil.copyfile(intraonly,'{}/gen{}.db'.format(driver_dir,g))
      option = {}
      option['master'] = '{}/gen{}.db'.format(driver_dir,g) 
      option['no_ct'] = False
      mc_keys = [ _ for _ in status["mc"].keys() if int(g) == status["mc"][_]["gen"][0] ] 
      if g==4:
         mc_keys = ['OFBQJSOFQDEBGM-UHFFFAOYSA-N','KVNRLNFWIYMESJ-UHFFFAOYSA-N','HPNMFZURTQLUMO-UHFFFAOYSA-N','UAOMVDZJSHZZME-UHFFFAOYSA-N','HQABUPZFAYXKJW-UHFFFAOYSA-N','ZMANZCXQSJIPKH-UHFFFAOYSA-N','PQPFFKCJENSZKL-UHFFFAOYSA-N','GVWISOJSERXQBM-UHFFFAOYSA-N'] 
  
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
               total_pc = get_pc_new(c_xyz,'{}/{}/data_{}'.format(run_dir,fitfolder,q))
               #total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(run_dir,fitfolder)
               write_fitinp(run_dir,args,fitfolder,total_pc,q,g,'_{}'.format(q),c_xyz)
               fit_exe(run_dir+'/'+fitfolder,'_{}'.format(q),'_{}'.format(int(q*10)))

      os.chdir(driver_dir) 
      for q in [-5,5]:
         substring = "python3 {}/Automation_Scripts/bundle_submit.py  -f fit_{}.in -p {} -t {} -ppn {} -q {} -sched {} -size {} -o fit{} -path_to_exe {} --silent "
         substring = substring.format(c["taffi_path"],q,1,c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_ppn"],q,c["orca_exe"])
         #substring = substring.format(c["taffi_path"],1,c["charges_qc_wt"],c["charges_qc_ppn"],'bsavoie',c["charges_qc_sched"],c["charges_qc_ppn"],c["orca_exe"])
         output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
         output = str(output,'utf-8')
         jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
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
         #xyz = '{}/{}/geo_opt.xyz'.format(driver_dir,inchi)
         element,geo = xyz_parse(xyz)
         adj_mat = Table_generator(element,geo)
         atom_types = sorted(list(set(id_types(element,adj_mat,2))),reverse=True)
         filename = '{}/chi2_Q.txt'.format(run_dir)
         f = open_chi(filename,atom_types)
         charge_dict = {}
         polar_dict = {}
         for q in [-0.5,0.5]:
            for count_i,i in enumerate(list(range(0,nconf))):
              fitfolder = 'fitcharge_{}'.format(i)
              title = 'fit_{}_{}'.format(q,i)
              c_xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
              total_pc = get_pc_new(c_xyz,'{}/{}/data_{}'.format(run_dir,fitfolder,q))
              #total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(run_dir,fitfolder)
              os.chdir(run_dir+'/'+fitfolder)
              #print(run_dir+'/'+fitfolder)
              charge_dict,polar_dict = write_chi(title,q,xyz,total_pc,atom_types,f,charge_dict,polar_dict)
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
    option = {}
    option['no_ct'] = True
    option['master'] = FinalFF
    for g in gens: 
      option['new_params'] = '{}/gen{}.db'.format(driver_dir,g)
      merge_FF_drude.fun(option)
   
    quit()
         
               

    # After charge and vdw parameterization
    # localFF has Intra, Drude charge, Drude VDW including cross term
    # this is the input for higher generation extract_vdw, so that same VDW wouldn't be refitted
    if tag == '':
      localFF = '{}/total_drude.db'.format(run_dir)
    else:
      localFF = '{}/total_drude_{}.db'.format(run_dir,tag)
    # FinalFF is all localFF combined, but w/o cross term
    # This is used to write fitcharge.inp
    if tag == '':
      FinalFF = '{}/final_drude.db'.format(driver_dir)
    else:
      FinalFF = '{}/final_drude_{}.db'.format(driver_dir,tag)

    # initial guess for charges (from origianl TAFFI)
    charge_IG = '{}/charge_only.db'.format(driver_dir)

    # copy original intra_only.db to total_drude.db
    if os.path.isfile(driver_dir+'/intra_only.db') is False:
         print("{}/intra_only.db doesn't exist! Exiting....".format(driver_dir))
         quit()
    print("Copy intra_only.db to {}.....".format(localFF))
    shutil.copyfile(driver_dir+"/intra_only.db",localFF)

    write_fitinp(c,run_dir,FinalFF,args,fitfolder,charge_IG)

    total_pc,pc_layer1,pc_layer2,pc_layer3 = place_perturbation(run_dir)
    
    element,geo = xyz_parse(run_dir+"/fitcharge/fitcharge.xyz")

    change_npert(run_dir,total_pc,fitfolder)

    #calc_CHELPG(run_dir,total_pc,pc_layer1,pc_layer2,pc_layer3,c,geo,element,user)

    #connolly_grid = get_connolly_grid(run_dir)

    #calc_connolly_pot(run_dir,total_pc,connolly_grid,c,user)

    write_charmm_data(run_dir,total_pc,geo,fitfolder)

    Drude_charge(run_dir,FinalFF,localFF,args,c,user,fitfolder)

    VDW_re(element,geo,FinalFF,localFF,run_dir,c,user,args,tag,fitfolder)

   
    return

def write_db(filename,charge_dict,polar_dict):
    with open(filename,'w') as f:

        f.write("\n# Charge definitions\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","Charge"))
        for key in charge_dict:
            f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(key),charge_dict[key]))
        f.write("\n# Atomic Polarizability\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","polarizability(A^3)"))
        for key in polar_dict:
            f.write("{:10s} {:40s} {:< 20.6f}\n".format("polar",str(key),polar_dict[key]))
    return


def read_db(dbfile):
    with open(dbfile,'r') as f:
      charge_dict = {}
      polar_dict = {}
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) == 0: continue
         if fields[0] == 'charge':
            charge_dict[fields[1]] = float(fields[2])
         if fields[0] == 'polar':
            polar_dict[fields[1]] = float(fields[2])
    return (charge_dict,polar_dict)

def read_polarout(filename):
    with open(filename,'r') as f:
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 2 and flag == 0: continue
          
            if flag == 0  and fields[0] == "Isotropic" and fields[1] == "polarizability" :
               polar = float(fields[3])*0.529177210**3 # A.U to A^3
               continue
            
            if flag == 0  and fields[0] == 'DIPOLE' and fields[1] == 'MOMENT' :
               flag = 1
               continue


            if flag == 1 and len(fields) == 4 and fields[0] == 'Magnitude' and fields[1] == '(Debye)':
               dipole = float(fields[3])
               flag = 0 
               continue
    return dipole,polar

def write_chi(title,q,xyz,total_pc,atomtypes,f,charge_dict,polar_dict):
    command = ('-logname fit_{}.out -o {}.db -xyz {}'.format(int(q*10),title,xyz)).split()
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

def fit_exe(fit_dir,tag,intag):
    with open('{}/fit{}.in'.format(fit_dir,intag),'w') as f:
         f.write("#!/bin/bash\n")
         f.write("charmm < fitcharge{}.inp > fit{}.out\n".format(tag,intag))
    #sublist = ('sh fit.sh').split()
    #output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]

    return

def get_pc_new(xyz,datafolder):
    
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


def get_pc(run_dir,calc_folder):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    if os.path.isdir('pointcharge_layer1'):
       pc_layer1 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer1') if os.path.isfile(os.path.join('pointcharge_layer1', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 1 not placed for {}".format(calc_folder))
    if os.path.isdir('pointcharge_layer2'):
       pc_layer2 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer2') if os.path.isfile(os.path.join('pointcharge_layer2', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 2 not placed for {}".format(calc_folder))
    if os.path.isdir('pointcharge_layer3'):
       pc_layer3 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer3') if os.path.isfile(os.path.join('pointcharge_layer3', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 3 not placed for {}".format(calc_folder))
   
    total_pc = pc_layer1 + pc_layer2 + pc_layer3

    print("Successfully place {} perturbations charges for {}, layer1: {}, layer2: {}, layer3: {}".format(total_pc,calc_folder,pc_layer1,pc_layer2,pc_layer3))
   
    os.chdir(ori_dir)

    return total_pc,pc_layer1,pc_layer2,pc_layer3


# Wrtie CHARMM fitcharge script
def write_fitinp(run_dir,args,fitfolder,npert,mag_c,current_gen,tag,c_xyz):

    ori_dir = os.getcwd()
    os.chdir(run_dir)

    # need intra_only, charge_only, VDW_only, lower gens' local Drude
    #genN.db has charges, VDW drude parameters for that generation (cross term included)

    if os.path.isfile(fitfolder+'/fitcharge{}.inp'.format(tag)) is False:
      #sublist = ('geo_opt.xyz  -o {} -q {} -gens {} -npert 0 -ascale {} -mixing_rule wh --fit -inpname fitcharge -FF '.format(fitfolder,c['charge'],c['gens'],args.alpha_scale)).split()
      sublist = ('python3 {}/FF_functions/gen_md_for_CHARMM_new.py {}  -o {}  -gens {} -npert {} -ascale {} -mixing_rule wh --fit -mag_c {} -inpname fitcharge{} -FF '.format(c["taffi_path"],c_xyz,fitfolder,c['gens'],npert,args.alpha_scale,mag_c,tag)).split()
      #sublist = ('{}  -o {}  -gens {} -npert {} -ascale {} -mixing_rule wh --fit -mag_c {} -inpname fitcharge{} -FF '.format(c_xyz,fitfolder,c['gens'],npert,args.alpha_scale,mag_c,tag)).split()
      FFstring = '{} {} {} '.format(intraonly,chargeonly,VDWonly)
      # append lower gens' local Drude FF
      for i in range(1,current_gen):
         FFstring += '{}/gen{}.db '.format(driver_dir,i)
      sublist.append(FFstring)
      #subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      process = subprocess.Popen(sublist)
      status = process.poll()
      while status == None:
          status = process.poll()
          time.sleep(0.1)
      #gen_md_for_CHARMM_new.main(sublist)
    os.chdir(ori_dir)

    return

def place_perturbation(run_dir):

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

def place_perturbation_submit(run_dir,calc_folder,xyz,dplus=0):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    jobids = []
    if os.path.isdir('pointcharge_layer1') is False:
       d = 700+int(dplus)
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 13 -c 1 -extra 0.1 -o layer1 {} -ofolder pointcharge_layer1'.format(run_dir,calc_folder,c['taffi_path'],d,xyz))
       sublist= sublist + ('pc_layer1  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
    else:
       if len(os.listdir('pointcharge_layer1')) == 0:
         shutil.rmtree('pointcharge_layer1')
    if os.path.isdir('pointcharge_layer2') is False:
       d = 500 +int(dplus)
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 22 -c 1 -extra 0.02 -o layer2 {} -ofolder pointcharge_layer2'.format(run_dir,calc_folder,c['taffi_path'],d,xyz))
       sublist= sublist + ('pc_layer2  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
    else:
       if len(os.listdir('pointcharge_layer2')) == 0:
         shutil.rmtree('pointcharge_layer2')
    if os.path.isdir('pointcharge_layer3') is False:
       d = 500 +int(dplus)
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 40 -c 1 -extra 0.003 -o layer3 {} -ofolder pointcharge_layer3'.format(run_dir,calc_folder,c['taffi_path'],d,xyz)) 
       sublist= sublist + ('pc_layer3  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
    else:
       if len(os.listdir('pointcharge_layer3')) == 0:
         shutil.rmtree('pointcharge_layer3')

    print("Submitting place_charge jobs for {}...".format(calc_folder))
   
    os.chdir(ori_dir)

    return jobids

def change_npert(run_dir,total_pc,fitfolder):

    inpname = '{}/{}/fitcharge.inp'.format(run_dir,fitfolder)
    # Write # of perturnation charges back to CHARMM fit script
    with open(inpname, 'r') as f:
       # read a list of lines into data
       lines = f.readlines()

    # change npert
    for count_i,i in enumerate(lines):
       fields = i.split() 
       if len(fields) < 2: continue
       if fields[0] == 'NPERT':
         lines[count_i] = '    NPERT {} -             ! number of perturbation ions used in QM ESP\n'.format(total_pc)
         break

    # and write everything back
    with open(inpname, 'w') as f:
       f.writelines( lines )

    return

# Calculating CHELPG potential
def calc_CHELPG(run_dir,total_pc,pc_layer1,pc_layer2,pc_layer3,c,geo,element,user):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/fitcharge')

    print("Start calculating CHELPG potential...")
    if os.path.isdir("configs") is False:
         os.makedirs("configs")
    os.chdir(run_dir+"/fitcharge/configs")
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else:
            scfpfile='{}/{}_charges.scfp'.format(i,i)
            if  os.path.isfile(scfpfile) is True: 
               continue
         os.chdir(str(i))
         shutil.copyfile(run_dir+"/fitcharge/fitcharge.xyz",'{}.xyz'.format(i))
         if i<pc_layer1:
            shutil.copyfile(run_dir+"/fitcharge/pointcharge_layer1/pointcharges.{}.pc".format(i),'pointcharges.pc')
         elif i< pc_layer1+pc_layer2:
            shutil.copyfile(run_dir+"/fitcharge/pointcharge_layer2/pointcharges.{}.pc".format(i-pc_layer1),'pointcharges.pc')
         else:
            shutil.copyfile(run_dir+"/fitcharge/pointcharge_layer3/pointcharges.{}.pc".format(i-(pc_layer1+pc_layer2)),'pointcharges.pc')
           
         # Write orca input file   
         name = "{}_charges".format(i)
         Write_charge_in(c,name,geo,element)
         scfpfile='{}_charges.scfp'.format(i)
         os.chdir(run_dir+'/fitcharge/configs')

    # Create directory for unperturbed charge potential
    if os.path.isdir("unperturbed") is False: os.makedirs("unperturbed") 
    os.chdir("unperturbed")
    shutil.copyfile("{}/fitcharge/fitcharge.xyz".format(run_dir),'unperturbed.xyz')
    Write_charge_in(c,'unperturbed_charges',geo,element)
    os.chdir(run_dir+'/fitcharge/configs')

    # Submit orca job
    jobids = []
    substring = "python3 {}/Automation_Scripts/orca_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],c["charges_qc_procs"],c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_size"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    monitor_jobs(jobids,user)

    print("COMPLETE: CHELPG potential calculation")
    os.chdir(ori_dir)

    return

def get_connolly_grid(run_dir):

    ori_dir = os.getcwd()
    # Genererating connolly surface grid
    os.chdir(run_dir+'/fitcharge')
    if os.path.isfile('map.xyz') is False:
         sublist = ('fitcharge.xyz -d 1.41').split()
         connolly_grid = generate_connolly_new.main(sublist)
    else:
         with open('map.xyz','r') as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0:
                  connolly_grid = int(fields[0]) 
                  break
    os.chdir(ori_dir)

    return connolly_grid

def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

# Revaluating potential at connolly grid
def calc_connolly_pot(run_dir,total_pc,connolly_grid,c,user):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/fitcharge')

    print("Start evaluation of potential at connolly grid...")
    if os.path.isdir("configs_connolly") is False:
         os.makedirs("configs_connolly")
    os.chdir(run_dir+"/fitcharge/configs_connolly")
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile(str(i)+'/tmp.output') is True: 
               # write an out file so that bundle_submit will leave it alone
               with open(str(i)+'/charge.out','w') as f:
                  f.write("complete")
               continue
         os.chdir(str(i))
         
         # copy necessary file from configs folder
         prepare(run_dir,i)
         if os.path.isfile('tmp.output') is False:
            with open('charge.in','w') as f:
               f.write("#!/bin/bash\n")
               f.write("{} {}".format((c['module_string'].split()[-2]).replace('\\n','\n'),(c['module_string'].split()[-1]).replace('\\n','\n')))
               f.write("/depot/bsavoie/apps/orca_4_1_2/orca_vpot {}_charges.gbw {}_charges.scfp {}/map.xyz tmp.output\n".format(i,i,run_dir))
         os.chdir(run_dir+'/fitcharge/configs_connolly')
   
    # Submit orca job
    jobids = []
    substring = "python3 {}/Automation_Scripts/bundle_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],1,c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_ppn"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    monitor_jobs(jobids,user)
   
    # Rewrite orca_vpot output to CHELPG output format
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile(str(i)+'_charges.vpot') is True: continue
         os.chdir(str(i))

         oldpot = '{}/fitcharge/configs/{}/{}_charges.vpot'.format(run_dir,i,i)
         newpot = '{}_charges.vpot'.format(i)
         rewrite_potential(oldpot,newpot,connolly_grid)
         os.chdir(run_dir+'/fitcharge/configs_connolly')

    if os.path.isdir('unperturbed') is False: os.makedirs('unperturbed')
    os.chdir('unperturbed')
    i='unperturbed'
    prepare(run_dir,'unperturbed')
    if os.path.isfile('tmp.output') is False:
       with open('charge.sh','w') as f:
            f.write("#!/bin/bash\n")
            f.write("{} {}".format((c['module_string'].split()[-2]).replace('\\n','\n'),(c['module_string'].split()[-1]).replace('\\n','\n')))
            f.write("/depot/bsavoie/apps/orca_4_1_2/orca_vpot {}_charges.gbw {}_charges.scfp map tmp.output\n".format(i,i))
       command = 'sh charge.sh'.split()
       process = subprocess.Popen(command)
       status = process.poll()
       while status == None:
          status = process.poll()
          time.sleep(3)
    oldpot = '{}/fitcharge/configs/{}/{}_charges.vpot'.format(run_dir,i,i)
    newpot = '{}_charges.vpot'.format(i)
    rewrite_potential(oldpot,newpot,connolly_grid)
    print("COMPLETE: revaluation of potential at connolly grid")
    os.chdir(run_dir+'/fitcharge/')
    os.chdir(ori_dir)

    return

def write_charmm_data(run_dir,total_pc,geo,fitfolder):
   
    data_dir = '{}/{}/data'.format(run_dir,fitfolder)
    # Rewrite everything to CHARMM compatible format
    if os.path.isdir(data_dir) is False:
       os.makedirs(data_dir)

    Ang2Bohr = 1.889725989 
    # ESP unit in CHARMM: e/A
    # ESP unit in ORCA: e/Bohr
    # length unit: CHARMM: Angstrom, ORCA: Bohr
    fitpot = "{}/{}/data/fitcharge.pot".format(run_dir,fitfolder)
    with open(fitpot,'w') as g:
       for i in range(0,total_pc):
          filename = "{}/fitcharge/configs_connolly/{}/{}_charges.vpot".format(run_dir,i,i) 
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
                  if flag == 2:
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
      
    fitpot0 = "{}/{}/data/fitcharge.pot0".format(run_dir,fitfolder)
    with open(fitpot0,'w') as g:
       i = 'unperturbed'
       filename = "{}/fitcharge/configs_connolly/{}/{}_charges.vpot".format(run_dir,i,i) 
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
               if flag == 2:
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
    qpos = "{}/{}/data/fitcharge.qpos".format(run_dir,fitfolder)
    with open(qpos,'w') as g:
       for i in range(0,total_pc):
          for count_j,j in enumerate(geo):
              g.write('{:< 15.6f} {:< 15.6f} {:< 15.6f}\n'.format(j[0],j[1],j[2]))
          filename = "{}/fitcharge/configs/{}/pointcharges.pc".format(run_dir,i) 
          with open(filename,'r') as f:
              for lc,lines in enumerate(f):
                  fields = lines.split()

                  # Initialize lists/arrays
                  if lc == 1: 
                      g.write('{:< 15.6f} {:< 15.6f} {:< 15.6f}\n'.format(float(fields[1]),float(fields[2]),float(fields[3])))

    
    write_datadir(run_dir,fitfolder)
    return

def Drude_charge(run_dir,FinalFF,localFF,args,c,user,fitfolder):
    
    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+fitfolder)
    print("Start parametrization....") 
    if os.path.isfile('fitcharge.log') is False:     
       jobids = []
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('module load openmpi/3.1.4 \n/depot/bsavoie/apps/charmm/exec/gnu/charmm < fitcharge.inp > fitcharge.log')
       sublist= sublist + ('charmm_fit  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
       monitor_jobs(jobids,user)
    if os.path.isfile('fitcharge.log') is False:
       print("ERROR: fitcharge.log not found, something went wrong with executing CHARMM")
       quit()
    else:
       print("COMPLETE: CHARMM parametrization")

    # Converting log file to db file   
    command = ('-o {}'.format(args.outFF)).split()
    log_to_db.main(command)

    os.chdir(run_dir)
    # Merge these Params with the parent
    # XXX can't merge with parent, if the other inchi is also merging then would have error
    option = {}
    #option['master'] = FinalFF
    option['new_params'] = './{}/{}'.format(fitfolder,args.outFF)
    #option['no_ct'] = True
    #option['replace'] = 1
    # Merge these with the local (total_drude.db)
    option['master'] = localFF 
    option['no_ct'] = False
    merge_FF_drude.fun(option)
    print("COMPLETE: Drude charge parametrization")
   
    os.chdir(ori_dir)

    return

def VDW_re(run_dir,args,tag,inchikey,current_gen,xyz):

    ori_dir = os.getcwd()
    print("Start vdw reparametrization...")
    # move to non_polar folder to relax Drude particles
    vdwfolder = nonpolar_dir + '/' + inchikey + '/vdw' 
    if os.path.isdir(vdwfolder) is False:
      print("ERROR: nonpolar vdw folder {} doesn't exist. Exiting...".format(vdw))
    os.chdir(vdwfolder)
    print("move to {} to reparametrize...".format(vdwfolder))
    # converting the new db file to charmm format FF which is needed when relaxing via CHARMM
    FFname = 'FF{}'.format(tag)
    if os.path.isdir(vdwfolder+'/'+FFname) is False:
      #sublist = ('{}/geo_opt.xyz  -o FF -q {} -gens {} -mixing_rule wh --polar -FF '.format(run_dir,c['charge'],c['gens'])).split()
      sublist = ('{}  -o {}  -gens {} -mixing_rule wh --polar -FF '.format(xyz,FFname,c['gens'])).split()
      FFstring = '{} {} {} '.format(intraonly,chargeonly,VDWonly)
      # append lower gens' local Drude FF
      for i in range(1,current_gen):
         FFstring += '{}/gen{}.db '.format(driver_dir,i)
      FFstring += '{}/drude_Q.db'.format(run_dir)
      sublist.append(FFstring)
      #gen_md_for_CHARMM.main(sublist)
      db_to_charmm.main(sublist)
    os.chdir('{}/configs'.format(vdwfolder))
    if tag == '':
      relax_drude_loop.main([])
    else:
      relax_drude_loop.main('-tag {}'.format(tag).split())

    # initial_param.db serves as extract_vdw's initial guess
    # normally this would be the FF from the cycle before the very last cycle
    # we want to use the fully parametrized FF as our initial guess
    os.chdir(vdwfolder)
    shutil.copyfile('{}/cycle-20/DFT-AA.db'.format(vdwfolder),'{}/initial_params.db'.format(vdwfolder))
    option = {}
    option['master'] = '{}/initial_params.db'.format(vdwfolder)
    option['no_ct'] = False
    option['replace'] = 1
    for i in range(1,current_gen):
      option['new_params'] = '{}/gen{}.db'.format(driver_dir,i)
      merge_FF_drude.fun(option)
    
    # mv to inchikey folder to execute extract_vdw
    os.chdir('{}/{}'.format(nonpolar_dir,inchikey))

    L2_s = "0.1" 
    L2_e = "0.5"
    jobids = []
    vdwdrude_dir = 'vdw_drude{}'.format(tag) 
    if os.path.isdir(vdwfolder+'/'+vdwdrude_dir) is False:
       # input for existing data base is genN.db (Drude vdw only)
       # initial_parms serves as initial guess, a mix of taffi VDW and Drude VDW
       FF_path = ' '.join(['{}/gen{}.db'.format(driver_dir,i) for i in range(1,current_gen)]) 
       charge_FF = '{}/drude_Q.db'.format(run_dir) 
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       if tag == '':
         sublist.append('{}/FF_functions/extract_vdw_drude.py -f vdw -FF_DFT \'{} {}\' -o {} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} '.format(c['taffi_path'],FF_path,charge_FF,vdwdrude_dir,L2_s,L2_e,c['functional']))
       else:
         sublist.append('{}/FF_functions/extract_vdw_drude.py -f vdw -FF_DFT \'{} {}\' -o {} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} -tag {}'.format(c['taffi_path'],FF_path,charge_FF,vdwdrude_dir,L2_s,L2_e,c['functional'],tag))
         
       sublist= sublist + ('vdw_drude_parse  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
       # Wait until jobs complete
       monitor_jobs(jobids,user)
    if os.path.isfile(vdwfolder+'/'+vdwdrude_dir+'/DFT-AA.db') is False:
       print("ERROR: vdw/{}/DFT-AA.db not found, something went wrong w/ vdw_parse_drude".format(vdwdrude_dir))
       quit()
    else:
       print("COMPLETE: reparametrize vdw")

    os.chdir(ori_dir)

    return
    

def write_datadir(run_dir,fitfolder):
    with open(run_dir+'/'+fitfolder+'/datadir.def','w') as f:
       f.write("* CHARMM testcase data directory assignment\n"+\
               "*\n"+\
               "faster on\n"+\
               "set pnode =\n"+\
               "if ?numnode .gt. 1 set pnode = node 0\n"+\
               "set 0 data/     ! input data directory\n"+\
               "set 9 scratch/  ! scratch directory\n"+\
               "set testcheck stream @0/test.str\n"+\
               "set qcheck stream @0/qtest.str\n"+
               "set testfail 0\n"+\
               "return\n")
    return

def prepare(run_dir,name):
   i = name
   shutil.copyfile('{}/fitcharge/configs/{}/{}_charges.gbw'.format(run_dir,i,i),'{}_charges.gbw'.format(i))
   shutil.copyfile('{}/fitcharge/configs/{}/{}_charges.scfp'.format(run_dir,i,i),'{}_charges.scfp'.format(i))
   #shutil.copyfile('{}/fitcharge/configs/{}/{}.xyz'.format(run_dir,i,i),'{}.xyz'.format(i))
   #if i != 'unperturbed':
   #   shutil.copyfile('{}/fitcharge/configs/{}/pointcharges.pc'.format(run_dir,i),'pointcharges.pc')
   shutil.copyfile('{}/fitcharge/map.xyz'.format(run_dir),'map')

   return
   

# This function grep FF_path from vdw_parse.submit
# it can also change the original folder path
# this can be handy if the original folder has been moved 
# or when the parametrization happened in a different folder
# it also change total.db to total_drude.db
# This function is only needed right now because the generation structure hasn't been coded

def get_FF_path(filename,changefolder=''):
    with open(filename,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split() 
         if len(fields) < 5: continue
         if fields[3] == '-FF_DFT':
            newfields = lines.split('\'')
            FF = newfields[1].split()[1:]
            for count_i,i in enumerate(FF):
               singleFF = i.split('/')
               if singleFF[-1] == 'total.db':
                  singleFF[-1] = 'total_drude.db'
                  FF[count_i] = '/'.join(singleFF)
            FF = ' '.join(FF)
    return FF
            

def Write_charge_in(c,name,geo,element):
   with open('charge.in','w') as f:
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


def GetInchi(Elements,Geometry): #convert geo2inchi key

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('inchi.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    # create inchikey
    p = subprocess.Popen("obabel -ixyz inchi.xyz -oinchikey".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
    inchikey, err = p.communicate()
    inchikey = str(inchikey, 'utf-8').split('\n')[0] ##might be a problem if there are characters utf-8 can't handle
    # remove xyz file
    os.remove("inchi.xyz")
   
    return inchikey
"""
def Get_Inchinew(Elements,Geometry):

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('inchi.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "inchikey")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "inchi.xyz")

    inchikey=obConversion.WriteString(mol).rstrip()
    # remove xyz file
    os.remove("inchi.xyz")

    return inchikey

def Get_Inchi_xyz(xyz):

    import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "inchikey")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, xyz)

    inchikey=obConversion.WriteString(mol).rstrip()

    return inchikey
"""


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

def clean_pcsubmit(run_dir,nconf):
   for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = '{}/fitcharge_{}'.format(run_dir,i)
      delete_wildcard(fitfolder,'pc_layer')

   return

def delete_wildcard(folder,pattern,addlayer=False):
    if addlayer:
      fileList = glob.glob('{}/*/*{}*'.format(folder,pattern))
    else:
      fileList = glob.glob('{}/*{}*'.format(folder,pattern))
    for filePath in fileList:
       try:
          os.remove(filePath)
       except:
          print("Error while deleting file : ", filePath)

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
