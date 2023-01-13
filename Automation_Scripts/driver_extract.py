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
import restart_scans
import codecs,json


def main(argv):
    


    user = getpass.getuser()
    run_dir = os.getcwd()

    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-o',dest='output',default='taffi.log',
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
    
    #print(c)
    frag_out = read_json('./out_inter.json')
    frag_out_gen = [frag_out[i]["gen"] for i in frag_out]
    
    # alljson has all frag_inter info
    alljson = read_alljson(run_dir+'/all.json')


    
    # loop over generation(major loop)
    for gen in range(1,max(frag_out_gen)+1):
      
      # prepare input dict for frag_gen
      geo_dict = {}
      for key in frag_out:
         if( frag_out[key]["gen"] == gen ):
            geo_dict[key] = frag_out[key]
     ######################################
     # VDW                                #
     ######################################
      process = {}
      L2_s ="0.1" 
      L2_e ="1.0"
      for key in geo_dict:
         os.chdir(key)
         sublist = ['{}/Automation_Scripts/shell_submit_extract.sh'.format(c['taffi_path'])]
         sublist.append('{}/FF_functions/extract_vdw_lst.py -f vdw -FF_DFT \'charges/CHELPG_calcs/charges/fit_charges.db {}\' -o outlier -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -method \'global Boltzmann Boltzmann_tot\' --remove_outlier'.format(c['taffi_path'],c['ff'],L2_s,L2_e))
         sublist= sublist + ('vdw_lst  -p 1 -t 4 -q {}'.format(c["vdw_qc_q"])).split()
         process[key] = subprocess.Popen(sublist)
         os.chdir(run_dir)
      # Update FF path
      for key in geo_dict:
         c['ff'] = c['ff'] +' ' +run_dir+'/'+key +'/total.db'

         
         
    quit()
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
