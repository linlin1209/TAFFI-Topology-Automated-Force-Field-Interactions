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
import numpy as np
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in CHARMM output pdb file and transform it to Gromacs compatible pdb which has time and volume ')

    #required (positional) arguments                                                                                                  

    parser.add_argument('-log', dest='logfile', default='md.log',
                        help = 'can also take wildcards like md.*.log for liquid simulation')

    parser.add_argument('-o', dest='outfile', default='energy.txt',
                        help = 'output log filename, default: energy.txt')

    parser.add_argument('-R', dest='recursive_opt', default=False, action='store_const',const=True,
                        help = 'When this flag is present, filenames and wildcards will be searched for recursively.')

    parser.add_argument('-N_start', dest='N_start', default=2,
                        help = 'start from md.2.log')

    parser.add_argument('--silent', dest='silent_flag', default=False, action='store_const',const=True,
                        help = 'When this flag is present, filenames and wildcards will be searched for recursively.')

    parser.add_argument('--gas', dest='gas_flag', default=False, action='store_const',const=True,
                        help = 'When this flag is present, volume will not be written')

    # Make parse inputs
    args=parser.parse_args(argv)
    
    args.N_start = int(args.N_start)

    # Find files
    args.logfile = args.logfile.split()
    wild_card_files  = [ i for i in args.logfile if "*" in i ]
    if args.recursive_opt == True:
        args.logfile = [ i for i in args.logfile if "*" not in i ]
        args.logfile += [ dp+"/"+files for i in wild_card_files for dp,dn,fn in os.walk('.') for files in fn if fnmatch.fnmatch(files,i) ] 
    else:
        args.logfile = [ i for i in args.logfile if "*" not in i ]
        for i in args.logfile:
            if os.path.isfile(i) is False:
               print("ERROR: {} file not found".format(args.logfile))
               quit()
        for i in wild_card_files:
            path = '/'.join(i.split('/')[:-1])
            if len(path) == 0:
                path = "."
            #args.logfile += [ path+"/"+files for files in os.listdir(path) if fnmatch.fnmatch(files,i) ] 
            args.logfile += [ path+"/"+files for files in os.listdir(path) if fnmatch.fnmatch(files,i.split('/')[-1]) ] 

    # Handle "./" condition
    for count_i,i in enumerate(args.logfile):
        if i[:2] == "./": args.logfile[count_i] = args.logfile[count_i][2:]    
    args.logfile = list(set(args.logfile))
    if len(args.logfile) > 1:
      args.logfile=sorted(args.logfile, key = log_num)   
      args.logfile = [ i for i in args.logfile if int(i.split('.')[-2]) >= args.N_start ]

    print("Processing: {}".format('  '.join(args.logfile)))

    
    print("START reading energy file {}".format(args.logfile))
    data = read_energy(args.logfile,args.silent_flag)
    print("Writing energy to {}".format(args.outfile))
    write_energy(data,args.outfile,args.gas_flag)
    #analysis(data)
    
    return

def log_num(x):
    return(int(x.split('.')[-2]))

def analysis(data):

    # These are identical to one in FOM_driver.py
    kb=0.001987204118                           # Boltzman constant in units of kcal/mol K
    kcal_mol_to_J = 4.184 * 1000.0 / 6.0221409E23 # employed for converting the isothermal compressibility to useful units.  
    kcal_to_kJ = 4.184
    pressure = 1.0                    # for thermal expansion coeff (atm)
    t = 298.15
    atm_to_pa = 101325.0              # atmospheres to pascals
    PV_const = pressure * 1E-30 * atm_to_pa / kcal_mol_to_J # A^-3 to m^-3, conversion to pa, conversion to kcal/mol
    Data = {}

    Data["U"] = [ data[key]['E_pot'] for key in data ] 
    Data["V"] = [ data[key]['Volume'] for key in data ] 

    # Convert to arrays
    Data["U"] = np.array(Data["U"])
    Data["V"] = np.array(Data["V"])

    # Calculate isothermal compressibility
    Data["iso_comp"] = np.std(Data["V"])**(2.0)/(np.mean(Data["V"])*kb*float(t))*1E-30*1E9/kcal_mol_to_J 

    # # Calculate thermal expansion coefficient (OLD)
    # tmp_dev_V = Data["V"] - mean(Data["V"])
    # tmp_dev_H = Data["H"] - mean(Data["H"])
    # Data["therm_exp"] = mean(tmp_dev_V*tmp_dev_H)/(kb*float(t)**(2.0)*mean(Data["V"]))*1E3

    # Calculate thermal expansion coefficient
    Data["therm_exp"] = ( np.mean(Data["V"]*Data["U"]) - np.mean(Data["V"])*np.mean(Data["U"]) + PV_const * np.std(Data["V"])**(2.0) ) \
                             / (kb*float(t)**(2.0)*np.mean(Data["V"])) * 1E3
    print("isothermal compressibility: {}".format(Data["iso_comp"]))
    print("thermal expansion coeff: {}".format(Data["therm_exp"]))



def write_energy(data,output,gas_flag=False):

    with open(output,'w') as f:
      f.write("{:<10s} {:>20s} {:>20s} {:>20s} {:>20s} {:>20s} {:>20s} {:>20s} {:>20s}\n".format('time','E_tot','E_k','E_pot','HFCTote','EHFCor','Volume','Custom1','Custom2'))
      for count_i,time in enumerate(data):
         if count_i == 0: 
            ener_key = list(data[time].keys())
            mean = {}
            for i in ener_key: mean[i] = []
            #mean['custom1'] = []
            #mean['custom2'] = []
            #mean['custom3'] = []
         for i in ener_key:
            if i not in list(data[time].keys()):
               data[time][i] = 0 # this happens sometime in gas simulation, where sometimes intermolecular are missing
            mean[i] += [data[time][i]]
         #custom1 = data[time]['EXTERN']+data[time]['INTERN'] 
         #custom2 = data[time]['E_k']+data[time]['E_pot']
         #custom3 = data[time]['HFCTote']+data[time]['EHFCor']
         #mean['custom1'] += [custom1]
         #mean['custom2'] += [custom2]
         #mean['custom3'] += [custom3]
         #f.write("{:>10.2f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f}\n".format(time,data[time]['E_tot'],data[time]['E_k'],data[time]['E_pot'],data[time]['HFCTote'],data[time]['EHFCor'],custom1,custom2,custom3))
         if gas_flag:
            f.write("{:>10.2f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} \n".format(time,data[time]['E_tot'],data[time]['E_k'],data[time]['E_pot'],data[time]['HFCTote'],data[time]['EHFCor']))
         else:
            f.write("{:>10.2f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f}\n".format(time,data[time]['E_tot'],data[time]['E_k'],data[time]['E_pot'],data[time]['HFCTote'],data[time]['EHFCor'],data[time]['Volume']))
      for i in mean:
         mean[i] = np.mean(np.array(mean[i]))
      #f.write("{:>10s} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f}\n".format('avg',mean['E_tot'],mean['E_k'],mean['E_pot'],mean['HFCTote'],mean['EHFCor'],mean['custom1'],mean['custom2'],mean['custom3']))
      if gas_flag:
         f.write("{:>10s} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} \n".format('avg',mean['E_tot'],mean['E_k'],mean['E_pot'],mean['HFCTote'],mean['EHFCor']))
      else:
         f.write("{:>10s} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f} {:>20.4f}\n".format('avg',mean['E_tot'],mean['E_k'],mean['E_pot'],mean['HFCTote'],mean['EHFCor'],mean['Volume']))

    return

def read_energy(logfile,silent_flag=False):

   data = {}
   for filename in logfile:
      time = -1
      if silent_flag is False:
         print("reading {}".format(filename))
      with open(filename,'r') as f:
         flag = 0
         for lines in f:
            fields = lines.split()
            if len(fields) < 3: continue
            if flag == 0 and fields[0] == 'DYNA>':
               time_pre = time
               time = float(fields[2])
               #if time_pre == time:
               #   print("ERROR: duplicate time at {:10.2f}, check file {:}".format(time,filename))
               #   quit()
               data[time] = {}
               data[time]['E_tot'] = float(fields[3])
               data[time]['E_k'] = float(fields[4])
               data[time]['E_pot'] = float(fields[5])
               data[time]['Temp'] = float(fields[6])
               if silent_flag is False:
                  print("Working on time: {:>10.2f}".format(time))
               flag = 1
               continue
            if flag == 1 and fields[0] == 'DYNA' and fields[1] == 'PROP>':
               # HFCTote: high frequency total energy
               data[time]['HFCTote'] = float(fields[3])
               data[time]['HFCKe'] = float(fields[4])
               # high frequency energy correction
               data[time]['EHFCor'] = float(fields[5])
               continue
            if flag == 1 and fields[0] == 'DYNA' and fields[1] == 'INTERN>':
               data[time]['bonds'] = float(fields[2])
               data[time]['angles'] = float(fields[3])
               data[time]['UREY-b'] = float(fields[4])
               data[time]['dihedrals'] = float(fields[5])
               data[time]['impropers'] = float(fields[6])
               data[time]['INTERN'] = np.sum(np.array([ float(i) for i in fields[2:7] ]))
               continue
            if flag ==1 and fields[0] == 'DYNA' and fields[1] == 'EXTERN>':
               data[time]['vdw'] = float(fields[2])
               data[time]['elec'] = float(fields[3])
               data[time]['EXTERN'] = float(fields[2])+float(fields[3])
               continue
            if flag == 3 and fields[0] == 'DYNA' and fields[1] == 'EWALD>':
               data[time]['EWKSum'] = float(fields[2])
               data[time]['EWSelf'] = float(fields[3])
               data[time]['EWEXcl'] = float(fields[4])
               data[time]['EWQCor'] = float(fields[5])
               data[time]['EWUTil'] = float(fields[6])
               continue
            if flag ==1 and fields[0] == 'DYNA' and fields[1] == 'PRESS>':
               data[time]['VIRE'] = float(fields[2])
               data[time]['VIRI'] = float(fields[3])
               data[time]['PRESSE'] = float(fields[4])
               data[time]['PRESSI'] = float(fields[5])
               data[time]['Volume'] = float(fields[6])
               flag = 0
               continue
   print("TOTAL {} frames read".format(len(data.keys())))
   return data
   


if __name__ == "__main__":
   main(sys.argv[1:])
