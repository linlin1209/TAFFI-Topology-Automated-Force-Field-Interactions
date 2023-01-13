#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *
import xyz_to_gro

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

    parser = argparse.ArgumentParser(description='take a number of molecules, use the input xyz file and density to output a initial structure in pdb file by gromacs') 
                                               


    #required (positional) arguments                                                                                                  
    parser.add_argument('nmol', help ='desired number of molecules')

    # optional arguments
    parser.add_argument('-xyz',dest='xyz',default='geo_opt.xyz',
                        help = 'xyz file which will be used to compute molecular weight, default: geo_opt.xyz')
    parser.add_argument('-molw',dest='molw',default=0,
                        help = 'molecular weight that will be used to compute volume, this is only needed for calc_vol mode, default: 0')
    parser.add_argument('-dens',dest='dens',default=750,
                        help = 'desired density in kg/m^3 , default: 750 kg/m^3')
    parser.add_argument('-Vmol',dest='Vmol',default=200,
                        help = 'desired molecular volume  in A^3 , default: 200 A^3')
    parser.add_argument('--calc_only', dest='calc_only', default=False, action='store_const', const=True,
                        help = 'When present, only output the supposed volume, but not generate the box, when this is on, -molw need to be specified')
    parser.add_argument('--use_Vmol', dest='use_Vmol', default=False, action='store_const', const=True,
                        help = 'When present, will use Vmol instead of density')
                                         


    # Make parse inputs
    args=parser.parse_args(argv)
    
    args.dens = float(args.dens)
    args.Vmol = float(args.Vmol)
    args.molw = float(args.molw)
    args.nmol = int(args.nmol)
    xyz = args.xyz

    if os.path.isfile(xyz) is False and args.molw == 0:
         print("ERROR: molecular weight has to be specified directly by -molw or calculated from  -xyz ! Exiting...")
         quit()
    if os.path.isfile(xyz):
         molw = get_molw(xyz)
    if args.molw != 0:
         molw = args.molw
   
    #       kg->g   g-->mol    mol -> #      m^3->A^3
    Const = 1000  * 1/molw   *  (6.02*1e23) *  1e-30    # unit: #/A^3
    
    if args.use_Vmol :
         vol = args.Vmol*args.nmol
    else:
         # convert density from kg/m^3 to #/A^3
         dens = args.dens*Const
         vol = args.nmol / dens
    boxl = vol**(1/3)
   
    if args.calc_only:
         print("{:} molecules in {:3.2f} x {:3.2f} x {:3.2f} A^3 box will have density of {:3.2f} kg/m^3, molecular volume of {:3.2f} A^3".format(args.nmol,boxl,boxl,boxl,args.dens,vol/args.nmol))
         quit()
    sublist = ('{} -o single.gro'.format(xyz)).split()
    xyz_to_gro.main(sublist)
    title = 'single'
    print("putting {:} molecules into {:3.2f} x {:3.2f} x {:3.2f} A^3 box".format(args.nmol,boxl+5,boxl+5,boxl+5))
    #increase the box for 2 A in case too dense for Grommacs to put all the mols
    command = "gmx_mpi insert-molecules -ci {}.gro -nmol {} -box {} {} {} -o  initial.pdb".format(title,args.nmol,(boxl/10)+0.5,(boxl/10)+0.5,(boxl/10)+0.5) # /10 cuz gromacs in nm
    print(command)
    output = subprocess.Popen(command.split(),stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()[-1]
    #print(output.split())

    return boxl+5 
    


    

def get_molw(filename):

    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}
    
    with open(filename,'r') as f:
         flag = 0
         atom_count = 0
         molw = 0
         for lc,lines in enumerate(f):
            fields = lines.split()  
            if lc == 0: 
               natoms = int(fields[0])
               continue
            if lc == 1: 
               flag =1   
               continue
            if flag == 1:
               if atom_count == natoms: break
               if fields[0] not in list(mass_dict.keys()):
                  print("Unknowm atom: {} in the {} file, Exiting...".format(fields[0],filename))
                  quit()
               molw += mass_dict[fields[0]]
               atom_count += 1
    print("{:} has molecular weight of {:3.5f} g/mol".format(filename,molw))

    return molw
               


if __name__ == "__main__":
   main(sys.argv[1:])
