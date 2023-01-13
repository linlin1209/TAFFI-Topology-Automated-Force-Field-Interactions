#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,fnmatch,getpass
import numpy as np
from matplotlib import pyplot as plt
# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
import random
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
import db_to_charmm,xyz_to_pdb,density_converter,add_time_pdb,read_charmm_energy,gen_itp_newdrude,crd_to_xyz
import add_drude_pdb
from monitor_jobs import *

def main(argv):


    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_file', help = 'A string with the name of the  *.xyz file to be geometry optimized.')

    parser.add_argument('-o', dest='ofolder', default='minimize',
                        help = 'Sets the output filename prefix (default: minimize)')

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'generations for taffi atom type')

    parser.add_argument('-FF', dest='FF', default='',
                        help = 'force field for minimization')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')

    parser.add_argument('--shake', dest='shake_flag', default=False, action='store_const', const=True,
                        help = 'When present, will use shake to contrain h bonds (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)

    args.coord_file = args.coord_file.split()
    args.gens = int(args.gens)

    # Check that the input is an .xyz file. 
    if len(args.coord_file) > 1:
        print("ERROR in charmm_minimize_molecule: Only one geometry can be supplied. Exiting...")
        quit()
    elif args.coord_file[0].split('.')[-1] != 'xyz':
        print("ERROR in charmm_minimize_molecule: Check to ensure that the input coordinate file is in .xyz format. Exiting...")
        quit()
    # Make sure xyz is in absolute path
    if split_word(args.coord_file[0])[0] != '/':
      args.coord_file[0] = '{}/{}'.format(os.getcwd(),args.coord_file[0])

    # check xyz exist
    if os.path.isfile(args.coord_file[0]) is False:
        print("ERROR: xyz file: {} not found, Exiting...".format(args.coord_file[0]))
        quit()

    # Make the output directory if it doesn't already exist.
    if os.path.exists(args.ofolder):
        print('ERROR: The desired output folder ({}) already exists. Exiting to avoid overwriting data...'.format(args.ofolder))
        return
    else:
        os.makedirs(args.ofolder)
        sys.stdout = Logger(args.ofolder)
        print("PROGRAM CALL: python {}\n".format(' '.join([ i for i in sys.argv])))

    charmm_minimize_molecule(args.coord_file[0],args.gens,args.polar_flag,args.ofolder,args.FF)


def charmm_minimize_molecule(xyz,gens,polar_flag,ofolder,FF):

   os.chdir(ofolder)

   sublist = ('{} -o single.pdb'.format(xyz)).split()
   xyz_to_pdb.main(sublist)
   
   if polar_flag:
      sublist = ('{}   -o ff  -gens {} -mixing_rule wh --polar -FF '.format(xyz,gens)).split()
      sublist.append('{}'.format(FF))
      db_to_charmm.main(sublist)
   else:
      sublist = ('{}   -o ff  -gens {} -mixing_rule wh  -FF '.format(xyz,gens)).split()
      sublist.append('{}'.format(FF))
      db_to_charmm.main(sublist)


   # this inp has been tested with tip3p, see /depot/bsavoie/apps/charmm/test/tip3p_test/polar
   # the result match with /depot/bsavoie/apps/charmm/test/tip3p_test/
   with open('min.inp','w') as f: 
      f.write("* setup structures\n")
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n")
      f.write("read rtf card name ff/lin.rtf\n")
      f.write("read para card name ff/lin.prm\n")
      f.write("open unit 91 card read name single.pdb\n")
      f.write("read sequ pdb unit 91\n")
      f.write("set residue LIN\n\n")

      f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n\n")

      f.write("rewind unit 91\n")
      f.write("read coor pdb unit 91\n")
      f.write("coor sdrude\n")
      f.write("coor shake\n")
      f.write("coor print\n\n")


      f.write("set cutoff 99\n")
      f.write("set ctofnb = @cutoff\n")
      f.write("Calc ctofnb = @ctofnb - 2.0\n")
      f.write("set ctonnb = @ctofnb\n")
      f.write("Calc ctonnb = @ctonnb - 1.0\n\n")

      f.write("! shake constraint\n")
      f.write("SHAKE bonh param nofast -\n")
      f.write("select .not. type D* end -\n")
      f.write("select .not. type D* end\n\n")

      f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("MINI SD nstep 5000 nprint 50 -\n")
      f.write("atom vatom vdistance cdie eps 1.0 vswitch fshift wmin 1.2 -\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("e14fac 0.0\n\n") 
   
      f.write('coor orient\n')
      f.write('coor print\n')
      f.write("! First with DSCF to eliminate too high freq from Drude\n")
      f.write("! this should be a more correct setting, but the result seems problematic\n")
      f.write("!VIBRAN\n")
      f.write("!DIAGONALIZE DSCF Finite\n")
      f.write("!END\n\n")
      f.write("! calculate molecular polarizability\n")
      f.write("! This won't work with Diagonalize option on\n")
      f.write("! This is from example in documentation\n")
      f.write("shake off\n")
      f.write("VIBRAN\n")
      f.write("DIAGONALIZE\n")
      f.write("PRINT NORMAL VECTORS DIPOLES SELECT ALL END\n")
      f.write("END\n\n")

      f.write("! calculate dipole\n")
      f.write("COOR DIPOLE\n")
      f.write("SET MU1 = ?RDIP\n\n")

      f.write("OPEN WRITE CARD UNIT 33 NAME mu.log\n")
      f.write("WRITE TITLE UNIT 33\n")
      f.write("* MONOMER:\n")
      f.write("*   MU    @MU1\n")
      f.write("*\n")


      f.write("OPEN WRITE CARD UNIT 10 NAME mm_min.crd\n")
      f.write("WRITE COOR CARD UNIT 10\n")

   shell_command('charmm < min.inp > min.log \n wait','exe')
   sublist = ('mm_min.crd -o minimized.xyz'.split())
   crd_to_xyz.main(sublist)

   return

def split_word(word):
    return [char for char in word]


def shell_command(command,shname):
    with open('{}.sh'.format(shname),'w') as f:
         f.write("#!/bin/bash\n")
         f.write(command)
    command = 'sh {}.sh'.format(shname).split()
    process = subprocess.Popen(command)
    status = process.poll()                                                                                                                                                                        
    while status == None:
       status = process.poll()
       time.sleep(3)

    return

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

if __name__ == "__main__":
    main(sys.argv[1:])
