#!/bin/env python                                                                                                                                                             
import sys,argparse,subprocess,os,time,math,shutil
from copy import deepcopy
from matplotlib import pyplot as plt
import numpy as np

def main(argv):

   #parser = argparse.ArgumentParser(description='post process ub jobs to calculate kTST')
    
    
   #parser.add_argument('ligand_name', type=str,
   #                     help='name of ligand')

   #args = parser.parse_args()

   with open('smiles_name','r') as f:
      for lc,lines in enumerate(f):
         fields  = lines.split()
         name = fields[1]
         smiles = fields[0]
         subprocess.call("obabel -:\"{}\" -oxyz -O {}.xyz --gen3d".format(smiles,name),shell=True)
   

               
# Create logger to save stdout to logfile
class Logger(object):

    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/exe.log", "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass
   

if __name__ == "__main__":
   main(sys.argv[1:])
