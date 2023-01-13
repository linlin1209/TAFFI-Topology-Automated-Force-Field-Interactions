#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re

# NOTE: the order of the imports matters until I resolve the namespace issues (XXX DON'T USE import *)
import random
import numpy as np
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

def main(argv):

    parser = argparse.ArgumentParser(description= 'calculate average and std dev based for each column')


    #required (positional) arguments                                                                                                  
    parser.add_argument('infile', help = 'concentration file to read')

    # default is 1 because the title 0 is repeqated
    parser.add_argument('-start', dest='t_start', default=1,
                        help = 'start time for calculating average')

    parser.add_argument('-block_freq', dest='block', default=1,
                        help = 'block average frequence')
   
    # Make parse inputs
    args=parser.parse_args(argv)
    args.t_start = int(args.t_start)
    args.block = int(args.block)

    with open(args.infile,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if lc == 0: 
            title = fields
            value = {k:[] for k in title}
            continue
         if lc < args.t_start: continue
         for count_i,i in enumerate(fields):
            value[title[count_i]].append(float(i))
    value_new = {}
    for key in value:
      value_new[key] = [] 
      for i in range(0,len(value[key])//args.block):
         value_new[key].append(np.mean(np.array(value[key][i*args.block:(i+1)*args.block])))
      if len(value[key])%args.block != 0:
         value_new[key].append(np.mean(np.array(value[key][len(value[key])//args.block*args.block:])))
    name = args.infile.split('.')[-2]+'_block.'+args.infile.split('.')[-1]
    with open(name,'w') as f:
      for key in value_new:
         f.write('{:<15s} '.format(key))
      f.write('\n')
      for i in range(0,len(value_new[list(value_new.keys())[0]])):
         for key in value_new:
            f.write('{:< 15.4f} '.format(value_new[key][i]))
         f.write('\n')
    avg = []
    print_s  = ""
    print_mean = ""
    print_std = ""
    # gather print title string
    for key in value:
      print_s += "{:<10s} ".format(key) 
    print(print_s)
    # calculate avg and std and put it in printout string
    for key in value:
      mean = np.mean(np.array(value[key]))
      std = np.std(np.array(value[key]))
      avg.append(mean)
      print_mean += "{:<10.3f} ".format(mean)
      print_std += "{:<10.3f} ".format(std)
    print(print_mean)
    print(print_std)
      
         
    

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/gro_to_lmp.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

    def flush(self):
        pass


if __name__ == "__main__":
   main(sys.argv[1:])
