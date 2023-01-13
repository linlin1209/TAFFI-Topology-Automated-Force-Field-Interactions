#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in charmm crd file then transform it to xyz file') 


    #required (positional) arguments                                                                                                  
    parser.add_argument('crdfile', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')

    parser.add_argument('-o', dest='outfile', default='',
                        help = 'output xyz filename, default: name of xyz file')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')

    # Make parse inputs
    args=parser.parse_args(argv)
   

    if args.outfile == '':
      args.outfile = args.crdfile.split('.')[-2].split('/')[-1] + '.xyz'

    elements,geo = read_crd(args.crdfile)
    write_xyz(elements,geo,args.outfile)
         
    print("Successfully convert {} to {}".format(args.crdfile,args.outfile))
    return


def write_xyz(elements,geo,filename):
   with open(filename,'w') as f:
      f.write('{}\n\n'.format(len(elements)))
      for count_i,i in enumerate(elements):
         f.write("{:3s} {:< 6.5f} {:< 6.5f} {:< 6.5f}\n".format(i,geo[count_i][0],geo[count_i][1],geo[count_i][2]))
   return
      

def read_crd(filename):
   elements = []
   Geo = []
   with open(filename,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields)>0 and split_string(fields[0])[0] == '*': continue
         if len(fields)<2: continue
         elements.append(split_string(fields[3])[0]) 
         Geo.append([float(i) for i in fields[4:7]])
   return elements,Geo

def split_string(word):
    return [char for char in word]
         


if __name__ == "__main__":
   main(sys.argv[1:])
