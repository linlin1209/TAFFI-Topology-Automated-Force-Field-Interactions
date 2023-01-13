#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,os,re


def main(argv):

   old_prm = 'total.db'
   new_prm = 'total_new.db'
   bond_dict = {}
   angle_dict = {}
   tors_dict = {}
   vdw_dict = {}
   bond_dict_old = {}
   angle_dict_old = {}
   tors_dict_old = {}
   vdw_dict_old ={}
   with open(new_prm,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) < 2 :continue
         if fields[0] == 'bond': 
            bond_name = '_'.join(fields[1:3])
            bond_dict[bond_name] = [float(fields[4]),float(fields[5])]
         if fields[0] == 'angle':
            angle_name = '_'.join(fields[1:4])
            angle_dict[angle_name] = [float(fields[5]),float(fields[6])]
         if fields[0] == 'torsion':
            tors_name = '_'.join(fields[1:5])
            if fields[5] == 'opls':
               tors_dict[tors_name] = [float(fields[6]),float(fields[7]),float(fields[8]),float(fields[9])]
            else:
               tors_dict[tors_name] = [float(fields[6]),float(fields[7]),float(fields[8])]
         if fields[0] == 'vdw' and fields[1] == fields[2]:
            type_name = fields[1]
            vdw_dict[type_name] = [float(fields[4]),float(fields[5])]
   with open(old_prm,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) < 2 :continue
         if fields[0] == 'bond': 
            bond_name = '_'.join(fields[1:3])
            bond_dict_old[bond_name] = [float(fields[4]),float(fields[5])]
         if fields[0] == 'angle':
            angle_name = '_'.join(fields[1:4])
            angle_dict_old[angle_name] = [float(fields[5]),float(fields[6])]
         if fields[0] == 'torsion':
            tors_name = '_'.join(fields[1:5])
            if fields[5] == 'opls':
               tors_dict_old[tors_name] = [float(fields[6]),float(fields[7]),float(fields[8]),float(fields[9])]
            else:
               tors_dict_old[tors_name] = [float(fields[6]),float(fields[7]),float(fields[8])]
         if fields[0] == 'vdw' and fields[1] == fields[2]:
            type_name = fields[1]
            vdw_dict_old[type_name] = [float(fields[4]),float(fields[5])]
   for bond in list(bond_dict.keys()):
      dif = [bond_dict[bond][count_i] - bond_dict_old[bond][count_i] for count_i,i in enumerate(bond_dict[bond])]
      if min(dif) > 0.03:
         val = ' '.join(['{:< 3.5f}'.format(i) for i in bond_dict[bond]])
         old_val = ' '.join(['{:< 3.5f}'.format(i) for i in bond_dict_old[bond]])
         print("{} {} {}".format(bond,val,old_val))
   for angle in list(angle_dict.keys()):
      dif = [angle_dict[angle][count_i] - angle_dict_old[angle][count_i] for count_i,i in enumerate(angle_dict[angle])]
      if min(dif) > 0.03:
         val = ' '.join(['{:< 3.5f}'.format(i) for i in angle_dict[angle]])
         old_val = ' '.join(['{:< 3.5f}'.format(i) for i in angle_dict_old[angle]])
         print("{} {} {}".format(angle,val,old_val))
   for tors in list(tors_dict.keys()):
      dif = [tors_dict[tors][count_i] - tors_dict_old[tors][count_i] for count_i,i in enumerate(tors_dict[tors])]
      if min(dif) > 0.03:
         val = ' '.join(['{:< 3.5f}'.format(i) for i in tors_dict[tors]])
         old_val = ' '.join(['{:< 3.5f}'.format(i) for i in tors_dict_old[tors]])
         print("{} {} {}".format(tors,val,old_val))
   for vdw in list(vdw_dict.keys()):
      dif = [vdw_dict[vdw][count_i] - vdw_dict_old[vdw][count_i] for count_i,i in enumerate(vdw_dict[vdw])]
      if min(dif) > 0 or min(dif) <0:
         val = ' '.join(['{:< 3.5f}'.format(i) for i in vdw_dict[vdw]])
         old_val = ' '.join(['{:< 3.5f}'.format(i) for i in vdw_dict_old[vdw]])
         print("{} {} {}".format(vdw,val,old_val))
      
               


if __name__ == "__main__":
   main(sys.argv[1:])
