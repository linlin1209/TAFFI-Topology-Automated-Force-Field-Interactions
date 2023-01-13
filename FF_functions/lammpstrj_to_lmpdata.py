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
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in lammps data file then reorganize the atom index in case of nonconsecutive')


    #required (positional) arguments                                                                                                  
    parser.add_argument('trj_file', help = 'lammps trj file')

    parser.add_argument('data_file',help = 'lammps data file')
                                            
    # optional arguments
    #parser.add_argument('-del_atom', dest='del_atom', default=0,
    #                    help = 'delete atom index if it was not deleted in the input file')


    # Make parse inputs
    args=parser.parse_args(argv)
   
    fun(args,argv)
    
    quit()


def fun(args,argv):

    
    # Check that the input is an .gro file. 
    if args.data_file.split('.')[-1] != 'data':
          print("ERROR: Check to ensure that the input coordinate file(s) are in .data format.")
          quit()
    elif os.path.isfile(args.data_file) != True:
          print("ERROR: Specified *.data file ({}) does not exit. Please check the path. Exiting...".format(args.data_file))
          quit()
      
    print("PROGRAM CALL: lmpdata_reindex.py {}\n".format(' '.join([ i for i in argv])))

    # Initialize mass_dict (used as default values in several places).
    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}
    with open(args.trj_file,'r') as f:
      flag = 0
      for lc,lines in enumerate(f):
         fields = lines.split()
         if flag == 0 and len(fields) < 3: continue
         if flag == 1 and fields[0] == 'ITEM:':
            break
         if fields[0] == 'ITEM:' and fields[1] == 'ATOMS' and fields[2] == 'id':
            flag  = 1
            coord = {}
            vel = {}
            continue
         if fields[0] == 'ITEM:' and fields[1] == 'BOX' and fields[2] == 'BOUNDS':
            flag  = 2
            box_lo  = []
            box_hi = []
            continue
         if flag == 1:
            index = int(fields[0])
            coord[index] = [ i for i in fields[2:5]]
            vel[index] = [ i for i in fields[5:8]]
            continue
         if flag == 2:
            if len(box_lo) == 3: 
               flag = 0
               continue
            box_lo.append(fields[0])
            box_hi.append(fields[1])
            continue

         
    Filename = 'data_from_lammpstrj.data'
    g = open(Filename,'w')
    with open(args.data_file,'r') as f:
      flag = 0
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) == 4 and fields[2] == 'xlo':
            fields[0] = box_lo[0]
            fields[1] = box_hi[0]
            g.write('{}\n'.format(' '.join(fields)))
            continue
         if len(fields) == 4 and fields[2] == 'ylo':
            fields[0] = box_lo[1]
            fields[1] = box_hi[1]
            g.write('{}\n'.format(' '.join(fields)))
            continue
         if len(fields) == 4 and fields[2] == 'zlo':
            fields[0] = box_lo[2]
            fields[1] = box_hi[2]
            g.write('{}\n'.format(' '.join(fields)))
            continue
         if flag == 0: 
            if len(fields)>0 and fields[0] == 'Velocities':
               flag = 2
               tag = 0
               g.write(lines)
               g.write('\n')
               continue
            else:
               g.write(lines)
            
         if len(fields)>0 and fields[0] == 'Atoms':
            flag = 1
            g.write('\n')
            tag = 0
            continue
         if flag == 1 :
            if len(fields) == 0: tag += 1
            else:
               index = int(fields[0])
               for count_i,i in enumerate(range(4,7)):
                  fields[i] = coord[index][count_i]
               g.write('{}\n'.format(' '.join(fields)))
               continue
            if tag == 2: 
               flag = 0
               g.write('\n')
            continue
         if flag == 2 :
            if len(fields) == 0: tag += 1
            else:
               index = int(fields[0])
               for count_i,i in enumerate(range(1,4)):
                  fields[i] = vel[index][count_i]
               g.write('{}\n'.format(' '.join(fields)))
               continue
            if tag == 2: 
               flag = 0
               g.write('\n')
            continue
    #print("Successfully convert {} to {}".format(args.data_file,Filename))

    return
# takes in a sorted array and find the missing number
def missing_number(a):
    return sum(range(a[0],a[-1]+1)) - sum(a)

def read_gro(filename):
    atom = {}
    id_atom = {} # key:mol_id value:atoms
    nmol = {} # key: resi name value: nmol for that residue
    with open(filename,'r') as f:
       for lc,lines in enumerate(f):
         fields = lines.split()
         if lc == 0 : continue
         if lc == 1: 
            tot_atom = int(fields[0])
            continue
         if lc == tot_atom+2:
            box = float(fields[0])*10
            continue
         info = read_gro_format(lines)
         atom_count = info['atom_count']-1
         atom[atom_count] = {}
         atom[atom_count]['mol_id'] = info['mol_id']-1
         atom[atom_count]['resi'] = info['resname']
         atom[atom_count]['coord'] = np.array([info['x']*10,info['y']*10,info['z']*10])
         atom[atom_count]['elem'] = info['element']
         atom[atom_count]['atom_name'] = info['atom_name']
         if atom[atom_count]['mol_id'] not in list(id_atom.keys()): 
            id_atom[atom[atom_count]['mol_id']] = [] 
         if info['resname'] not in nmol.keys():
            nmol[info['resname']] = 0
         if atom_count == 0 or atom[atom_count]['mol_id'] != atom[atom_count-1]['mol_id']:
            nmol[info['resname']] += 1
         id_atom[atom[atom_count]['mol_id']].append(atom_count)
    return atom,id_atom,nmol,box

def read_gro_format(line):
   info = {}
   info['mol_id'] = int(line[0:5]) 
   info['resname'] = str(''.join([ i for i in line[5:10] if i != ' ']))
   info['atom_name'] = str(''.join([ i for i in line[10:15] if i != ' ']))
   info['element'] = split_string(info['atom_name'])[0]
   info['atom_count'] = int(line[15:20])
   info['x'] = float(line[20:28])
   info['y'] = float(line[28:36])
   info['z'] = float(line[36:44])

   return info

def Write_data_sphere(Filename,atom,id_atom,nmol,box):
    # Write the data file
    #with open(Filename+'/'+Filename.split('/')[-1]+'.data','w') as f:
    with open(Filename+'/in.data','w') as f:
        
        # Write system properties
        f.write("LAMMPS data file via vdw_self_gen.py, on {}\n\n".format(datetime.datetime.now()))

        f.write("{} atoms\n".format(len(list(atom.keys()))))
        f.write("{} atom types\n".format(len(nmol.keys())))
        f.write("\n")

        # Write box dimensions
        f.write("{:< 20.16f} {:< 20.16f} xlo xhi\n".format(-box/2,box/2))
        f.write("{:< 20.16f} {:< 20.16f} ylo yhi\n".format(-box/2,box/2))
        f.write("{:< 20.16f} {:< 20.16f} zlo zhi\n\n".format(-box/2,box/2))

        # Write Masses
        mass = {'HYD':1.00794*2,'NIT':28.01348,'AMM':14.00674+1.00794*3}
        mass = {'NIO':30.01,'NIT':28.01348,'OXY':32.0}
        resi = {'HYD':1,'NIT':2,'AMM':3}
        resi = {'NIO':1,'NIT':2,'OXY':3}
        f.write("Masses\n\n")
      
        for i in resi.keys():
            f.write("{:} {:< 8.6f}\n".format(resi[i],mass[i])) # count_i+1 bc of LAMMPS 1-indexing
         
        #for count_i,i in enumerate(list(nmol.keys())):
        #    f.write("{} {:< 8.6f}\n".format(resi[i],mass[i])) # count_i+1 bc of LAMMPS 1-indexing
        #    #resi[i] = count_i+1 
        f.write("\n")

        # Write Atoms
        f.write("Atoms\n\n")
        for i in atom:
            f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
            .format(i+1,i,resi[atom[i]['resi']],0,atom[i]['coord'][0],atom[i]['coord'][1],atom[i]['coord'][2],0,0,0))

    return

def Write_data_simple(Filename,atom,id_atom,nmol,box):
    # Write the data file
    #with open(Filename+'/'+Filename.split('/')[-1]+'.data','w') as f:
    with open(Filename+'/in.data','w') as f:
        
        # Write system properties
        f.write("LAMMPS data file via vdw_self_gen.py, on {}\n\n".format(datetime.datetime.now()))

        f.write("{} atoms\n".format(len(list(atom.keys()))))
        f.write("{} atom types\n".format(len(nmol.keys())))
        f.write("\n")

        # Write box dimensions
        f.write("{:< 20.16f} {:< 20.16f} xlo xhi\n".format(-box/2,box/2))
        f.write("{:< 20.16f} {:< 20.16f} ylo yhi\n".format(-box/2,box/2))
        f.write("{:< 20.16f} {:< 20.16f} zlo zhi\n\n".format(-box/2,box/2))

        # Write Masses
        mass = {'HYD':1.00794*2,'NIT':28.01348,'AMM':14.00674+1.00794*3}
        resi = {'HYD':1,'NIT':2,'AMM':3}
        f.write("Masses\n\n")
      
        #for i in resi.keys():
        #    f.write("{} {:< 8.6f}\n".format(resi[i],mass[i])) # count_i+1 bc of LAMMPS 1-indexing
         
        #for count_i,i in enumerate(list(nmol.keys())):
        #    f.write("{} {:< 8.6f}\n".format(resi[i],mass[i])) # count_i+1 bc of LAMMPS 1-indexing
        #    #resi[i] = count_i+1 
        f.write("\n")

        # Write Atoms
        f.write("Atoms\n\n")
        atom_count = 1
        for i in atom:
            if atom[i]['resi'] == 'NTG':
               charges   = {1:-0.40505,2:-0.40505,3:0.8101}
               mol_id = atom[i]['mol_id']
               if i%3 == 2:
                  f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
                  .format(atom_count,mol_id,2,charges[i%3+1],atom[i]['coord'][0],atom[i]['coord'][1],atom[i]['coord'][2],0,0,0))
                  atom_count += 1
               else:
                  f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
                  .format(atom_count,mol_id,1,charges[i%3+1],atom[i]['coord'][0],atom[i]['coord'][1],atom[i]['coord'][2],0,0,0))
                  atom_count += 1
            else:
               mol_id = atom[i]['mol_id']
               f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
               #.format(i+1,mol_id,3,0,atom[i]['coord'][0],atom[i]['coord'][1],atom[i]['coord'][2],0,0,0))
               .format(atom_count,mol_id,3,0,atom[i]['coord'][0],atom[i]['coord'][1],atom[i]['coord'][2],0,0,0))
               atom_count += 1
               # dummy atom for hydrogen, this is for rigid body fix
               f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
               .format(atom_count,mol_id,6,0,atom[i]['coord'][0]+0.001,atom[i]['coord'][1],atom[i]['coord'][2],0,0,0))
               atom_count += 1
        
        
        f.write("\nBonds\n\n")
        bond_count = 1
        atom_count = 1
        for i in atom:
            if atom[i]['resi'] == 'NTG':
               if i%3 == 2:
                  f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(bond_count,1,i-1,i+1))
                  bond_count += 1
                  f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(bond_count,1,i,i+1))
                  bond_count += 1
            else: 
                  f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(bond_count,3,atom_count,atom_count+1))
                  atom_count += 1
                  bond_count += 1
            atom_count += 1
                  

    return


# Wrapper function for writing the lammps data file which is read by the .in.init file during run initialization
def Write_data(Filename,Atom_types,Sim_Box,Elements,Geometry,Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,Dihedrals,Dihedral_types,Dihedral_params,\
               Impropers,Improper_types,Improper_params,Charges,VDW_params,Masses,Molecule,Improper_flag=False):

    # Write an xyz for easy viewing
    #with open(Filename+'/'+Filename.split('/')[-1]+'.xyz','w') as f:
    with open(Filename+'/in.xyz','w') as f:
        f.write('{}\n\n'.format(len(Geometry)))
        for count_i,i in enumerate(Geometry):
            f.write('{:20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(Elements[count_i],i[0],i[1],i[2]))

    # Create type dictionaries (needed to convert each atom,bond,angle, and dihedral type to consecutive numbers as per LAMMPS convention)
    # Note: LAMMPS orders atomtypes, bonds, angles, dihedrals, etc as integer types. Each of the following dictionaries holds the mapping between
    #       the true type (held in the various types lists) and the lammps type_id, obtained by enumerated iteration over the respective set(types).
    Atom_type_dict = {}
    for count_i,i in enumerate(sorted(set(Atom_types))):
        for j in Atom_types:
            if i == j:
                Atom_type_dict[i]=count_i+1
            if i in list(Atom_type_dict.keys()):
                break
    Bond_type_dict = {}
    for count_i,i in enumerate(sorted(set(Bond_types))):
        for j in Bond_types:
            if i == j:
                Bond_type_dict[i]=count_i+1
            if i in list(Bond_type_dict.keys()):
                break
    Angle_type_dict = {}
    for count_i,i in enumerate(sorted(set(Angle_types))):
        for j in Angle_types:
            if i == j:
                Angle_type_dict[i]=count_i+1
            if i in list(Angle_type_dict.keys()):
                break
    Dihedral_type_dict = {}
    for count_i,i in enumerate(sorted(set(Dihedral_types))):
        for j in Dihedral_types:
            if i == j:
                Dihedral_type_dict[i]=count_i+1
            if i in list(Dihedral_type_dict.keys()):
                break
    Improper_type_dict = {}
    for count_i,i in enumerate(sorted(set(Improper_types))):
        for j in Improper_types:
            if i == j:
                Improper_type_dict[i]=count_i+1
            if i in list(Improper_type_dict.keys()):
                break

    # Write the data file
    #with open(Filename+'/'+Filename.split('/')[-1]+'.data','w') as f:
    with open(Filename+'/in.data','w') as f:
        
        
        # Write system properties
        f.write("LAMMPS data file via vdw_self_gen.py, on {}\n\n".format(datetime.datetime.now()))

        f.write("{} atoms\n".format(len(Elements)))
        f.write("{} atom types\n".format(len(set(Atom_types))))
        if len(Bonds) > 0:
            f.write("{} bonds\n".format(len(Bonds)))
            f.write("{} bond types\n".format(len(set(Bond_types))))
        if len(Angles) > 0:
            f.write("{} angles\n".format(len(Angles)))
            f.write("{} angle types\n".format(len(set(Angle_types))))
        if len(Dihedrals) > 0:
            f.write("{} dihedrals\n".format(len(Dihedrals)))
            f.write("{} dihedral types\n".format(len(set(Dihedral_types))))
        if Improper_flag and len(Impropers) > 0:
            f.write("{} impropers\n".format(len(Impropers)))
            f.write("{} improper types\n".format(len(set(Improper_types))))
        f.write("\n")

        # Write box dimensions
        f.write("{:< 20.16f} {:< 20.16f} xlo xhi\n".format(Sim_Box[0],Sim_Box[1]))
        f.write("{:< 20.16f} {:< 20.16f} ylo yhi\n".format(Sim_Box[2],Sim_Box[3]))
        f.write("{:< 20.16f} {:< 20.16f} zlo zhi\n\n".format(Sim_Box[4],Sim_Box[5]))

        # Write Masses
        f.write("Masses\n\n")
        for count_i,i in enumerate(sorted(set(Atom_types))):
            for j in set(Atom_types):
                if Atom_type_dict[j] == count_i+1:
                    f.write("{} {:< 8.6f}\n".format(count_i+1,Masses[str(j)])) # count_i+1 bc of LAMMPS 1-indexing
        f.write("\n")

        # Write Atoms
        f.write("Atoms\n\n")
        for count_i,i in enumerate(Atom_types):
            f.write("{:<8d} {:< 4d} {:< 4d} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:< 20.16f} {:d} {:d} {:d}\n"\
            .format(count_i+1,Molecule[count_i],Atom_type_dict[i],Charges[count_i],Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],0,0,0))

        # Write Bonds
        if len(Bonds) > 0:
            f.write("\nBonds\n\n")
            for count_i,i in enumerate(Bonds):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Bond_type_dict[Bond_types[count_i]],i[0]+1,i[1]+1))

        # Write Angles
        if len(Angles) > 0:
            f.write("\nAngles\n\n")
            for count_i,i in enumerate(Angles):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Angle_type_dict[Angle_types[count_i]],i[0]+1,i[1]+1,i[2]+1))

        # Write Dihedrals
        if len(Dihedrals) > 0: 
            f.write("\nDihedrals\n\n")
            for count_i,i in enumerate(Dihedrals):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Dihedral_type_dict[Dihedral_types[count_i]],i[0]+1,i[1]+1,i[2]+1,i[3]+1))

        # Write Impropers
        if Improper_flag and len(Impropers) > 0: 
            f.write("\nImpropers\n\n")
            for count_i,i in enumerate(Impropers):
                f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Improper_type_dict[Improper_types[count_i]],i[0]+1,i[1]+1,i[2]+1,i[3]+1))

    # Write the settings file
    fixed_modes = {'bonds':[], 'angles':[]}
    #with open(Filename+'/'+Filename.split('/')[-1]+'.in.settings','w') as f:
    with open(Filename+'/in.settings','w') as f:

        # Write non-bonded interactions (the complicated form of this loop is owed
        # to desire to form a nicely sorted list in terms of the lammps atom types
        # Note: Atom_type_dict was initialize according to sorted(set(Atom_types) 
        #       so iterating over this list (twice) orders the pairs, too.
        f.write("     {}\n".format("# Non-bonded interactions (pair-wise)"))
        for count_i,i in enumerate(sorted(set(Atom_types))):     
            for count_j,j in enumerate(sorted(set(Atom_types))): 

                # Skip duplicates
                if count_j < count_i:
                    continue

                # Conform to LAMMPS i <= j formatting
                if Atom_type_dict[i] <= Atom_type_dict[j]:
                    f.write("     {:20s} {:<10d} {:<10d} ".format("pair_coeff",Atom_type_dict[i],Atom_type_dict[j]))
                else:
                    f.write("     {:20s} {:<10d} {:<10d} ".format("pair_coeff",Atom_type_dict[j],Atom_type_dict[i]))

                # Determine key (ordered by initialize_VDW function such that i > j)
                if i > j:
                    key = (i,j)
                else:
                    key = (j,i)

                # Write the parameters
                for k in VDW_params[key]:
                    if type(k) is str:
                        f.write("{:20s} ".format(k))
                    if type(k) is float:
                        f.write("{:< 20.6f} ".format(k))
                f.write("\n")                        

        # Write stretching interactions
        # Note: Bond_type_dict was initialized by looping over sorted(set(Bond_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Stretching interactions"))
        for i in sorted(set(Bond_types)):
            f.write("     {:20s} {:<10d} ".format("bond_coeff",Bond_type_dict[i]))
            for j in Bond_params[i]:
                if j == "fixed":
                    continue
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
            f.write("\n")

            # populate fixed_modes
            if Bond_params[i][0] == "fixed":
                fixed_modes["bonds"] += [Bond_type_dict[i]]

        # Write bending interactions
        # Note: Angle_type_dict was initialized by looping over sorted(set(Angle_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Bending interactions"))
        for i in sorted(set(Angle_types)):
            f.write("     {:20s} {:<10d} ".format("angle_coeff",Angle_type_dict[i]))
            for j in Angle_params[i]:
                if j == "fixed":
                    continue
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
            f.write("\n")

            # populate fixed_modes
            if Angle_params[i][0] == "fixed":
                fixed_modes["angles"] += [Angle_type_dict[i]]

        # Write dihedral interactions
        # Note: Dihedral_type_dict was initialized by looping over sorted(set(Dihedral_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        f.write("\n     {}\n".format("# Dihedral interactions"))
        for i in sorted(set(Dihedral_types)):
            f.write("     {:20s} {:<10d} ".format("dihedral_coeff",Dihedral_type_dict[i]))
            for j in Dihedral_params[i]:
                if type(j) is str:
                    f.write("{:20s} ".format(j))
                if type(j) is float:
                    f.write("{:< 20.6f} ".format(j))
                if type(j) is int:
                    f.write("{:< 20d} ".format(j))
            f.write("\n")

        # Write improper interactions
        # Note: Improper_type_dict was initialized by looping over sorted(set(Improper_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        if Improper_flag:
            f.write("\n     {}\n".format("# Improper interactions"))
            for i in sorted(set(Improper_types)):
                f.write("     {:20s} {:<10d} ".format("improper_coeff",Improper_type_dict[i]))
                for j in Improper_params[i][1:]:
                    if type(j) is str:
                        f.write("{:20s} ".format(j))
                    if type(j) is float:
                        f.write("{:< 20.6f} ".format(j))
                    if type(j) is int:
                        f.write("{:< 20d} ".format(j))
                f.write("\n")

    return Atom_type_dict,fixed_modes

# Writing lammps molecule file, for inserting molecules
#def Write_mol(Filename,Atom_types,Sim_Box,Elements,Geometry,Bonds,Bond_types,Angles,Angle_types,Dihedrals,Dihedral_types,\
#               Impropers,Improper_types,Charges,molname,Improper_flag=False):
def Write_mol(Filename,Data,Improper_flag=False):

    
    Atom_types     = []
    Bonds          = []
    Bond_types     = []
    Angles         = []
    Angle_types    = []
    Dihedrals      = []
    Dihedral_types = []
    Impropers      = []
    Improper_types = []
    Charges        = []
    Masses         = []
    atom_index = 0

    for mol in Data:
        
       Atom_types     = Atom_types + Data[mol]["Atom_types"]
       Bonds          = Bonds + [ (j[0]+atom_index,j[1]+atom_index) for j in Data[mol]["Bonds"] ]
       Bond_types     = Bond_types + Data[mol]["Bond_types"]
       Angles         = Angles + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index) for j in Data[mol]["Angles"] ]
       Angle_types    = Angle_types + Data[mol]["Angle_types"]
       Dihedrals      = Dihedrals + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index,j[3]+atom_index) for j in Data[mol]["Dihedrals"] ]
       Dihedral_types = Dihedral_types + Data[mol]["Dihedral_types"]
       Impropers      = Impropers + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index,j[3]+atom_index) for j in Data[mol]["Impropers"] ]
       Improper_types = Improper_types + Data[mol]["Improper_types"]
       Charges        = Charges + Data[mol]["Charges"]
       Masses         = Masses + [ Data[mol]["Masses"][j] for j in Data[mol]["Atom_types"] ] 
       atom_index += len(Data[mol]['Geometry'])

    # Create type dictionaries (needed to convert each atom,bond,angle, and dihedral type to consecutive numbers as per LAMMPS convention)
    # Note: LAMMPS orders atomtypes, bonds, angles, dihedrals, etc as integer types. Each of the following dictionaries holds the mapping between
    #       the true type (held in the various types lists) and the lammps type_id, obtained by enumerated iteration over the respective set(types).
    Atom_type_dict = {}
    for count_i,i in enumerate(sorted(set(Atom_types))):
        for j in Atom_types:
            if i == j:
                Atom_type_dict[i]=count_i+1
            if i in list(Atom_type_dict.keys()):
                break
    Bond_type_dict = {}
    for count_i,i in enumerate(sorted(set(Bond_types))):
        for j in Bond_types:
            if i == j:
                Bond_type_dict[i]=count_i+1
            if i in list(Bond_type_dict.keys()):
                break
    Angle_type_dict = {}
    for count_i,i in enumerate(sorted(set(Angle_types))):
        for j in Angle_types:
            if i == j:
                Angle_type_dict[i]=count_i+1
            if i in list(Angle_type_dict.keys()):
                break
    Dihedral_type_dict = {}
    for count_i,i in enumerate(sorted(set(Dihedral_types))):
        for j in Dihedral_types:
            if i == j:
                Dihedral_type_dict[i]=count_i+1
            if i in list(Dihedral_type_dict.keys()):
                break
    Improper_type_dict = {}
    for count_i,i in enumerate(sorted(set(Improper_types))):
        for j in Improper_types:
            if i == j:
                Improper_type_dict[i]=count_i+1
            if i in list(Improper_type_dict.keys()):
                break
    for mol in Data:
       # Write an xyz for easy viewing
       with open(Filename+'/'+mol+'.xyz','w') as f:
           f.write('{}\n\n'.format(len(Data[mol]['Geometry'])))
           for count_i,i in enumerate(Data[mol]['Geometry']):
               f.write('{:20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(Data[mol]['Elements'][count_i],i[0],i[1],i[2]))

       # Write the txt file
       with open(Filename+'/'+mol+'.txt','w') as f:
           
           # Write system properties
           f.write("# molecule file for {}\n\n".format(mol))

           f.write("{} atoms\n".format(len(Data[mol]['Elements'])))
           if len(Data[mol]['Bonds']) > 0:
               f.write("{} bonds\n".format(len(Data[mol]['Bonds'])))
           if len(Data[mol]['Angles']) > 0:
               f.write("{} angles\n".format(len(Data[mol]['Angles'])))
           if len(Data[mol]['Dihedrals']) > 0:
               f.write("{} dihedrals\n".format(len(Data[mol]['Dihedrals'])))
           if Improper_flag and len(Impropers) > 0:
               f.write("{} impropers\n".format(len(Impropers)))
           f.write("\n")

           # Write Coordinates
           f.write("Coords\n\n")
           for count_i in range(len(Data[mol]['Elements'])):
               f.write("{:<8d} {:< 20.16f} {:< 20.16f} {:< 20.16f} \n"\
               .format(count_i+1,Data[mol]['Geometry'][count_i,0],Data[mol]['Geometry'][count_i,1],Data[mol]['Geometry'][count_i,2]))
            
           # Write types
           f.write("\nTypes\n\n")
           for count_i in range(len(Data[mol]['Elements'])):
               f.write("{:<8d} {:< 4d} \n".format(count_i+1,Atom_type_dict[Data[mol]['Atom_types'][count_i]]))

           # Write charges
           f.write("\nCharges\n\n")
           for count_i in range(len(Data[mol]['Elements'])):
               f.write("{:<8d} {:< 20.16f}\n".format(count_i+1,Data[mol]['Charges'][count_i]))


           # Write Bonds
           if len(Data[mol]['Bonds']) > 0:
               f.write("\nBonds\n\n")
               for count_i,i in enumerate(Data[mol]['Bonds']):
                   f.write("{:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Bond_type_dict[Data[mol]['Bond_types'][count_i]],i[0]+1,i[1]+1))

           # Write Angles
           if len(Data[mol]['Angles']) > 0:
               f.write("\nAngles\n\n")
               for count_i,i in enumerate(Data[mol]['Angles']):
                   f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Angle_type_dict[Data[mol]['Angle_types'][count_i]],i[0]+1,i[1]+1,i[2]+1))

           # Write Dihedrals
           if len(Data[mol]['Dihedrals']) > 0: 
               f.write("\nDihedrals\n\n")
               for count_i,i in enumerate(Data[mol]['Dihedrals']):
                   f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Dihedral_type_dict[Data[mol]['Dihedral_types'][count_i]],i[0]+1,i[1]+1,i[2]+1,i[3]+1))

           # Write Impropers
           if Improper_flag and len(Data[mol]['Impropers']) > 0: 
               f.write("\nImpropers\n\n")
               for count_i,i in enumerate(Data[mol]['Impropers']):
                   f.write("{:<8d} {:<8d} {:<8d} {:<8d} {:<8d} {:<8d}\n".format(count_i+1,Improper_type_dict[Data[mol]['Improper_types'][count_i]],i[0]+1,i[1]+1,i[2]+1,i[3]+1))
      
    return 


def Num2Element(num):
    # Initialize periodic table

    periodic = {1: 'H', 2: 'HE', 3: 'LI', 4: 'BE', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'NE', 11: 'NA', 12: 'MG', 13: 'AL', 14: 'SI', 15: 'P', 16: 'S', 17: 'CL', 18: 'AR', 19: 'K', 20: 'CA', 21: 'SC', 22: 'TI', 23: 'V', 24: 'CR', 25: 'MN', 26: 'FE', 27: 'CO', 28: 'NI', 29: 'CU', 30: 'ZN', 31: 'GA', 32: 'GE', 33: 'AS', 34: 'SE', 35: 'BR', 36: 'KR', 37: 'RB', 38: 'SR', 39: 'Y', 40: 'ZR', 41: 'NB', 42: 'MO', 43: 'TC', 44: 'RU', 45: 'RH', 46: 'PD', 47: 'AG', 48: 'CD', 49: 'IN', 50: 'SN', 51: 'SB', 52: 'TE', 53: 'I', 54: 'XE', 55: 'CS', 56: 'BA', 72: 'HF', 73: 'TA', 74: 'W', 75: 'RE', 76: 'OS', 77: 'IR', 78: 'PT', 79: 'AU', 80: 'HG', 81: 'TL', 82: 'PB', 83: 'BI', 84: 'PO', 85: 'AT', 86: 'RN'}

    element = periodic[num]
    return element



# A wrapper for the commands to parse the bonds, angles, and dihedrals from the adjacency matrix.
# Returns:   list of atomtypes, bond_types, bond instances, angle_types, angle instances, dihedral_types,
#            diehdral instances, charges, and VDW parameters.
def Find_parameters(Adj_mat,Bond_mats,Geometry,Atom_types,FF_db="FF_file",Improper_flag = False, force_read=False,remove_multi=False):

    # List comprehension to determine bonds from a loop over the adjacency matrix. Iterates over rows (i) and individual elements
    # ( elements A[count_i,count_j] = j ) and stores the bond if the element is "1". The count_i < count_j condition avoids
    # redudant bonds (e.g., (i,j) vs (j,i) ). By convention only the i < j definition is stored.
    print("Parsing bonds...")
    Bonds          = [ (count_i,count_j) for count_i,i in enumerate(Adj_mat) for count_j,j in enumerate(i) if j == 1 ]
    Bond_types     = [ (Atom_types[i[0]],Atom_types[i[1]]) for i in Bonds ]

    # List comprehension to determine angles from a loop over the bonds. Note, since there are two bonds in every angle, there will be
    # redundant angles stored (e.g., (i,j,k) vs (k,j,i) ). By convention only the i < k definition is stored.
    print("Parsing angles...")
    Angles          = [ (count_j,i[0],i[1]) for i in Bonds for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in i ]
    Angle_types     = [ (Atom_types[i[0]],Atom_types[i[1]],Atom_types[i[2]]) for i in Angles ]

    # List comprehension to determine dihedrals from a loop over the angles. Note, since there are two angles in every dihedral, there will be
    # redundant dihedrals stored (e.g., (i,j,k,m) vs (m,k,j,i) ). By convention only the i < m definition is stored.
    print("Parsing dihedrals...")
    Dihedrals      = [ (count_j,i[0],i[1],i[2]) for i in Angles for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in i ]
    Dihedral_types = [ (Atom_types[i[0]],Atom_types[i[1]],Atom_types[i[2]],Atom_types[i[3]]) for i in Dihedrals ]

    # List comprehension to determine dihedrals from a loop over the angles. Note, since there are two angles in every dihedral, there will be
    # redundant dihedrals stored (e.g., (i,j,k,m) vs (m,k,j,i) ). By convention only the i < m definition is stored.
    if Improper_flag:    print("Parsing impropers...")
    Impropers      = [ (i[1],i[0],i[2],count_j) for i in Angles for count_j,j in enumerate(Adj_mat[i[1]]) if j == 1 and count_j not in i ]
    Improper_types = [ (Atom_types[i[0]],Atom_types[i[1]],Atom_types[i[2]],Atom_types[i[3]]) for i in Impropers ]

    # Canonicalize the modes
    for i in range(len(Bonds)):
        Bond_types[i],Bonds[i] = canon_bond(Bond_types[i],ind=Bonds[i])
    for i in range(len(Angles)):
        Angle_types[i],Angles[i] = canon_angle(Angle_types[i],ind=Angles[i])
    for i in range(len(Dihedrals)):
        Dihedral_types[i],Dihedrals[i] = canon_dihedral(Dihedral_types[i],ind=Dihedrals[i])
    for i in range(len(Impropers)):
        Improper_types[i],Impropers[i] = canon_improper(Improper_types[i],ind=Impropers[i])        

    # Remove redundancies    
    if len(Bonds) > 0: Bonds,Bond_types = list(map(list, list(zip(*[ (i,Bond_types[count_i]) for count_i,i in enumerate(Bonds) if count_i == [ count_j for count_j,j in enumerate(Bonds) if j == i or j[::-1] == i ][0]  ]))))
    if len(Angles) > 0: Angles,Angle_types = list(map(list, list(zip(*[ (i,Angle_types[count_i]) for count_i,i in enumerate(Angles) if count_i == [ count_j for count_j,j in enumerate(Angles) if j == i or j[::-1] == i ][0]  ]))))
    if len(Dihedrals) > 0: Dihedrals,Dihedral_types = list(map(list, list(zip(*[ (i,Dihedral_types[count_i]) for count_i,i in enumerate(Dihedrals) if count_i == [ count_j for count_j,j in enumerate(Dihedrals) if j == i or j[::-1] == i ][0]  ]))))
    if len(Impropers) > 0: Impropers,Improper_types = list(map(list, list(zip(*[ (i,Improper_types[count_i]) for count_i,i in enumerate(Impropers) if count_i == [ count_j for count_j,j in enumerate(Impropers) if j[0] == i[0] and len(set(i[1:]).intersection(set(j[1:]))) ][0] ]))))

    # Remove bimodal angles (e.g., if the same angle is present in
    # more than one configuration, then only the most abundant is kept)
    if remove_multi is True:
        for i in set(Angle_types):

            # Calculate most populated bin
            bins = [ 20.0*_ for _ in arange(9) ]
            hist = zeros(len(bins))
            angles = []
            for count_j,j in enumerate(Angle_types):
                if j == i:
                    v1 = Geometry[Angles[count_j][0]] - Geometry[Angles[count_j][1]]
                    v2 = Geometry[Angles[count_j][2]] - Geometry[Angles[count_j][1]]
                    v1 = v1/norm(v1)
                    v2 = v2/norm(v2)
                    angles += [ acos(dot(v1,v2)) * 180/pi ]
                    hist[int(angles[-1]/20.0)] += 1            
                else:
                    angles += [0.0]
            angle = bins[where(hist == max(hist))[0][0]]

            # Determine the angles to be removed
            del_inds = []
            for count_j,j in enumerate(Angle_types):                
                if j == i and abs(angles[count_j] - angle) > 20.0:
                    del_inds += [count_j]
            Angle_types = [ j for count_j,j in enumerate(Angle_types) if count_j not in del_inds ]
            Angles = [ j for count_j,j in enumerate(Angles) if count_j not in del_inds ]

    # Add the opls/harmonic type as a fifth element in the Dihedral_types tuples
    for count_i,i in enumerate(Dihedrals):
        if 2 in [ j[i[1],i[2]] for j in Bond_mats ]:
            Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["harmonic"])        
        # elif ( "R" in Dihedral_types[count_i][1] or "E" in Dihedral_types[count_i][1] or "Z" in Dihedral_types[count_i][1] ) and \
        #      ( "R" in Dihedral_types[count_i][2] or "E" in Dihedral_types[count_i][2] or "Z" in Dihedral_types[count_i][2] ):
        #     Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["harmonic"])        
        else:
            Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["opls"])                        

    ##############################################################
    # Read in parameters: Here the stretching, bending, dihedral #
    # and non-bonding interaction parameters are read in from    #
    # parameters file. Mass and charge data is also included.    #
    # The program looks for a simple match for the first entry   #
    # in each line with one of the bond or angle types.          #
    # INPUT: param_file, BOND_TYPES_LIST, ANGLE_TYPES_LIST       #
    #        DIHEDRAL_TYPES_LIST, ELEMENTS                       #
    # OUTPUT: CHARGES, MASSES, BOND_PARAMS, ANGLE_PARAMS,        #
    #         DIHERAL_PARAMS, PW_PARAMS                          #
    ##############################################################

    # Initialize dictionaries

    # Read in masses and charges
    Masses = {}
    with open(FF_db,'r') as f:
        content=f.readlines()
        
    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() == 'atom':
            Masses[fields[1]] = float(fields[3]) 

    # Read in bond parameters
    Bond_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0] == 'bond':
            if fields[3] == "harmonic":
                Bond_params[(fields[1],fields[2])] = [float(fields[4]),float(fields[5])]
            elif fields[3] == "fixed":
                Bond_params[(fields[1],fields[2])] = ["fixed",0.0,float(fields[4])]
            else:
                print("ERROR: only harmonic bond definitions are currently supported by gen_md_for_sampling.py. Exiting...")
                quit()

    # Read in angle parameters
    Angle_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() == 'angle':
            if fields[4] == "harmonic":
                Angle_params[(fields[1],fields[2],fields[3])] = [float(fields[5]),float(fields[6])]
            elif fields[4] == "fixed":
                Angle_params[(fields[1],fields[2],fields[3])] = ["fixed",0.0,float(fields[5])]
            else:
                print("ERROR: only harmonic angle definitions are currently supported by gen_md_for_sampling.py. Exiting...")
                quit()

    # Read in dihedral parameters
    Dihedral_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() in ['dihedral','torsion']:            
            if fields[5] == "opls":
                Dihedral_params[(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
            elif fields[5] == "harmonic":
                Dihedral_params[(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ]
            else:
                print("ERROR: Only opls and harmonic dihedral types are currently supported by gen_md_for_sampling.py. Exiting...")
                quit()
        
    # Read in improper parameters
    Improper_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() in ['improper']:
            if fields[5] == "harmonic":
                Improper_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),float(fields[7])]
            else:
                print("ERROR: Only opls type dihedral definitions are currently supported by gen_md_for_vdw.py. Exiting...")
                quit()
                
    # Search for charges based on atom type
    with open(FF_db,'r') as f:
        content=f.readlines()

    Charges = zeros(len(Atom_types))
    for i in range(len(Atom_types)):
        for lines in content:
            fields=lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['charge'] and Atom_types[i] == fields[1]:
                Charges[i] = float(fields[2])

    # Search for polarizability based on atom type
    with open(FF_db,'r') as f:
        content=f.readlines()

    alphas = zeros(len(Atom_types))
    for count_i,i in enumerate(alphas):
        alphas[count_i] = -1
    for i in range(len(Atom_types)):
        for lines in content:
            fields=lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['polar'] and Atom_types[i] == fields[1]:
                alphas[i] = float(fields[2])

    # Search for LP parameters
    LP_params = {}
    with open(FF_db,'r') as f:
        for lines in f:
            fields = lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['charge']:

                tag = split_string(fields[1])
                if tag[0] == 'L' and tag[1] == 'P':
                     LP_params[fields[1]] = float(fields[2])

    # Search for VDW parameters
    VDW_params = {}
    with open(FF_db,'r') as f:
        for lines in f:
            fields = lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['vdw']:

                # Only two parameters are required for lj types
                if fields[3] == "lj":
                    if fields[1] > fields[2]:
                        VDW_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                    else:
                        VDW_params[(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]
                elif fields[3] == "buck":
                    if fields[1] > fields[2]:
                        VDW_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5]),float(fields[6])]
                    else:
                        VDW_params[(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5]),float(fields[6])]

    # Check for missing parameters
    Missing_masses = [ i for i in Atom_types if str(i) not in list(Masses.keys()) ] 
    Missing_charges = [ count_i for count_i,i in enumerate(Charges) if i == -100.0 ]; Missing_charges = [ Atom_types[i] for i in Missing_charges ]
    Missing_bonds = [ i for i in Bond_types if (i[0],i[1]) not in list(Bond_params.keys()) ]
    Missing_angles = [ i for i in Angle_types if (i[0],i[1],i[2]) not in list(Angle_params.keys()) ]
    Missing_dihedrals = [ i for i in Dihedral_types if (i[0],i[1],i[2],i[3],i[4]) not in list(Dihedral_params.keys()) ]
    Missing_impropers = []
    if Improper_flag is True: Missing_impropers = [ i for i in Improper_types if (i[0],i[1],i[2],i[3]) not in list(Improper_params.keys()) ]

    # When the force_read option is set to True, the script will attempt to replace the missing dihedral types with whatever (matching)
    # parameters are in the force-field file. For example, if an opls type is expected but a harmonic type is supplied, the harmonic type
    # will be used IF THE FLAG IS TRUE, otherwise the script will print an error and exit.
    if force_read:

        # Assemble lists of types and keys present in the force field dictionary
        bt,bk = ( [ tuple(i[:2]) for i in list(Bond_params.keys()) ],     [ i for i in list(Bond_params.keys()) ] )
        at,ak = ( [ tuple(i[:3]) for i in list(Angle_params.keys()) ],    [ i for i in list(Angle_params.keys()) ] )
        dt,dk = ( [ tuple(i[:4]) for i in list(Dihedral_params.keys()) ], [ i for i in list(Dihedral_params.keys()) ] )

        # Replace missing dihedral types with whatever matching type is present in the supplied force-field
        for i in set(Missing_dihedrals):
            if i[:4] in dt:
                new_type = dk[dt.index(i[:4])]
                print("using dihedral type {} in place of {}...".format(new_type,i))
                r_inds = [ count_j for count_j,j in enumerate(Dihedral_types) if j == i ]
                for j in r_inds:
                    Dihedral_types[j] = new_type
                Missing_dihedrals = [ j for j in Missing_dihedrals if j != i ]


    # Print diagnostics on missing parameters and quit if the prerequisites are missing.
    if ( len(Missing_masses) + len(Missing_charges) + len(Missing_bonds) + len(Missing_angles) + len(Missing_dihedrals) + len(Missing_impropers) ) > 0:
        print("\nUh Oh! There are missing FF parameters...\n")

        if Missing_masses:
            print("Missing masses for the following atom types: {}".format([ i for i in set(Missing_masses) ]))
        if Missing_charges:
            print("Missing charges for the following atom types: {}".format([ i for i in set(Missing_charges) ]))
        if Missing_bonds:
            print("Missing bond parameters for the following bond types: {}".format([ i for i in set(Missing_bonds) ]))
        if Missing_angles:
            print("Missing angle parameters for the following angle types: {}".format([ i for i in set(Missing_angles) ]))
        if Missing_dihedrals:
            print("Missing dihedral parameters for the following dihedral types: {}".format([ i for i in set(Missing_dihedrals) ]))
        if Improper_flag and Missing_impropers:
            print("Missing improper parameters for the following improper types: {}".format([ i for i in set(Missing_impropers) ]))
        
        print("\nEnsure the specification of the missing parameters. Exiting...")
        quit()

    return list(Bonds),list(Bond_types),Bond_params,list(Angles),list(Angle_types),Angle_params,list(Dihedrals),list(Dihedral_types),Dihedral_params,list(Impropers),list(Improper_types),Improper_params,Charges.tolist(),Masses,VDW_params,alphas.tolist(),LP_params

def split_string(word): 
    return [char for char in word]  

# Description: finds the number of disconnected subnetworks in the 
#              adjacency matrix, which corresponds to the number of 
#              separate molecules.
#
# Inputs:      adj_mat: numpy array holding a 1 in the indices of bonded
#                        atom types. 
#
# Returns:     mol_count: scalar, the number of molecules in the adj_mat
def mol_count(adj_mat):
    
    # Initialize list of atoms assigned to molecules and a counter for molecules
    placed_idx = []    
    mol_count = 0

    # Continue the search until all the atoms have been assigned to molecules
    while len(placed_idx)<len(adj_mat):

        # Use sequential elements of the adj_mat as seeds for the spanning network search
        for count_i,i in enumerate(adj_mat):

            # Only proceed with search if the current atom hasn't been placed in a molecule
            if count_i not in placed_idx:

                # Increment mol_count for every new seed and add the seed to the list of placed atoms
                mol_count += 1               
                placed_idx += [count_i]
                
                # Find connections
                idx = [ count_j for count_j,j in enumerate(i) if j==1 and count_j not in placed_idx ]
                
                # Continue until no new atoms are found
                while len(idx) > 0:
                    current = idx.pop(0)
                    if current not in placed_idx:
                        placed_idx += [current]
                        idx += [ count_k for count_k,k in enumerate(adj_mat[current]) if k == 1 and count_k not in placed_idx ]
    return mol_count




# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def init_VDW(Data,fit_pairs,read_pairs,Pair_min_dict,fit,FF_db,VDW_dict_prev,UA=0,avoid_read_flag=0,mixing_rule="lb",verbose=True):

    if verbose is True:
        if avoid_read_flag == 0 and mixing_rule in ["none","lb"]:
            print("\nReading VDW parameters or initializing based on UFF values and Lorentz-Berthelot mixing rules\n")
        elif avoid_read_flag == 0 and mixing_rule in ["wh"]:
            print("\nReading VDW parameters or initializing based on UFF values and Waldman-Hagler mixing rules\n")
        elif avoid_read_flag == 1 and mixing_rule in ["none","lb"]:
            print("\nInitializing VDW parameters based on UFF values and Lorentz-Berthelot mixing rules\n")
        elif avoid_read_flag == 1 and mixing_rule in ["wh"]:
            print("\nInitializing VDW parameters based on UFF values and Waldman-Hagler mixing rules\n")

        print("\t{:85s}  {:22s} {:20s}".format("pair_interaction","params","origin"))

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # Initialize initial_guess dictionary based on the initial_params.db file (if present). 
    # The default is to use these parameters in place of UFF if the file is present.     
    INIT_dict = {}
    if os.path.isfile('initial_params.db') and avoid_read_flag == 0:
        with open('initial_params.db','r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "vdw":
                    INIT_dict[(fields[1],fields[2])] = [ float(i) for i in fields[4:] ]
                    INIT_dict[(fields[2],fields[1])] = [ float(i) for i in fields[4:] ]

    # Collect types that will be converted into UA
    if UA==1:
        UA_types = {}
        UA_H_types = []
        for i in list(Data.keys()):

            # Iterate over the molecule a adj_mat and identify hydrogens attached to carbons. 
            for count_j,j in enumerate(Data[i]["adj_mat_a"]):
                if int(Data[i]["types_a"][count_j].split('[')[1].split(']')[0]) != 6: continue
                UA_H_types += [ Data[i]["types_a"][count_k] for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_a'][count_k].split('[')[1].split(']')[0]) == 1 ]

            # Iterate over the molecule b adj_mat and identify hydrogens attached to carbons.
            for count_j,j in enumerate(Data[i]["adj_mat_b"]):
                if int(Data[i]["types_b"][count_j].split('[')[1].split(']')[0]) != 6: continue
                UA_H_types += [ Data[i]["types_b"][count_k] for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_b'][count_k].split('[')[1].split(']')[0]) == 1 ]

    # Initialize VDW_dict with fit types, first guess based on element types and Lorentz-Berthelot mixing rules
    VDW_dict = {}
    for i in fit_pairs:

        # Grab the base atomic number from the atomtype
        type_1 = int(i[0].split('[')[1].split(']')[0])
        type_2 = int(i[1].split('[')[1].split(']')[0])
        
        # for UA fits, the eps of Hydrogen containing pairs are set to zero
        # so that they do not participate in the fitting. 
        if UA==1 and ( i[0] in UA_H_types or i[1] in UA_H_types ):
            eps    = 0.0
            sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0
            VDW_dict[i] = (eps,sigma)
            VDW_dict[(i[1],i[0])] = (eps,sigma)

            # Print diagnostic
            if verbose is True:
                print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"Zeroed"))

        # Non-hydrogen UA types read their values from the AA-fit
        elif UA==1:
            eps    = VDW_dict_prev[i][0]
            sigma  = VDW_dict_prev[i][1]
            VDW_dict[i] = (eps,sigma)
            VDW_dict[(i[1],i[0])] = (eps,sigma)

            # Print diagnostic
            if verbose is True:
                print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"AA-Fit"))

        # For all atom fits, initial guesses for eps and sigma are based on scaled UFF parameters with a 
        # check to make sure that the radii aren't too small for the fit set. 
        else:
            if i in INIT_dict:
                eps    = INIT_dict[i][0]
                sigma  = INIT_dict[i][1]
                if sigma < Pair_min_dict[i]*2.**(-1./6.): sigma = Pair_min_dict[i]*2.**(-1./6.)
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)
                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"initial_params.db"))
            else:
                if mixing_rule in ["lb",'none']:
                    eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) 
                    sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0
                elif mixing_rule == "wh":
                    sigma  = ((UFF_dict[type_1][1]**(6.0)+UFF_dict[type_2][1]**(6.0))/2.0)**(1.0/6.0)
                    eps    = (UFF_dict[type_1][0]*UFF_dict[type_1][1]**(6.0) * UFF_dict[type_2][0]*UFF_dict[type_2][1]**(6.0) )**(0.5) / sigma**(6.0)
                else: 
                    print("ERROR in init_VDW: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
                    quit()
                if sigma < Pair_min_dict[i]*2.**(-1./6.): sigma = Pair_min_dict[i]*2.**(-1./6.)
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)
        
                # Print diagnostic
                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"UFF"))

    # For diagnostic purposes, the fixed_params.db and user_supplied FF files are reparsed here in order to 
    # print the proper origin information for the user.
    # Read in fixed parameters from the fixed_params.db file
    FIXED_dict = {"vdws":{}}
    if os.path.isfile('fixed_params.db'):         
        FIXED_dict = get_FF_data(['fixed_params.db'])

    # If a FF_db file is supplied then the program avoids fitting any of the
    # pairs that are pairs are in the database
    READ_dict = {"vdws":{}}
    if FF_db != []:        
        READ_dict = get_FF_data(FF_db,mixing_rule=mixing_rule)

    # Initialize read_pairs with entries from the FF_db.
    for i in read_pairs:

        # Grab the base atomic number from the atomtype
        type_1 = int(i[0].split('[')[1].split(']')[0])
        type_2 = int(i[1].split('[')[1].split(']')[0])
        
        # for UA fits, the eps of Hydrogen containing pairs are set to zero
        # so that they do not participate in the fitting. 
        if UA==1 and ( i[0] in UA_H_types or i[1] in UA_H_types ):
            VDW_dict[i] = (0.0,1.0)
            VDW_dict[(i[1],i[0])] = (0.0,1.0)
            continue
        else:        

            # Mixing rule protocols
            if mixing_rule != "none": 

                # Use mixing rule protocol for read parameters
                if (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):

                    eps_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for fixed parameters
                elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):

                    eps_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for fixed/read parameters
                elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):

                    eps_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for read/fixed parameters
                elif (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):

                    eps_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][4])

                # Apply mixing rule
                if mixing_rule in ["lb"]:
                    eps    = (eps_1*eps_2)**(0.5)
                    sigma  = (sigma_1+sigma_2)/2.0
                elif mixing_rule == "wh":
                    sigma  =  ( ( sigma_1**(6.0) + sigma_2**(6.0) ) / 2.0 )**(1.0/6.0)
                    eps    =  ( eps_1*sigma_1**(6.0) * eps_2*sigma_2**(6.0) )**(0.5) / sigma**(6.0)    
                else: 
                    print("ERROR in init_VDW: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
                    quit()
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)

                # Print commands
                if verbose is True:
                    if (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"supplied_file"))
                    elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed_params.db"))
                    elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed/supplied mixing"))
                    elif (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed/supplied mixing"))



            # Add the pair type (along with the reverse combination) using the parameters found in the READ_dict
            elif (i[0],i[1],'lj') in list(READ_dict["vdws"].keys()) or (i[1],i[0],'lj') in list(READ_dict["vdws"].keys()):
                if (i[0],i[1],'lj') in list(READ_dict["vdws"].keys()):
                    VDW_dict[i] = (float(READ_dict["vdws"][(i[0],i[1],'lj')][3]),float(READ_dict["vdws"][(i[0],i[1],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(READ_dict["vdws"][(i[0],i[1],'lj')][3]),float(READ_dict["vdws"][(i[0],i[1],'lj')][4]))            

                else:
                    VDW_dict[i] = (float(READ_dict["vdws"][(i[1],i[0],'lj')][3]),float(READ_dict["vdws"][(i[1],i[0],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(READ_dict["vdws"][(i[1],i[0],'lj')][3]),float(READ_dict["vdws"][(i[1],i[0],'lj')][4]))

                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"supplied_file"))

            # Add the pair type (along with the reverse combination) using the parameters found in the FIXED_dict
            elif (i[0],i[1],'lj') in list(FIXED_dict["vdws"].keys()) or (i[1],i[0],'lj') in list(FIXED_dict["vdws"].keys()):

                if (i[0],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                    VDW_dict[i] = (float(FIXED_dict["vdws"][(i[0],i[1],'lj')][3]),float(FIXED_dict["vdws"][(i[0],i[1],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(FIXED_dict["vdws"][(i[0],i[1],'lj')][3]),float(FIXED_dict["vdws"][(i[0],i[1],'lj')][4]))            

                else:
                    VDW_dict[i] = (float(FIXED_dict["vdws"][(i[1],i[0],'lj')][3]),float(FIXED_dict["vdws"][(i[1],i[0],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(FIXED_dict["vdws"][(i[1],i[0],'lj')][3]),float(FIXED_dict["vdws"][(i[1],i[0],'lj')][4]))

                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed_params.db"))

            else:
                print("ERROR in init_VDW: Information for pair type {} was expected in the the fixed_params.db".format((i[0],i[1],'lj')))
                print("                   and/or the supplied FF files, but it could not be found. Please check that the parameters exist. Exiting...")
                quit()


    # Perform conversions to buckingham style. Conversions follow the approach of Rappe et al. 
    # with a "form factor" of 12.0 which reproduces the long-range LJ behavior
    if fit == 'buck':
        for i in list(VDW_dict.keys()):
            A = VDW_dict[i][0]*0.5*exp(12.0)
            B = VDW_dict[i][1]/12.0
            C = VDW_dict[i][0] * 2.0 * VDW_dict[i][1]**(6.0)
            VDW_dict[i] = (A,B,C)

    return VDW_dict

# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def initialize_VDW(atomtypes,sigma_scale=1.0,eps_scale=1.0,VDW_FF={},Force_UFF=0,mixing_rule='lb'):

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # NEW: Initialize VDW_dict first guess based on element types and Lorentz-Berthelot mixing rules
    # Order of operations: (1) If the parameters are in the supplied FF database then they are used as is
    # (2) If the self-terms are in the supplied FF datebase then they are used to generate the mixed
    # interactions (3) UFF parameters are used. 
    VDW_dict = {}
    origin = {}
    VDW_styles = []
    for count_i,i in enumerate(atomtypes):
        for count_j,j in enumerate(atomtypes):
            if count_i < count_j:
                continue

            # Check for parameters in the database
            if (i,j) in VDW_FF and Force_UFF != 1:

                # Determine appropriate lammps style
                if VDW_FF[(i,j)][0] == "lj":
                    VDW_type = "lj/cut/coul/long"
                elif VDW_FF[(i,j)][0] == "buck":
                    VDW_type = "buck/coul/long"
                else:
                    print("ERROR in initialize_VDW: only lj and buck pair types are supported. Exiting...")
                    quit()

                # Assign style
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type] + VDW_FF[(i,j)][1:]
                else:
                    VDW_dict[(j,i)] = [VDW_type] + VDW_FF[(i,j)][1:]
                origin[(i,j)] = origin[(j,i)] = "read"

            # Check for reverse combination
            elif (j,i) in VDW_FF and Force_UFF != 1:

                # Determine appropriate lammps style
                if VDW_FF[(j,i)][0] == "lj":
                    VDW_type = "lj/cut/coul/long"
                elif VDW_FF[(j,i)][0] == "buck":
                    VDW_type = "buck/coul/long"
                else:
                    print("ERROR in initialize_VDW: only lj and buck pair types are supported. Exiting...")
                    quit()

                # Assign style
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type] + VDW_FF[(j,i)][1:]
                else:
                    VDW_dict[(j,i)] = [VDW_type] + VDW_FF[(j,i)][1:]
                origin[(i,j)] = origin[(j,i)] = "read"

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1 and mixing_rule == 'lb':

                # Check compatibility with mixing rules
                if VDW_FF[(i,i)][0] != "lj" or VDW_FF[(j,j)][0] != "lj":
                    print("ERROR in initialize_VDW: only lj styles support mixing rules. Exiting...")
                    quit()

                # Apply mixing rules and assign
                VDW_type = "lj/cut/coul/long"
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(j,j)][1])**(0.5) * eps_scale
                sigma  = (VDW_FF[(i,i)][2]+VDW_FF[(j,j)][2])/2.0 * sigma_scale
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]
                origin[(i,j)] = origin[(j,i)] = "lb"

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1 and mixing_rule == 'wh':

                # Check compatibility with mixing rules
                if VDW_FF[(i,i)][0] != "lj" or VDW_FF[(j,j)][0] != "lj":
                    print("ERROR in initialize_VDW: only lj styles support mixing rules. Exiting...")
                    quit()

                # Apply mixing rules and assign
                VDW_type = "lj/cut/coul/long"
                sigma  = ((VDW_FF[(i,i)][2]**(6.0)+VDW_FF[(j,j)][2]**(6.0))/2.0)**(1.0/6.0)
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(i,i)][2]**(6.0) * VDW_FF[(j,j)][1]*VDW_FF[(j,j)][2]**(6.0) )**(0.5) / sigma**(6.0)
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]
                origin[(i,j)] = origin[(j,i)] = "wh"

            # Last resort: Use UFF parameters.
            else:
                VDW_type = "lj/cut/coul/long"
                type_1 = int(i.split('[')[1].split(']')[0])
                type_2 = int(j.split('[')[1].split(']')[0])
                eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) * eps_scale
                sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0 * sigma_scale
                if i > j:
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]
                origin[(i,j)] = origin[(j,i)] = "UFF"
                
            # Collect a list of the LAMMPS styles used in the simulation
            VDW_styles += [VDW_type]

    # Print summary
    print("\n{}".format("*"*177))
    print("* {:^173s} *".format("Initializing VDW parameters for the simulation (those with * were read from the FF file(s))"))
    print("*{}*".format("-"*175))
    print("* {:<50s} {:<50s} {:<20s}  {:<18s} {:<18s} {:<8s}    *".format("Type","Type","VDW_type","eps (kcal/mol)","sigma (angstroms)","origin"))
    print("{}".format("*"*177))
    for j in list(VDW_dict.keys()):
        print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f}  {:<18s}".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2],origin[j]))
    print("")


#    return VDW_dict,sorted(set(VDW_styles))
    return VDW_dict
# Description: A wrapper for the commands to generate a cubic box and array of molecules for the lammps run
def Pack_Box(Data,Box,atom,Box_offset=0,Density=0.0,N_Density=0.0):
    
    print("\nPacking simulation cell...")

    low = -Box/2
    high = Box/2
    # Define sim box
    Sim_Box = array([low,high,low,high,low,high])

    # Intialize lists for iterating over molecules and keeping track of how many have been placed
    keys = list(Data.keys())             # list of unique molecule keys
    placed_num = [0]*len(keys)     # list for keeping track of how many of each molecule have been placed
    atom_index = 0                 # an index for keeping track of how many atoms have been placed
    mol_index = 0                  # an index for keeping track of how many molecules have been placed

    # Initialize various lists to hold the elements, atomtype, molid labels etc.
    # Create a list of molecule ids for each atom        
    Geometry_sim       = zeros([sum([ Data[i]["N_mol"]*len(Data[i]["Geometry"]) for i in list(Data.keys()) ]),3])
    Adj_mat_sim        = zeros([sum([ Data[i]["N_mol"]*len(Data[i]["Geometry"]) for i in list(Data.keys()) ]),sum([ Data[i]["N_mol"]*len(Data[i]["Geometry"]) for i in list(Data.keys()) ])])
    Molecule_sim       = []
    Molecule_files     = []
    Elements_sim       = [] 
    Atom_types_sim     = []
    Bonds_sim          = []
    Bond_types_sim     = []
    Angles_sim         = []
    Angle_types_sim    = []
    Dihedrals_sim      = []
    Dihedral_types_sim = []
    Impropers_sim      = []
    Improper_types_sim = []
    Charges_sim        = []
    Masses_sim         = []

    for count_k,k in enumerate(keys):
        
         for i in range(0,Data[k]["N_mol"]):
             Adj_mat_sim[atom_index:(atom_index+len(Data[k]["Geometry"])),atom_index:(atom_index+len(Data[k]["Geometry"]))] = Data[k]["Adj_mat"]
             Molecule_sim       = Molecule_sim + [mol_index]*len(Data[k]["Elements"])
             Molecule_files     = Molecule_files + [k]
             Elements_sim       = Elements_sim + Data[k]["Elements"]
             Atom_types_sim     = Atom_types_sim + Data[k]["Atom_types"]
             Bonds_sim          = Bonds_sim + [ (j[0]+atom_index,j[1]+atom_index) for j in Data[k]["Bonds"] ]
             Bond_types_sim     = Bond_types_sim + Data[k]["Bond_types"]
             Angles_sim         = Angles_sim + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index) for j in Data[k]["Angles"] ]
             Angle_types_sim    = Angle_types_sim + Data[k]["Angle_types"]
             Dihedrals_sim      = Dihedrals_sim + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index,j[3]+atom_index) for j in Data[k]["Dihedrals"] ]
             Dihedral_types_sim = Dihedral_types_sim + Data[k]["Dihedral_types"]
             Impropers_sim      = Impropers_sim + [ (j[0]+atom_index,j[1]+atom_index,j[2]+atom_index,j[3]+atom_index) for j in Data[k]["Impropers"] ]
             Improper_types_sim = Improper_types_sim + Data[k]["Improper_types"]
             Charges_sim        = Charges_sim + Data[k]["Charges"]
             Masses_sim         = Masses_sim + [ Data[k]["Masses"][j] for j in Data[k]["Atom_types"] ] 
             for j in range(atom_index,atom_index+len(Data[k]["Geometry"])): 
               Geometry_sim[j] = atom[j]['coord']
             atom_index += len(Data[k]["Geometry"])
             # Increment atom_index based on the number of atoms in the current geometry
             mol_index += 1



    return Elements_sim,Atom_types_sim,Geometry_sim,Bonds_sim,Bond_types_sim,Angles_sim,Angle_types_sim,Dihedrals_sim,Dihedral_types_sim,Impropers_sim,Improper_types_sim,Charges_sim,Molecule_sim,Molecule_files,Adj_mat_sim,Sim_Box

# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def initialize_VDW_old(atomtypes,sigma_scale=1.0,eps_scale=1.0,VDW_type="lj/cut/coul/long",VDW_FF={},Force_UFF=0):

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # NEW: Initialize VDW_dict first guess based on element types and Lorentz-Berthelot mixing rules
    # Order of operations: (1) If the parameters are in the supplied FF database then they are used as is
    # (2) If the self-terms are in the supplied FF datebase then they are used to generate the mixed
    # interactions (3) UFF parameters are used. 
    VDW_dict = {}
    for count_i,i in enumerate(atomtypes):
        for count_j,j in enumerate(atomtypes):
            if count_i < count_j:
                continue

            # Check for parameters in the database
            if (i,j) in VDW_FF and Force_UFF != 1:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
            elif (j,i) in VDW_FF and Force_UFF != 1:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1:
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(j,j)][1])**(0.5) * eps_scale
                sigma  = (VDW_FF[(i,i)][2]+VDW_FF[(j,j)][2])/2.0 * sigma_scale
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]

            # Last resort: Use UFF parameters.
            else:
                type_1 = int(i.split('[')[1].split(']')[0])
                type_2 = int(j.split('[')[1].split(']')[0])
                eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) * eps_scale
                sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0 * sigma_scale
                if i > j:
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]

    # Print summary
    print("\n{}".format("*"*167))
    print("* {:^163s} *".format("Initializing VDW parameters for the simulation (those with * were read from the FF file(s))"))
    print("*{}*".format("-"*165))
    print("* {:<50s} {:<50s} {:<20s}  {:<18s} {:<18s}   *".format("Type","Type","VDW_type","eps (kcal/mol)","sigma (angstroms)"))
    print("{}".format("*"*167))
    for j in list(VDW_dict.keys()):
        if j in VDW_FF and Force_UFF != 1:
            print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f} *".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2]))
        else:
            print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f}".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2]))
    print("")


    return VDW_dict

# Description:   A wrapper for the commands to create a data dictionary holding the FF, mode, and geometry information for each molecule
# Inputs:  Filename:
# Returns: Data: a dictionary keyed to the list of coord_files
#                a keyed entry is itself a dictionary holding the geometry, 
#                elements, charges, bonds, bond_types, angles, angle_types,
#                dihedrals, dihedral_types and masses.
#                subdictionary keys (e.g., Data[*]["Geometry"])
#                Geometry: array of the molecule's geometry
#                Elements: list of the molecule's elements (indexed to the rows of "Geometry")
#                Atom_types: list of the molecule's atom_types (indexed to the rows of "Geometry")
#                Charges: list of the molecule's charges (indexed to the rows of "Geometry")
#                Bonds: a list of the molecule's bonds, each element being a tuple of the atom indices that are bonded
#                Angles: a list of the molecule's angles, each element being a tuple of the atom indices that form an angle mode
#                Dihedrals: a list of the molecule's dihedrals, each element being a tuple of the atom indices that form a dihedral
#                Bond_types: a list of the molecule's bond_types, indexed to Bonds
#                Angle_types: a list of the molecule's angle_types, indexed to Angles
#                Dihedral_types: a list of the molecule's dihedral_types, indexed to Dihedrals
#                Impropers: a list of the molecule's improper dihedrals by atom tindex
#                Improper_types: a list of the molecule's improper_types, indexed to Impropers
#
def get_data(FF_all,N_mol,atom,q_list,gens,Improper_flag=False,force_read_opt=False,remove_multi=False):

    # Initialize dictionary to hold the FF, mode, and geometry information about each molecule
    Data = {}

    # Iterate over all molecules being simulated and collect their FF, mode, and geometry information.
    for count_i,i in enumerate(list(N_mol.keys())):

        # Initialize dictionary for this geometry and set number of molecules to place in the simulated system
        Data[i] = {}      
        Data[i]["N_mol"] = N_mol[i]

        # Extract Element list and Coord list from the file
        Data[i]["Elements"] = []
        start = False
        for a in atom.keys():
            if atom[a]['resi'] == i:
               if start is False:
                  mol_id_0 = atom[a]['mol_id']
                  Data[i]["Elements"].append(atom[a]['elem'])
                  Data[i]["Geometry"] = atom[a]['coord']
                  start = True
                  continue
               else:
                  if atom[a]['mol_id'] != mol_id_0: break
                  Data[i]["Elements"].append(atom[a]['elem'])
                  Data[i]["Geometry"] = np.vstack((Data[i]["Geometry"],atom[a]['coord']))

        # Generate adjacency table
        Data[i]["Adj_mat"] = Table_generator(Data[i]["Elements"],Data[i]["Geometry"])

        # Automatically parse the atom types based on the adjacency matrix and gens argument. 
        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        Hybridizations = Hybridization_finder(Data[i]["Elements"],Data[i]["Adj_mat"],force_opt=True)
    
        # Find atom types
        print("Determining atom types based on a {}-bond deep search...".format(i,gens))
        Data[i]["Atom_types"] = id_types(Data[i]["Elements"],Data[i]["Adj_mat"],gens)

        # Determine bonding matrix for the compound
        #if q_list[count_i] == "none":
        #    q_tmp = 0.0
        #else:
        #    q_tmp = q_list[count_i]
        q_tmp = 0.0
        Data[i]["Bond_mats"] = find_lewis(Data[i]["Elements"], Data[i]["Adj_mat"], q_tot=int(q_tmp),b_mat_only=True,verbose=False)        

        # Check the number of molecules
        mol_in_in = mol_count(Data[i]["Adj_mat"])
        if mol_in_in > 1:
            print("ERROR: {} molecules were discovered in geometry {}. Check the geometry of the input file. Exiting...".format(mol_in_in,i))
            quit()

        # Generate list of bonds angles and dihedrals    
        print("\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Parsing Modes and FF Information for Molecule {}".format(i),"*"*167))
        Data[i]["Bonds"],Data[i]["Bond_types"],Data[i]["Bond_params"],Data[i]["Angles"],Data[i]["Angle_types"],Data[i]["Angle_params"],\
        Data[i]["Dihedrals"],Data[i]["Dihedral_types"],Data[i]["Dihedral_params"],Data[i]["Impropers"],Data[i]["Improper_types"],Data[i]["Improper_params"],Data[i]["Charges"],Data[i]["Masses"],Data[i]["VDW_params"],Data[i]['alphas'],Data[i]['LP_params'] =\
            Find_parameters(Data[i]["Adj_mat"],Data[i]["Bond_mats"],Data[i]["Geometry"],Data[i]["Atom_types"],FF_db=FF_all,Improper_flag = Improper_flag, force_read=force_read_opt,remove_multi=remove_multi)

        # Print System characteristics
        print("\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Mode Summary for Molecule {}".format(i),"*"*167))
        print("\nAtom_types ({}):\n".format(len(set(Data[i]["Atom_types"]))))
        for j in sorted(set(Data[i]["Atom_types"])):
            print("\t{}".format(j))
        print("\nBond types ({}):\n".format(len(set(Data[i]["Bond_types"]))))
        for j in sorted(set(Data[i]["Bond_types"])):
            print("\t{}".format(j))
        print("\nAngle types ({}):\n".format(len(set(Data[i]["Angle_types"]))))
        for j in sorted(set(Data[i]["Angle_types"])):
            print("\t{}".format(j))
        print("\nDihedral types ({}):\n".format(len(set(Data[i]["Dihedral_types"]))))
        for j in sorted(set(Data[i]["Dihedral_types"])):
            print("\t{}".format(j))
        if Improper_flag:
            print("\nImproper types ({}):\n".format(len(set(Data[i]["Improper_types"]))))
            for j in sorted(set(Data[i]["Improper_types"])):
                print("\t{}".format(j))
        
        # Subtract off residual
        print("\n{:40s} {}".format("Residual Charge:",sum(Data[i]["Charges"])))
        #if q_list[count_i] == "round":
        #    q_tot = int(round(sum(Data[i]["Charges"])))
        #elif q_list[count_i] != "none":
        #    q_tot = int(q_list[count_i])

        # Avoid subtracting off the residual if q_list is none
        #if q_list[count_i] != "none":
        #    correction = (float(q_tot)-sum(Data[i]["Charges"]))/float(len(Data[i]["Atom_types"]))
        #    print("{:40s} {}".format("Charge added to each atom:",correction))
        #    for j in range(len(Data[i]["Atom_types"])):
        #        Data[i]["Charges"][j] += correction

        print("{:40s} {}".format("Final Total Charge:",sum(Data[i]["Charges"])))
        print("\n{}".format("*"*192))
        print("* {:^188s} *".format("System Characteristics"))
        print("*{}*".format("-"*190))
        print("* {:<87s} {:<24s} {:<24s} {:<25s} {:<25s}*".format("Type","Element","Mass","Charge","Polarizability"))
        print("{}".format("*"*192))
        for j in range(len(Data[i]["Atom_types"])):
            print(" {:<88s} {:<23s} {:< 24.6f} {:< 24.6f} {:< 24.6f}".format(Data[i]["Atom_types"][j],Data[i]["Elements"][j],Data[i]["Masses"][str(Data[i]["Atom_types"][j])],Data[i]["Charges"][j],Data[i]["alphas"][j]))

    return Data

# A simple wrapper function for generating a fix statement for scaled charges
def scale_charges(charge_scale,Atom_type_dict,Atom_types,Charges):

    fixes = []
    if charge_scale != 1.0:
        fixes = '# Setting up fixes to scale the charges by {}\n'.format(charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):
            fixes += 'group charge_{} type {}\n'.format(i,i)        
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):
            index = Atom_types.index([ k for k in list(Atom_type_dict.keys()) if Atom_type_dict[k] == i ][0])
            fixes += 'variable charge_{} equal {}\n'.format(i,Charges[index]*charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):        
            index = Atom_types.index([ k for k in list(Atom_type_dict.keys()) if Atom_type_dict[k] == i ][0])
            fixes += 'fix {} charge_{} adapt 0 atom charge v_charge_{} reset yes\n'.format(i,i,i)
        fixes += '\n'
    return fixes

                    

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

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
