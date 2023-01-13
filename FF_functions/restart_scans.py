#!/bin/env python

import sys,argparse,os,ast,re,fnmatch,matplotlib,subprocess,shutil,itertools,time

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *
from paramgen import *

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a .xyz file and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('folder', help = 'Folder for recursive search of REDO files in need of resubmission.')
   
    # optional arguments
    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'Controls the bond search depth for identifying unique atomtypes (default: 2)')

    parser.add_argument('-c', dest='constraints', default=" ", 
                        help = 'Designates the contraints to apply during the mode scans. "bond" "angle" and "dihedral" are valid options. The -c argument should be supplied as a '+\
                               'space delimited string. Based of the supplied argument(s), the corresponding modes are constrained during the mode scans.') 
    parser.add_argument('--local_relaxation', dest='local_relaxation', default=0, action='store_const',const=1,
                        help = 'This flag toggles the omission of constraints coincident with the dihedral modes being parametrized. For example, when performing the scan of a dihedral with -c '+\
                               '"bond angle dihedral" enabled, ordinarily all degrees of freedom are constrained. But with this flag enabled, bonds and angles involving any of the atoms in the '+\
                               'dihedral, and any coincident dihedrals, are allowed to relax. (Default: off)')

    parser.add_argument('-f', dest='functional', default='B3LYP',
                        help = 'Sets the functional for the calculation. (default: B3LYP; other typical options are WB97X-D3, and M062X)')

    parser.add_argument('-b', dest='basis', default='def2-TZVP',
                        help = 'Sets the basis set for the calculation. (default: def2-TZVP)')

    parser.add_argument('--no_D3', dest='D3_flag', default=1, action='store_const', const=0,
                        help = 'When this flag is present, the D3-dispersion correction is disabled. By default the dispersion correction is used.')

    parser.add_argument('--frag', dest='frag_opt', default=0, action='store_const',const=1,
                        help = 'This flag toggles the use of molecular fragments for the parametrizations. The smallest fragment consistent with each mode is retained for each parametrization. '+\
                               'For bonds and angles, the fragmentation is based upon the -gens flag. For dihedrals, one heavy atom beyond the 1 and 4 atoms is retained. Hydrogenation is used to '+\
                               'correct any dangling bonds. (Default: off)')

    # Make relevant inputs lowercase
    args=parser.parse_args(argv)    
    args.gens = int(args.gens)
    args.constraints = args.constraints.split()

    # Automatically turn bond, angle, and dihedral into bonds, angles, and dihedrals.
    for count_i,i in enumerate(args.constraints):
       if args.constraints[count_i][-1] != "s": args.constraints[count_i] += "s"

    # Consistency checks
    if False in [ i in ["bonds","angles","dihedrals"] for i in args.constraints ]: print("ERROR: only bonds, angles, and dihedrals are valid constraint options. Exiting..."); quit()

    # Find the REDO folders
    REDOS = [ dp+'/'+dirs for dp,dn,fn in os.walk(args.folder) for dirs in dn if dirs == "REDO" and os.path.isfile(dp+'/'+dirs+'/geo_opt/geo_opt.xyz') and os.path.isfile(dp+'/'+dirs+'/geo_opt/geo_opt.out') and os.path.isfile(dp+'/'+dirs+'/constraints.txt') ] 

    # Remove jobs that already have the resubmitted scan jobs present
    del_ind = set([ count_i for count_i,i in enumerate(REDOS) if "0" in [ j for j in os.listdir(i) if os.path.isdir(i+'/'+j) ] ])
    REDOS = [ i for count_i,i in enumerate(REDOS) if count_i not in del_ind ]

    # Generate resubmission jobs
    for i in REDOS:
       
        # Check that the geometry optimization ran to completion
        comp_flag = 0
        with open(i+'/geo_opt/geo_opt.out','r') as f:
            for lines in f:
                if "****ORCA TERMINATED NORMALLY****" in lines:
                    comp_flag = 1
        if comp_flag == 0:
            print("{} has an unoptimized geometry. Skipping...".format(i))
            continue

        # Parse the geometry, generate adj_mat, and id types
        elements,geo = xyz_parse(i+'/geo_opt/geo_opt.xyz')
        adj_mat = Table_generator(elements,geo)

        # Parse orca_input
        orca_dict = orca_in_parse(i+'/geo_opt/geo_opt.in')

        # Get constraint info
        cons = []
        with open(i+'/constraints.txt','r') as f:
            for lc,lines in enumerate(f):
                fields = lines.split()
                if lc == 0:               
                    N_steps = (int(fields.pop(-1))-1)/2   # NOTE: N_steps is defined in the gen_*_files programs as the number of steps in each direction not the total number of steps
                    step = float(fields.pop(-1))*float(N_steps)
                    cons += [[ int(j) for j in fields ]]
                else:
                    cons += [[ int(j) for j in fields ]]               

        # Grab atom types for this scan from the parent frag-*.modes file
        #mode_file = os.path.abspath(i)
        #frag = next( j for j in os.path.abspath(i).split('/')[::-1] if "frag" in j).split("_")[0]+".modes" # -> frag-0.modes
        #frag = "/".join(os.path.abspath(i).split('/')[:-1])+'/../../../'+frag # -> angles/[ ...]/../../../frag-0.modes
        
        frag = "out_"+ next( j for j in os.path.abspath(i).split('/')[::-1] if "-N" in j).split("_")[0]+".modes" # may not work if the last two string of the inchi key is not '-N' 
        mode_files = []
        inchi_keys = [ j for j in os.listdir('..') if os.path.isdir('../'+j) == True ]
        for inchi in inchi_keys:
            for filename in os.listdir('../{}'.format(inchi)):
                if fnmatch.fnmatch(filename, frag):
                    mode_files.append( os.path.join(os.path.abspath('../{}'.format(inchi)), frag))
        
        print("Looking for mode in files: {}".format(mode_files))
        
        sought = tuple(os.path.abspath(i).split('/')[-2].split('_'))
        modes= {}
        bond_flag = 0
        angle_flag = 0
        dihedral_flag = 0
        if len(cons[0]) == 2: bond_flag = 1
        if len(cons[0]) == 3: angle_flag = 1
        if len(cons[0]) == 4: dihedral_flag = 1        
        for mode_file in mode_files:
            #with open(frag,'r') as f:
            with open(mode_file,'r') as f:
                parse_mode = 0
                counter = 0
                for lines in f:
                    fields = lines.split()
                    if len(fields) == 0:
                        continue

                    # Look for the start of a relevant mode
                    if bond_flag == 1 and len(fields) == 2 and fields[0] == "bond" and fields[1] == "start":
                        parse_mode = 1
                        continue
                    if angle_flag == 1 and len(fields) == 2 and fields[0] == "angle" and fields[1] == "start":                
                        parse_mode = 1
                        continue
                    if dihedral_flag == 1 and len(fields) == 2 and fields[0] == "harmonic_dihedral" and fields[1] == "start":                                    
                        parse_mode = 1
                        continue

                    # Look for the end of a relevant mode
                    if bond_flag == 1 and len(fields) == 2 and fields[0] == "bond" and fields[1] == "end":
                        parse_mode = 0
                        continue
                    if angle_flag == 1 and len(fields) == 2 and fields[0] == "angle" and fields[1] == "end":                
                        parse_mode = 0
                        continue
                    if dihedral_flag == 1 and len(fields) == 2 and fields[0] == "harmonic_dihedral" and fields[1] == "end":                                    
                        parse_mode = 0
                        continue
                    
                    # Parse mode(s)
                    if parse_mode == 1:
                        ind = list(range(int(fields[0])))
                        parse_mode += 1
                        continue
                    if parse_mode == 2:
                        start_ind = len(list(modes.keys()))
                        for j in ind: modes[start_ind+j] = {"con":tuple([int(k) for k in fields[j].split("_")]),"atomtypes":[]}
                        parse_mode += 1
                        continue
                    if parse_mode == 3:
                        for j in ind:
                            modes[start_ind+j]["atomtypes"] += [fields[j]]
        
        # Grab the atomtypes that match this mode scan
        for j in list(modes.keys()):
            scanned = tuple([ modes[j]["atomtypes"][k] for k in modes[j]["con"] ])
            if scanned[::-1] == sought or scanned == sought:
                atom_types = modes[j]["atomtypes"]
                break
        
        # Generate scan (NOTE: by using step as the disp the restarted scans are narrower than the original submissions)
        if len(cons[0]) == 2:
            gen_bond_files(i,args.gens,adj_mat,cons[0],elements,geo,atom_types,bond_step=step/N_steps,bond_disp=step,\
                           charge=int(orca_dict["0"]["charge"]),multiplicity=int(orca_dict["0"]["multiplicity"]),procs=orca_dict["0"]["N_proc"],\
                           theory="DFT",constraints=args.constraints,local_relaxation=args.local_relaxation,frag_opt=args.frag_opt,functional=args.functional,basis=args.basis,D3_flag=args.D3_flag)

        if len(cons[0]) == 3:      
            gen_angle_files(i,args.gens,adj_mat,cons[0],elements,geo,atom_types,angle_step=step/N_steps,angle_disp=step,\
                            charge=int(orca_dict["0"]["charge"]),multiplicity=int(orca_dict["0"]["multiplicity"]),procs=orca_dict["0"]["N_proc"],\
                            theory="DFT",constraints=args.constraints,local_relaxation=args.local_relaxation,frag_opt=args.frag_opt,functional=args.functional,basis=args.basis,D3_flag=args.D3_flag)

        if len(cons[0]) == 4:      
            gen_harmonic_dihedral_files(i,args.gens,adj_mat,cons[0],elements,geo,atom_types,charge=int(orca_dict["0"]["charge"]),multiplicity=int(orca_dict["0"]["multiplicity"]),procs=orca_dict["0"]["N_proc"],\
                                        theory="DFT",constraints=args.constraints,local_relaxation=args.local_relaxation,frag_opt=args.frag_opt,\
                                        dihedral_disp=step,dihedral_step=step/N_steps,functional=args.functional,basis=args.basis,D3_flag=args.D3_flag)
        
        # Move the new scan
        #src  = "/".join(i.split('/')[:-1]+[i.split('/')[-2]])
        src  = "/".join(i.split('/')[:-1]+[i.split('/')[-2]])
        dst  = "/".join(i.split('/')[:-2]+[".tmp.tmp.tmp"])
        dst2 = "/".join(i.split('/')[:-1])
        
        shutil.move(src,dst)
        shutil.rmtree(dst2)
        shutil.move(dst,dst2)
        
        # Rename the files to match TAFFI paramgen.py formatting
        base = "-".join(dst2.split("/")[-1].split("_"))
        parent = dst2.split("/")[-3]
        for dp,dn,fn in os.walk(dst2):           
           for count_i,i in enumerate(fn):
               tail = i.split('.')[-1] 
               num=".".join(i.split(".")[:-1])
               num="_".join(num.split("_")[1:])              
               os.rename(dp+'/'+i,dp+'/{}_{}_{}.{}'.format(parent,base,num,tail))

# main sentinel
if __name__ == "__main__":
   main(sys.argv[1:])
