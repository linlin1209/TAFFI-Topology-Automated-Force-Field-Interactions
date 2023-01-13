#!/bin/env python                                                                                                                                                             
import sys,os,argparse
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from file_parsers import xyz_parse
from adjacency import Table_generator
from id_types import Hybridization_finder




def main(argv):

    parser = argparse.ArgumentParser(description='Reads in CHARMM fitcharge output log file and transform it db file')

    #required (positional) arguments                                                                                                  

    parser.add_argument('-logname', dest='logfile', default='fitcharge.log',
                        help = 'input log filename, default: fitcharge.log')

    parser.add_argument('-xyz', dest='xyz', default='fitcharge.xyz',
                        help = 'xyz file used to determine atom type, default: fitcharge.xyz')


    parser.add_argument('-o', dest='outfile', default='fit_CHARMM.db',
                        help = 'output db filename, default: fit_CHARMM.db')

    parser.add_argument('-gens', dest='gens',type=int, default=2,
                        help = 'atomtypes generation')

    # Make parse inputs
    args=parser.parse_args(argv)

    if os.path.isfile(args.logfile) is False:
         print("ERROR: {} file not found".format(args.logfile))
         quit()
    if os.path.isfile(args.xyz) is False:
         print("ERROR: {} file not found".format(args.xyzfile))
         quit()

    element,geo = xyz_parse(args.xyz)

    with open(args.logfile,'r') as f:
        charge_polar_dict = {}
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 2 and flag == 0: continue
      
            if flag == 0  and fields[0] == 'chi2' and fields[1] == '(potential)':
               try:
                  chi2 = float(lines.split('=')[-1])  
               except:
                  chi2 = 0
          
            if flag == 0  and fields[0] == 'Final' and fields[1] == 'atomic' and fields[2] == 'charges':
               flag = 1
               continue

            if flag == 0  and fields[0] == 'TOTAL' and fields[1] == 'POLARIZABILITY' :
               flag = 2
               continue

            if flag == 1 :
               if fields[0] == 'Total':
                  flag = 0
                  continue 
               atom_name = fields[0]
               charge_polar_dict[atom_name] = {'charge':None,'polar':'NA'}
               charge_polar_dict[atom_name]['charge'] = float(fields[1])
               if len(fields) == 3:
                  charge_polar_dict[atom_name]['polar'] = float(fields[2])
               continue

            if flag == 2 and len(fields) == 5:
               polar = float(fields[4])
               flag = 0 
               break
    print("Molecular polarizability: {:10.6f} (Angstrom^3)".format(polar))
    print("chi2: {:10.6f} ".format(chi2))
               
    adj = Table_generator(element,geo)
    hybrid = Hybridization_finder(element,adj,force_opt=True)
    # one bond depth atom type (this makes making neighbor list very easy)
    atomT = id_types(element,adj,args.gens)
    with open(args.outfile,'w') as f:

        f.write("\n# Charge definitions\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","Charge"))
        for i in set(atomT):
            idx = [ count_k for count_k,k in enumerate(atomT) if k == i ]
            atomnames = list(charge_polar_dict.keys()) 
            f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(i),charge_polar_dict[atomnames[idx[0]]]['charge']))
        f.write("\n# Atomic Polarizability\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","polarizability(A^3)"))
        for i in set(atomT):
            idx = [ count_k for count_k,k in enumerate(atomT) if k == i ]
            atomnames = list(charge_polar_dict.keys()) 
            if charge_polar_dict[atomnames[idx[0]]]['polar'] != 'NA':
               f.write("{:10s} {:40s} {:< 20.6f}\n".format("polar",str(i),charge_polar_dict[atomnames[idx[0]]]['polar']))
            else:
               f.write("{:10s} {:40s} {:< 20.6f}\n".format("polar",str(i),0.0))


    return chi2,polar


if __name__ == "__main__":
   main(sys.argv[1:])
