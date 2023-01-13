#!/bin/env python

import sys,os,argparse
import numpy as np

def main(argv):

    # Define argument parser
    parser = argparse.ArgumentParser(description='This program rescales the coordinates and simulation box of a lammps data file.')

    parser.add_argument('datafile', help = 'Name of the datafile to being rescaled.')

    parser.add_argument('-o', dest='output_name', default=None,
                        help = 'Optional output name for the rescaled data file. (default: overwrite original)')

    parser.add_argument('-d_N', dest='d_N', default=0.05,
                        help = 'Number density value in units of atoms/A^3 for the rescaled simulation. (default: 0.05, a low condensed phase density)')

    # Parse arguments
    args=parser.parse_args(argv)    
    args.d_N = float(args.d_N)

    # Check that the file exists
    if os.path.isfile(args.datafile) is False:
        print("ERROR: {} does not exist. Exiting...".format(args.datafile))
        quit()

    # Parse the save path
    if len(args.datafile.split("/")) == 1:
        save_path = ""
    else:
        save_path = "/".join(args.datafile.split("/")[:-1])

    # Parse and rescale the coords with the output directed to .tmp.rescale 
    with open(args.datafile,'r') as f:
        box=np.zeros([3,2])
        atom_flag = 0
        with open('.tmp.rescale','w') as o:            
            for lines in f:
                fields = lines.split()

                # Parse the number of atoms in the simulation
                if len(fields) == 2 and fields[1] == "atoms":
                    N_atoms  = int(fields[0])

                # Parse the box information
                if len(fields) == 4 and fields[2] == 'xlo':
                    box[0,:] = np.array([float(fields[0]),float(fields[1])])
                    continue

                if len(fields) == 4 and fields[2] == 'ylo':
                    box[1,:] = np.array([float(fields[0]),float(fields[1])])
                    continue

                if len(fields) == 4 and fields[2] == 'zlo':
                    box[2,:] = np.array([float(fields[0]),float(fields[1])])
                    
                    # Print the new box
                    scale_factor = ( ( float(N_atoms)/( (box[0,1]-box[0,0]) * (box[1,1]-box[1,0]) * (box[2,1]-box[2,0]) ) ) / args.d_N )**(1.0/3.0)
                    box = box * scale_factor
                    o.write("{:< 16.8f} {:< 16.8f} xlo xhi\n".format(box[0,0],box[0,1]))
                    o.write("{:< 16.8f} {:< 16.8f} ylo yhi\n".format(box[1,0],box[1,1]))
                    o.write("{:< 16.8f} {:< 16.8f} zlo zhi\n".format(box[2,0],box[2,1]))
                    continue

                # Check for the start of the coords section
                if len(fields) > 0 and fields[0] == "Atoms":
                    atom_flag = 1

                # Check for the end of the coords section (non-essential but just in case
                if len(fields) > 0 and fields[0] == "Bonds":
                    atom_flag = 0

                # Rescale the coords
                if atom_flag == 1 and len(fields) >= 7:
                    o.write("{:8} {:8} {:8} {:< 16.8f} {:< 16.8f} {:< 16.8f} {:< 16.8f} {}\n".format(fields[0],fields[1],fields[2],float(fields[3]),\
                                                                                             float(fields[4])*scale_factor,float(fields[5])*scale_factor,float(fields[6])*scale_factor,\
                                                                                             " ".join([ "{:>4}".format(j) for j in fields[7:] ])))
                # Write lines as they are found
                else:
                    o.write(lines)

    # Move the temporary datafile to its final location
    if args.output_name is None:
        os.rename(".tmp.rescale",args.datafile)
    else:
        os.rename(".tmp.rescale",save_path+'/'+args.output_name)

if __name__ == "__main__":
   main(sys.argv[1:])
