#!/bin/env python                                                                                                                                                              
import sys,os,argparse,shutil,getpass

def main(argv):
    

    user = getpass.getuser()

    parser = argparse.ArgumentParser(description='takes in the dimer xyz file in configs folder, do drude relaxation, this program assumes there is a FF folder that contains lin.prm and lin.rtf in vdw folder')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-xyz', dest='xyz', default='fitcharge.xyz',
                        help = 'xyz file name, default: fitcharge.xyz')
    parser.add_argument('--silent',dest='silent', default=False, const=True, action='store_const',
                        help = 'Shortcircuits all script related print statements') 
    parser.add_argument('-tag',dest='tag',default='',
                        help = 'tag for FF folder and relax drude, this assumes FF name is FF_tag and output drude: _drude_tag.xyz, default:None')

    args=parser.parse_args(argv)    

    title = args.xyz.split('.')[0]

    if not args.silent: 
      print("PROGRAM CALL: {} \n".format(' '.join([ i for i in sys.argv]))) 

    current_dir = os.getcwd()
    FFpath = '{}/FF{}'.format('/'.join(current_dir.split('/')[:-2]),args.tag)
    drude_dir = '{}/relax_drude{}'.format(current_dir,args.tag)
    xyzname = '{}/{}_drude{}.xyz'.format(current_dir,title,args.tag)

    if not os.path.isdir(FFpath):
         print("FF folder not found. Exiting....")
         quit()

    # check if drude has already been relaxed
    if os.path.isfile(xyzname):
         print("Drude for {} has been relaxed !".format(args.xyz))
         return
    if not os.path.isdir(drude_dir):  os.makedirs(drude_dir)
    shutil.copy2('{}/lin.prm'.format(FFpath),'{}/lin.prm'.format(drude_dir))
    shutil.copy2('{}/lin.rtf'.format(FFpath),'{}/lin.rtf'.format(drude_dir))


    os.chdir(drude_dir)

    # read in nonpolar config dimer xyz file
    prev_layer = '/'.join(os.getcwd().split('/')[:-1]) 
    elem_a,elem_b,geo_a,geo_b,types_a,types_b = scrape_xyz(prev_layer+'/'+args.xyz,parse_charges=0)
    
    # write CHARMM relaxation inp
    Write_inp(elem_a,elem_b,geo_a,geo_b)

    # Run charmm drude relaxation
    Write_exe()
    command = 'sh relax.sh'
    output = subprocess.Popen(command.split(),stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()

    # Read relax geo and charge and write it in xyz format
    charges = read_charge()
    prev_layer = '/'.join(os.getcwd().split('/')[:-1]) 
    Write_drude_xyz(charges,xyzname,elem_a,elem_b,types_a,types_b)
   
    os.chdir(prev_layer)


    if not os.path.isfile(xyzname):
         print("{} relaxation failed".format(args.xyz))
         quit()
   
    # if successful, rm relax folder
    else:
         shutil.rmtree(drude_dir)
    
    return

# Return the coordinates and atom types from an xyz file
# atom types are expected as a fifth column and molecule
# index is expected in the sixth column. These are atypical
# an xyz file, but automatically outputed by the vdw_gen.py 
# program.
def scrape_xyz(name,parse_charges=0):

    atom_count = 0
    with open(name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()

            # Grab the number of atoms in the geometry and initialize the parsing lists
            if lc == 0:
                elem = ["X"]*int(fields[0])
                geo = np.zeros([int(fields[0]),3]) 
                types = ["X"]*int(fields[0])
                a_ind = []
                b_ind = []
                
                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges =np.zeros([int(fields[0])])

            # Within the geometry block parse the geometry, charges, mol_ids, and atom types from the enhanced xyz file
            if lc > 1 and len(fields) >= 7:
                geo[atom_count,:] = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                if fields[4] == '0':
                    a_ind += [atom_count]
                else:
                    b_ind += [atom_count]
                elem[atom_count] = fields[0]

                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges[atom_count] = float(fields[5])
                types[atom_count] = fields[6]
                atom_count += 1    

    # If flag is enabled return the charges
    if parse_charges == 1:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               [ i for count_i,i in enumerate(types) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(types) if count_i in b_ind ],\
               charges[a_ind],charges[b_ind]

    # Don't return the charges list if the flag isn't enabled
    if parse_charges == 0:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               [ i for count_i,i in enumerate(types) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(types) if count_i in b_ind ]

def Write_inp(elem_a,elem_b,geo_a,geo_b):
      
    with open('relax.inp','w') as f:
         f.write('! this is necessart to work with non-interger total charge\n'+\
                 'bomlev -1\n\n'+\
                 'read rtf card name lin.rtf\n'+\
                 'read para card name lin.prm\n'+\
                 'read sequence lin 2\n\n'+
      
                 'generate lin first none last none setup warn drude dmass 0.4\n\n'+\

                 'read coor card free\n'+\
                 '* lin.qpos.1.crd\n'+\
                 '*\n')
         f.write('  {}\n'.format(len(elem_a)*2))
         for count_i,i in enumerate(geo_a):
           atomname = elem_a[count_i]+str(count_i+1)
           f.write("{:< 4d} {:< 4d} LIN {:} {:< 10.6f} {:< 10.6f} {:< 10.6f} LIN {:d} {:7.5f}\n"\
            .format(count_i+1,1,atomname,i[0],i[1],i[2],1,0.00000))
         for count_i,i in enumerate(geo_b):
           atomname = elem_a[count_i]+str(count_i+1)
           f.write("{:< 4d} {:< 4d} LIN {:} {:< 10.6f} {:< 10.6f} {:< 10.6f} LIN {:d} {:7.5f}\n"\
            .format(count_i+1+len(elem_a),2,atomname,i[0],i[1],i[2],1,0.00000))
         f.write('\n\ncoor sdrude\n'+\
                 'coor shake\n\n\n'+\

                 'update inbfrq -1 ihbfrq 0 -'+\
                 '  switch atom vatom vswitch cutnb 990. ctofnb 950. ctonnb 900.\n'+\
                 'cons fix sele .not. type D* end ! fix atoms other than Drude(relax Drude only)\n'+\
                 'mini ABNR nstep 1000 nprint 20\n'+\
                 'coor print\n\n'+\

                 'scalar charge show\n'+\
                 'open write unit 10 card name drude.crd\n'+\
                 'write coor card unit 10\n')
    return

def Write_exe():
    with open('relax.sh','w') as f:
         f.write('#!/bin/bash\n')
         f.write('charmm < relax.inp > relax.log\n')
         f.write('wait')
    return

def read_charge():
    with open('relax.log','r') as f:
        flag = 0
        charges = []
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 4: continue
            if fields[1] == 'scalar' and fields[2] == 'charge' and fields[3] == 'show':
               flag = 1
               continue
            if flag == 1 and fields[1] != 'LIN':
               flag = 0
               break 
            if flag == 1:
               charges.append(float(fields[6]))
               continue
    return charges

def Write_drude_xyz(charges,filename,elem_a,elem_b,types_a,types_b):
    with open('drude.crd','r') as f:
        count = 0
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if lc < 3: continue
            if lc == 3:
               natoms = int(fields[0])
               ele = []
               mol = []
               geos = np.zeros([natoms,3])
               flag = 1
               continue
            if flag == 1:
               geos[count] = np.array([ float(i) for i in fields[4:7] ])
               tmpele = []
               tmp = split_str(fields[3])
               for i in tmp:
                  if i.isnumeric(): break
                  else: tmpele += i
               ele.append("".join(tmpele))
               mol.append(int(fields[1])-1)
               count += 1
               if count == natoms: break
    with open(filename,'w') as f:
         atom_count = 0
         f.write('{} \n\n'.format(natoms))
         for count_i,i in enumerate(ele):
            f.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {:< 5d} {:< 12.6f} \n'.format(i,geos[count_i,0],geos[count_i,1],geos[count_i,2],mol[count_i],charges[count_i]))
    return

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder, "a",buffering = 1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
