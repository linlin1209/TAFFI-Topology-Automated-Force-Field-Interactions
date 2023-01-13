#!/bin/env python                                                                                                                                                             

import sys,argparse,subprocess,os,time,math
import jenerate_resp_points as GRC # Generate Resp Charge 
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from file_parsers import xyz_parse
from subprocess import PIPE
import numpy as np


#### output map.xyz : in Bohr (for ORCA)
#### output connoly_vis.xyz (in Anstrom)


def main(argv):

    parser = argparse.ArgumentParser(description='Reads in xyz files and generate connolly grid for ORCA to reevaluate and a xyz file for visualization')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    parser.add_argument('-d', dest='density',type=float, default=6,
                        help = 'Sets density for generating the grid points (default: 6)')

    parser.add_argument('-radii_factor', dest='size_factor',type=float, default=10,
                        help = 'radii factor for VDW radii, this is only used for test default: 10 (10 is VDW_radii*1.0)')

    parser.add_argument('-gen', dest='gens',type=int, default=2,
                        help = 'taffi atom-type generation, default: 2')

    parser.add_argument('-qtot', dest='qtot',type=int, default=0,
                        help = 'total charge on molecule, default: 0 (neutral species)')

    parser.add_argument('-o', dest='file_name', default='charge',
                        help = 'output xyz name (default: charge.xyz)')

    parser.add_argument('-extra', dest='extra',type=float, default=0.1,
                        help = 'extra density factor for placing charges in gap this is only used for test(default: 0.1)')

    parser.add_argument('--test', dest='test', default=False, action='store_const', const=True,
                        help = 'when this falg is on, test mode is True, will report grid pts generated only, this is for getting extra_scale for placing gap charges')


    global args
    args=parser.parse_args(argv)    

    # Find files
    args.coord_files = args.coord_files.split()
    wild_card_files  = [ i for i in args.coord_files if "*" in i ]
    args.coord_files = [ i for i in args.coord_files if "*" not in i ]
    for i in wild_card_files:
        path = '/'.join(i.split('/')[:-1])
        if len(path) == 0:
          path = "."
        args.coord_files += [ path+"/"+files for files in os.listdir(path) if fnmatch(files,i) ] 

    # Handle "./" condition
    for count_i,i in enumerate(args.coord_files):
        if i[:2] == "./": args.coord_files[count_i] = args.coord_files[count_i][2:]    
    args.coord_files = list(set(args.coord_files))

    # Perform some consistency checks
    if False in [ i.split('.')[-1] == 'xyz' for i in args.coord_files ]:
        print( "ERROR: Check to ensure that the input file(s) are in .xyz format.")
        return
    if False in [ os.path.isfile(i) for i in args.coord_files ]:
        print( "ERROR: Could not find file {}. Exiting...".format(next([ i for i in args.coord_files if os.path.isfile(i) == False ])))


    # Create geo_dict for function input
    geo_dict_in = {}
    for f in args.coord_files:
    
      elements,geo = xyz_parse(f)
      inchikey = GetInchi(elements,geo)
      geo_dict_in[inchikey] = {}
      geo_dict_in[inchikey]["elements"] = elements
      geo_dict_in[inchikey]["geo"] = geo
    
    total_pts = fun(geo_dict_in)
    
    return total_pts 



def fun(geo_dict):
    
    for xyz_file in geo_dict:        
        inchikey_in = xyz_file

        # Print the name
        print( "  {}".format(xyz_file))


        if(args.test):
            test_pt = grid_gen(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["geo"],args.size_factor,args.extra)
            return

        scale_factors = [ 13, 22, 30, 50, 60 ]
        density_scales = [ 6.0, 1.1, 1.3, 0.6, 0.2]
        grid_points = {}
        total_points = 0
        for count_i,i in enumerate(scale_factors):
            grid_points[i] = grid_gen(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["geo"],scale_factors[count_i],density_scales[count_i])
            total_points += len(grid_points[i])

        # additional density scale factor, this is by trial and error to get an approximate grid pts that are similar to the final # of charges placed reported by Roux's group
        #convert to Bohr
        xyzscale = 1.889725989
        with open('map.xyz','w') as f:
            f.write("{}\n".format(total_points))
            for key in grid_points:
               for i in grid_points[key]:
                  f.write("{} {} {}\n".format(i[0]*xyzscale,i[1]*xyzscale,i[2]*xyzscale))

        with open('connolly_vis.xyz','w') as f:
            f.write("{} \n\n".format(total_points+len(geo_dict[xyz_file]["elements"])))
            for count_i,i in enumerate(geo_dict[xyz_file]["geo"]):
               f.write("{} {} {} {}\n".format(geo_dict[xyz_file]["elements"][count_i],i[0],i[1],i[2]))
            tmp = ['Li','Na','F', 'Cl', 'Br']
            for count_k,key in enumerate(grid_points):
               for i in grid_points[key]:
                  f.write("{} {} {} {}\n".format(tmp[count_k],i[0],i[1],i[2]))
            
    return total_points 

def grid_gen(elements,geo,size_factor,density_scale):
 
    conf = [] 
    for count_i,i in enumerate(elements):
        conf.append(dict(element=i,xyz=tuple(geo[count_i]))) 
    
    # Use Atom class in jenerate_resp.py 
    atoms = [GRC.Atom(**a) for a in conf]
    x, y, z = np.array([c['xyz'] for c in conf]).T
    #elements = [c['element'] for c in conf]
    # array of points for each scale of molecular surface
    #try different density, 6 is a good number
    points = []
    
    # Size factors (multiplication constants) multiplied by the VDW radii of corresponding atoms
    # sizefactor and density scale are set to be the same as 2005 Anisimov et al.
    # doi:10.1021/ct049930p
    for atom in atoms:
        points.extend(atom.radius_points(density=float(args.density)*density_scale, scale=size_factor/10, check_against=atoms))

    points = np.array(points) 
    rows = points.shape[0]
    print("Total of {} grid pts generated for Radii scale at {}".format(rows,float(size_factor/10)))

    
    return points

# XXX need to switch this to a more robust method that calls open babel's python binding 
def GetInchi(Elements,Geometry): #convert geo2inchi key

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('inchi.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    # create inchikey
    p = subprocess.Popen("obabel -ixyz inchi.xyz -oinchikey".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
    inchikey, err = p.communicate()
    inchikey = str(inchikey, 'utf-8').split('\n')[0] ## XXX might be a problem if there are characters utf-8 can't handle
    # remove xyz file
    os.remove("inchi.xyz")
   
    return inchikey


# Logger object redirects standard output to a file.
class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/place_charge.log", "a",buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass
   
            
if __name__ == "__main__":
   main(sys.argv[1:])



