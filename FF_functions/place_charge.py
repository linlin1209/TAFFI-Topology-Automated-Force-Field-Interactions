#!/bin/env python                                                                                                                                                             

import sys,argparse,subprocess,os,math
from subprocess import PIPE
from scipy.spatial.distance import cdist
import numpy as np
import jenerate_resp_points as GRC # Generate Resp Charge 
# Add TAFFI Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from file_parsers import xyz_parse
from adjacency import Table_generator,find_lewis,Find_modes_baonly
from id_types import Hybridization_finder


def main(argv):

    parser = argparse.ArgumentParser(description='place perturbation charges for input xyz file according to Roux\'s group method')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'name of xyz file')

    parser.add_argument('-d', dest='density',type=float, default=6,
                        help = 'Sets density for generating the grid points (default: 6)')

    parser.add_argument('-radii_factor', dest='size_factor',type=float, default='10',
                        help = 'radii factor for VDW radii default: 10 (10 is VDW_radii*1.0)')

    parser.add_argument('-d_factor', dest='density_scale_factor',type=float, default=1,
                        help = 'density scale factor (default: 1)')
    
    parser.add_argument('-o', dest='file_name', default='charge',
                        help = 'output xyz name (default: charge.xyz)')

    parser.add_argument('-ofolder', dest='folder_name', default='pointcharge',
                        help = 'output folder name (default: pointcharge)')

    parser.add_argument('-c', dest='criteria',type=float, default=5.0,
                        help = 'criteria for overlap intersection vector (default: 5 degree)')

    parser.add_argument('-extra', dest='extra',type=float, default=0.1,
                        help = 'extra density factor for placing charges in gap (default: 0.1)')

    parser.add_argument('-qtot', dest='qtot',type=int, default=0,
                        help = 'total charge on molecule, default: 0 (neutral species)')

    parser.add_argument('-gens', dest='gens',type=int, default=2,
                        help = 'generation for TAFFI atomtype, default: 2')

    parser.add_argument('--test', dest='test', default=False, action='store_const', const=True,
                        help = 'when this flag is on, test mode is True, will report grid pts generated only, this is for getting extra_scale for placing gap charges')


    # Make relevant inputs lowercase
    args=parser.parse_args(argv)    

    if os.path.isfile(args.coord_files) is False:
         print("ERROR: Could not find file {}. Exiting...".format(args.coord_files))

    elements,geo = xyz_parse(args.coord_files)

    # Generate adjacency table
    adj_mat  = Table_generator(elements,geo)

    # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
    hybridizations = Hybridization_finder(elements,adj_mat)

    # Find atom types
    atom_types = id_types(elements,adj_mat,args.gens)


    # Determine bonding matrix for the compound
    bond_mat = find_lewis(elements, adj_mat, q_tot=args.qtot, verbose=False)

    # Find modes (returns all modes in canonicalized format)
    bonds_tmp,angles_tmp,bond_types_tmp,angle_types_tmp = Find_modes_baonly(adj_mat,atom_types,bond_mat,return_all=1)

    if(args.test):
         test_pt = grid_gen(elements,geo,6,args.size_factor,args.extra)
         return

    grid_points = grid_gen(elements,geo,args.density,args.size_factor,args.density_scale_factor)
    #additional density scale factor, this is by trial and error to get an approximate grid pts that are similar to the final # of charges placed reported by Roux's group
    total_charges=place_all(args.criteria,grid_points,geo,bonds_tmp,hybridizations,elements,args.density,args.size_factor,args.extra,args.file_name,args.folder_name)

    return total_charges

def grid_gen(elements,geo,density,size_factor,density_scale):
 
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
        points.extend(atom.radius_points(density=float(density)*density_scale, scale=size_factor/10, check_against=atoms))

    points = np.array(points) 
    rows = points.shape[0]
    print("Total of {} grid pts generated for Radii scale at {}".format(rows,float(size_factor/10)))

    return points

def place_all(criteria,points,geo,bonds,hybridizations,element,density,size_factor,extra_d_scale,name,folder_name):

    if(os.path.isdir(folder_name)):
         print("ERROR: output folder {}  exist, exiting avoid overwriting....".format(folder_name))
         quit()
    os.makedirs(folder_name) 
    os.chdir(folder_name)
    # place charges along the bonds
    print("Placing charges along the bonds...")
    intersec1 = place_bonds(criteria,points,geo,bonds,element)
    # place charges for lone pairs
    print("Placing charges along lone pairs...")
    intersec2 = place_lone_pair(criteria,points,geo,bonds,hybridizations,element)
    # rm charges placed for lone_pair that are too close to charges placed for bonds
    intersec2 = remove_contact(intersec2,intersec1)
    # place charges in the gap
    intersec3 = place_gap(element,geo,size_factor,extra_d_scale,intersec1,intersec2)


    all_intersec = intersec1 + intersec2 + intersec3
    for count_i,i in enumerate(all_intersec):
         filename = 'pointcharges.{}.pc'.format(count_i)
         with open(filename,'w') as f:
            f.write("1\n")
            f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))

    total_charges = len(all_intersec)
      
    print("Total {} charges placed!".format(total_charges))

    # write original xyz with all the placed charge for visualization
    name = name + '.xyz'
    with open(name,'w') as f:
         f.write("{} \n\n".format(len(geo)+total_charges))
         for count_i,i in enumerate(geo):
            f.write("{} {} {} {}\n".format(element[count_i],i[0],i[1],i[2]))
         for i in intersec1:
            f.write("F {} {} {} \n".format(i[0],i[1],i[2]))
         for i in intersec2:
            f.write("Cl {} {} {} \n".format(i[0],i[1],i[2]))
         for i in intersec3:
            f.write("Br {} {} {} \n".format(i[0],i[1],i[2]))
 
    os.chdir("..")

    return total_charges

def place_bonds(criteria,points,geo,bonds,element):


    # Get bonds vectors
    intersec = []
    for pair in bonds:
         # Get bonds vector and bond central coordinate
         v = np.zeros(3)
         bc = np.zeros(3)
         for i in range(3):
            v[i] = geo[pair[0]][i] - geo[pair[1]][i] 
            bc[i] = (geo[pair[0]][i] + geo[pair[1]][i]) / 2.0 
         # Get the intersection of bondvector and connolly surface
         intersec_tmp = {}
         for p in points: 
            v2 = np.zeros(3)
            for i in range(3):
               v2[i] = p[i] - bc [i]
            angle = angle_between(v,v2)
            if angle < criteria :
               intersec_tmp[angle] = p
         # if more than 2 grid points fits the criteria for one bond, delete the one that's off most
         if len(intersec_tmp) > 2:
            del intersec_tmp[max(intersec_tmp.keys())] 
         elif len(intersec_tmp) == 1:
            print("WARNING: Only 1 charge placed for the bond, try increasing density")
         elif len(intersec_tmp) == 0:
            print("WARNING: 0 charge placed for the bond, try increasing density")
         for key in intersec_tmp:
            intersec.append(intersec_tmp[key])
    # if distance between two pts <1.5 A then remove one
    shape = np.array(intersec).shape
    if len(shape) == 1 and shape[0] == 0 :
      return intersec 
    elif len(shape) == 1 and shape[0] != 0:
      intersec = np.array(intersec).reshape(-1,1)
    Dist_Mat = np.triu(cdist(intersec,intersec))
    rm_ind,y_ind = np.where( (Dist_Mat < 1.5) & (Dist_Mat > 0.0))
    # remove pts that are repeatedly chosen ( because it sneak through Dist_Mat > 0.0)
    tmp_x,tmp_y = np.where( (Dist_Mat == 0.0) )
    for i in range(0,len(tmp_x)):
      if tmp_y[i] > tmp_x[i]:
         rm_ind = np.append(rm_ind,tmp_y[i])
         
    rm_ind = list(dict.fromkeys(rm_ind))
    intersec = [ i for count_i,i in enumerate(intersec) if count_i not in rm_ind] 
      
    print("Placing {} charges along {} bonds....".format(len(intersec),len(bonds)))  


    return intersec

def place_lone_pair(criteria,points,geo,bonds,hybridizations,element):

    # Get bonds vectors
    intersec = []

    # atoms that have special contraints (sp3 O and sp2 N)
    special = {}
    for count_i,hybrid in enumerate(hybridizations): 
         if hybrid == 'sp3' and element[count_i] == 'O':
            link = []
            for pair in bonds:
               if pair[0] == count_i:
                  link.append(pair[1])
               elif pair[1] == count_i:
                  link.append(pair[0])
            if len(link) == 2:
               special[count_i] = {}
               special[count_i]['element'] = 'N'
               special[count_i]['atom_num'] = count_i 
               special[count_i]['link'] = link
         elif hybrid == 'sp2' and element[count_i] == 'N':
            link = []
            for pair in bonds:
               if pair[0] == count_i:
                  link.append(pair[1])
               elif pair[1] == count_i:
                  link.append(pair[0])
            if len(link) == 2:
               special[count_i] = {}
               special[count_i]['element'] = 'N'
               special[count_i]['atom_num'] = count_i 
               special[count_i]['link'] = link

    for atom in special:
         v1 = np.zeros(3)
         v2 = np.zeros(3)
         for i in range(0,3):
            v1[i] = geo[special[atom]['link'][0]][i] - geo[atom][i] 
            v2[i] = geo[special[atom]['link'][1]][i] - geo[atom][i] 

         # A-P directions for P1, P2', P2'', P3', P3'' (see paper)
         P_d = []
         # points (for visual check if the direction is correct)
         P_d_p = []
         # A-P1 's direction (see paper)
         P_d.append(-bisector(v1,v2))
         # A-P2' 's direction
         P_d.append(np.cross(v1,v2))
         # A-P2'' 's directiono
         P_d.append(np.cross(v2,v1))
         # A-P3' 's direction
         P_d.append(bisector(P_d[1],P_d[0]))
         # A-P3 '' 's direction
         P_d.append(bisector(P_d[2],P_d[0]))
         for i in range(0,5):
            P_d_p.append(geo[atom]+P_d[i])
         #intersec_tmp = {}
         intersec_tmp_a  = []
         intersec_tmp_p = []
         count = [0] * 5
         for p in points: 
            v3 = p - geo[atom]
            for i in range(0,5):
               angle = angle_between_normal(P_d[i],v3,normal=True)
               if angle < criteria :
                  intersec_tmp_a.append(angle)
                  intersec_tmp_p.append(p)
                  count[i] += int(1)

         intersec = []
         # if more than 1 grid points fits the criteria for line pair, only keep one that has min angle
         for i in range(0,5):
            begin = int(np.sum(count[:i]))
            end = int(np.sum(count[:i+1])) 
            if count[i] >1:
               index = [intersec_tmp_a.index(min(intersec_tmp_a[begin:end]))]
               #key = [min(list(intersec_tmp.keys())[begin:end])]

            elif count[i] == 0:
               print("WARNING: 0 charge placed for one of the direction, try increasing density")
               index = []
            else:
               index = [begin]

            if index != []:
               intersec.append(intersec_tmp_p[index[0]])
             
    print("Placing {} charges along {} special atoms....".format(len(intersec),len(special.keys())))  

    return intersec

def place_gap(element,geo,size_factor,extra_d_scale,intersec1,intersec2):
   
    density = 6 # this doesn't matter cuz extra_d_scale controls it
    grid_points = grid_gen(element,geo,density,size_factor,extra_d_scale)
    print("Generting initial gap grids... {} pts".format(len(grid_points)))
    grid_points = remove_contact(grid_points,intersec1)
    grid_points = remove_contact(grid_points,intersec2)


    print("Placing {} charges in gap...".format(len(grid_points)))
    return grid_points
    

def remove_contact(grid_points,intersec):
    shape1 = np.array(grid_points).shape
    shape2 = np.array(intersec).shape
    if (len(shape1) ==1 and shape1[0] == 0) or (len(shape2) == 1 and shape2[0] == 0):
      return grid_points
    if len(shape1) == 1 and shape1[0] != 0:
      grid_points = np.array(grid_points).reshape(-1,1)
    if len(shape2) == 1 and shape2[0] != 0:
      intersec = np.array(intersec).reshape(-1,1)
    Dist_Mat = np.triu(cdist(grid_points,intersec))
    rm_ind,y_ind = np.where( (Dist_Mat < 1.5) & (Dist_Mat > 0.0))
    rm_ind = list(dict.fromkeys(rm_ind))
    
    print("{} charges removed due to close contact w/ already placed charge".format(len(rm_ind)))
    grid_points = [ i for count_i,i in enumerate(grid_points) if count_i not in rm_ind]


    return grid_points


def angle_between(v1, v2,normal=False):
# Return acute angle
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))/math.pi*180.0

    if not normal:
       if angle>90.0:
         angle = 180.0 - angle
    return angle

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def bisector(v1,v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    bi = v1_u + v2_u

    return bi

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



