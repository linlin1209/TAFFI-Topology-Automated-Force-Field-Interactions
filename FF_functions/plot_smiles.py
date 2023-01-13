#Nick Iovanac
#Please let me know if you find any bugs or ways to improve this
#niovanac@purdue.edu

#Accepts a file containing SMILES strings and can draw nice PDF images of them. 
# You can adjust the thickness of bonds, the size of non-carbon atom designations, 
# and a few other things. If given multiple SMILES in a file, it will draw them 
# in a grid. It can also display numerical values below the SMILES in case you are 
# looking at properties. Give it a try with the provided example files. 
# Try 'python plot_smiles.py -h' for help on the command line


import matplotlib
matplotlib.use("Agg") #Switch this to "Agg" to run on the compute nodes
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import sys
import argparse
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix,AddHs,GetShortestPath




def main(argv):
    parser = argparse.ArgumentParser(description='Plot compounds on grid along with property values if desired')
    parser.add_argument('-s',dest='smiles',default=None,help='List containing smiles strings. The number must be equal to squared integer values. Required')
    parser.add_argument('-p',dest='props',default=None,help='List containing corresponding property values. Optional')
    parser.add_argument('-pn',dest='prop_name',default='Property',help='Property name (e.g. logP, SAS, pKa). Optional, defaults to a blank')
    parser.add_argument('-o',dest='output_name',default='test.pdf',help='Name to save output under. Default: "test.pdf"')
    parser.add_argument('-f',dest='font_size',default = 8,help = 'Size of non-carbon atom designations. Default: 12')
    parser.add_argument('-w',dest='bond_width',default = 1.0,help ='Thickness with which bonds are drawn. Default: 1.0')
    parser.add_argument('-ps',dest='property_size',default = 12, help='Fontsize for displaying properties. Default: 12')

    args = parser.parse_args()

    smiles = args.smiles
    props = args.props
    prop_name = args.prop_name
    output_name = args.output_name
    font_size = float(args.font_size)
    bond_width = float(args.bond_width)
    property_size = float(args.property_size)

    grid_plot(smiles,props,output_name,prop_name,font_size,bond_width,property_size)


def grid_plot(smiles,props,output_name,prop_name,font_size,bond_width,property_size):

    smile_list =[]
    with open(smiles,'r') as f:
        for line in f:
            smile_list.append(line.strip())

    props_list = []
    if props:
        with open(props,'r') as f:
            for line in f:
                props_list.append(line.strip())



    n_smiles = len(smile_list)
    if not (n_smiles**(1.0/2.0)==int(n_smiles**(1.0/2.0))):
        print('Cant create a grid unless number of SMILES is a perfect square')
        quit()
    root = int(n_smiles**(1.0/2.0))

    #Generate Grid

    fig = plt.figure(figsize=(8,8))

    wspace = 0.5
    grid = gridspec.GridSpec(root,root,wspace=wspace,hspace=0.0)

    #For each box, create a subgrid with two entries. One for the SMILES, one for a property if supplied
    for i in range(n_smiles):
        inner_grid = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec = grid[i],wspace =0.0,hspace = 0.0, height_ratios = [4,1])
        
        #Does not apply any coloration if no properties supplied. Otherwise, will color the grid based on the corresponding
        # values.
        color = 'w'
        if props is not None:
            try:
                t = round(float(props_list[i]),2)
            except ValueError:
                t=props_list[i]
                
        #Plot smile
        ax = plt.Subplot(fig,inner_grid[0],facecolor=color)
        ax=plot_smiles(smile_list[i],ax,color,font_size,bond_width)
        ax.set_yticks([])
        ax.set_xticks([])
        #Plot prop
        fig.add_subplot(ax)
        if props is not None:
            ax = plt.Subplot(fig,inner_grid[1],facecolor=color)
            draw_scale = 1
            ax.text(0.5,0.5,t,color='black',fontweight ='bold',fontsize=property_size,ha='center',va='center',bbox={'facecolor':'none', 'edgecolor':'None','pad':1.0 * draw_scale})
            ax.set_yticks([])
            ax.set_xticks([])
            fig.add_subplot(ax)
            
    #Remove all but the main grids. 


    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
        if ax.is_first_row():
            ax.spines['top'].set_visible(False)
        if ax.is_last_row():
            ax.spines['bottom'].set_visible(False)
        if ax.is_first_col():
            ax.spines['left'].set_visible(False)
        if ax.is_last_col():
            ax.spines['right'].set_visible(False)



    plt.savefig(output_name,bbox_inches='tight')
    return


   


def plot_smiles(smiles,ax,bg,font_size,bond_width):
    """
    Main function for plotting SMILES
    Inputs:

    smiles: A valid SMILES string
    ax: Matplotlib 'axes' object
    bg: Background color
    font_size: How big to draw non carbons
    bond_width: How thicc to draw the bonds

    Output:

    ax: Matplotlib 'axes' object, now with SMILES string plotted
    """

    if not bg:
        bg = 'none'
    box = ax.get_position()

    def get_coords(mol):
        try:
            n_atom = mol.GetNumAtoms()
        except:
            print('Failed to read the following SMILES:')
            print(smiles)
            quit()
        geo_2D = []
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]

        #Note that hydrogens are always kept implicit. Including them seems to mess
        # up the geometry calculation

        AllChem.EmbedMolecule(mol)
        AllChem.Compute2DCoords(mol)
        for c in mol.GetConformers():
            for atom,symbol in enumerate(symbols):
                p = c.GetAtomPosition(atom)
                geo_2D.append([p.x,p.y,p.z])

        geo_2D = np.asarray(geo_2D)
        return geo_2D,symbols


    #You may need to add to this for additional atom types
    color_dict = {
    'N': (0,0,1),
    'O': (1,0,0),
    'Cl': (0,1,0),
    'F': (0,1,1),
    'Br': (1,0,1),
    'S': (1,1,0),
    'OH': (1,0,0)
    }


    #Add '0.2' for each successive ring. Could be done more pythonically
    ring_size_dict = {
        10:1.8,
        9: 1.6,
        8: 1.4,
        7: 1.2,
        6: 1,
        5: 0.8,
        4: 0.6,
        3: 0.4
        }


    other = (0.792,0.792,0.792)

    

    m = Chem.MolFromSmiles(str(smiles))
    geo_2D,symbols = get_coords(m)
    adj_mat = GetAdjacencyMatrix(m,useBO=1)
    
    #Calculates smallest set of smallest rings
    ssr = Chem.GetSymmSSSR(m)

        
        #Draw lines
    line_width = 1.0


    #Kludge for finding alcohols
    alcohol_list = []

    #Loops over adjacency matrix and draws corresponding bonds.
    # Currently has support for single, double, triple bonds, rings,
    # and aromatic rings. NO SUPPORT FOR STEROCHEMISTRY
    for count_i,i in enumerate(adj_mat):
        alcohol_list.append(1)
        if symbols[count_i] == 'O':
            #Kludge for finding alcohols
            if np.sum(adj_mat[count_i]) > 1.0:
                alcohol_list[count_i] = 0
        for count_j,j in enumerate(i):
            if count_j > count_i:
                if j and symbols[count_j] != 'H': #Might want to add something to make H explicit if desired
                    bond_start = np.asarray([geo_2D[count_i][0],geo_2D[count_j][0]])
                    bond_end =np.asarray([geo_2D[count_i][1],geo_2D[count_j][1]])
                    bond_vector = [bond_start[1]-bond_start[0],bond_end[1]-bond_end[0]]

                    ax.add_line(matplotlib.lines.Line2D(bond_start,bond_end,color=(0,0,0),linewidth=bond_width,alpha=1.0))

                    #Draws double and triple bonds
                    
                    if j > 1 and j.is_integer():
                        v_angle = np.arctan2(bond_vector[1],bond_vector[0])
                        orthonormalx = -np.sin(v_angle)*0.2
                        orthonormaly = np.cos(v_angle)*0.2
                        d_start = bond_start+orthonormalx
                        d_end = bond_end + orthonormaly
                        ax.add_line(matplotlib.lines.Line2D(d_start,d_end,color=(0,0,0),linewidth=bond_width,alpha=1.0))
                        if j > 2:
                            t_start = bond_start - orthonormalx
                            t_end = bond_end - orthonormaly
                            ax.add_line(matplotlib.lines.Line2D(t_start,t_end,color=(0,0,0),linewidth=bond_width,alpha=1.0))
                        
        
        #Add aromatic ring designation
    for i in ssr:
        rings = list(i)
        order = adj_mat[rings[0],rings[1]]
        if not order.is_integer():
            aromatics = []
            for a in rings:
                aromatics.append(geo_2D[a])
            ring_size = ring_size_dict[len(aromatics)]
            aromatics = np.asarray(aromatics)
            center = tuple(np.mean(aromatics,axis=0))
            circle = plt.Circle(center,ring_size,color='gray',fill=False)
            ax.add_artist(circle)
            
        


    for count_i,i in enumerate(geo_2D):
        if symbols[count_i] == 'C' or symbols[count_i]=='H':
            continue
        if symbols[count_i] == 'O':
            if alcohol_list[count_i]:
                symbols[count_i] = 'OH'
        
        try:
            color = color_dict[symbols[count_i]]
        except:
            color = other
        ax.text(i[0], i[1],symbols[count_i],style='normal',color=color,fontweight ='bold',fontsize=font_size,ha='center',va='center',bbox={'facecolor':bg, 'edgecolor':'None','boxstyle':'circle,pad=0.1'})

    ax.axis('equal')
    return ax

if __name__ == '__main__':
    main(sys.argv[1:])

