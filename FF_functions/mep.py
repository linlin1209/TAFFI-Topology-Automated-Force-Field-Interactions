#!/bin/env python
import os
import sys
import subprocess
import numpy as np


def read_xyz(xyz):
    atoms = list()
    x = list()
    y = list()
    z = list()

    with open(xyz) as fp:
        # Skip the first two lines.
        next(fp)
        next(fp)
        for line in fp:
            data = line.split()
            atoms.append(data[0])
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))

    return atoms, np.array(x), np.array(y), np.array(z)


def read_vpot(vpot):
    v = list()

    with open(vpot) as fp:
        next(fp)
        for line in fp:
            data = line.split()
            v.append(float(data[3]))

    return np.array(v)

def main(argv):


    ang_to_au = 1.0 / 0.5291772083

    elements = [None,
         'H', 'He',
         'Li', 'Be',
         'B', 'C', 'N', 'O', 'F', 'Ne',
         'Na', 'Mg',
         'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
         'K', 'Ca',
         'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
         'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
         'Rb', 'Sr',
         'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
         'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
         'Cs', 'Ba',
         'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
         'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
         'Fr', 'Ra',
         'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
         'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub']

    basename = argv[0]
    xyz = basename + '.xyz'

    if not os.path.isfile(xyz):
        sys.exit('Could not find the .xyz. To quickly generate one for '
                 'your molecule run: echo 11 | orca_plot {}.gbw -i.'.format(basename))

    atoms, x, y, z = read_xyz(xyz)


    natoms = len(atoms)

    extent = 7.0
    xmin = x.min() * ang_to_au - extent
    xmax = x.max() * ang_to_au + extent
    ymin = y.min() * ang_to_au - extent
    ymax = y.max() * ang_to_au + extent
    zmin = z.min() * ang_to_au - extent
    zmax = z.max() * ang_to_au + extent


    vpot = read_vpot('result.out')

    with open(basename + '_mep.cube', 'w') as fp:
        fp.write('Generated with ORCA\n')
        fp.write('Electrostatic potential for ' + basename + '\n')
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            len(atoms), xmin, ymin, zmin))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, (xmax - xmin) / float(npoints - 1), 0.0, 0.0))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, 0.0, (ymax - ymin) / float(npoints - 1), 0.0))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, 0.0, 0.0, (zmax - zmin) / float(npoints - 1)))
        for i, atom in enumerate(atoms):
            index = elements.index(atom)
            fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n'.format(
                index, 0.0, x[i] * ang_to_au, y[i] * ang_to_au, z[i] * ang_to_au))

        m = 0
        n = 0
        vpot = np.reshape(vpot, (npoints, npoints, npoints))
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    fp.write('{0:14.5e}'.format(vpot[ix][iy][iz]))
                    m += 1
                    n += 1
                    if (n > 5):
                        fp.write('\n')
                        n = 0
                if n != 0:
                    fp.write('\n')
                    n = 0
if __name__ == '__main__':
    main(sys.argv[1:])


