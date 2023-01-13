#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,getpass,json,fnmatch

# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
import codecs,json


def main(argv):
    contour_x = []
    contour_y = []
    contour_z = []
    iteration = []
    z = {}
    title = []
    num_pair = 5
    angledir = [ f for f in os.listdir('.') if (os.path.isdir(f) and f.split('_')[0]=='angle')  ]
    dihedir = [ f for f in os.listdir('.') if (os.path.isdir(f) and f.split('_')[0]=='dihe')  ]
    #for i in angledir:
    for i in dihedir:
      os.chdir(i)
      distdir = [ f for f in os.listdir('.') if os.path.isdir(f) ]
      for j in distdir:
         os.chdir(j)
         with open('fit_charges.log','r') as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0:
                  dist = [ float(fields[count_k+1]) for count_k,k in enumerate(fields) if k == '-w_dist']
                  #angle = [ float(fields[count_k+1]) for count_k,k in enumerate(fields) if k == '-w_angle']
                  angle = [ float(fields[count_k+1]) for count_k,k in enumerate(fields) if k == '-w_dihedral']
                  contour_x += dist
                  contour_y += angle
                  continue
               if len(fields) == 0: continue
               if fields[0] == '<xhi_pot^2>':
                  contour_z .append(float(fields[3]))
                  #if contour_y[-1] > 0.03 or contour_x[-1] > 0.04: 
                  #   contour_x.pop()
                  #   contour_y.pop()
                  #   contour_z.pop()

         os.chdir('..')
      os.chdir('..')
    contour_xhi('dihe_xhi2',contour_x,contour_y,contour_z)
    #contour_xhi('angle_xhi2',contour_x,contour_y,contour_z)
    #contour_vdw(contour_x,contour_y,z,title)
    quit()
    #contour_xhi('eps',contour_x,contour_y,z[1]['eps'])
    #quit()
    x = contour_x
    y = contour_y
    for i in z:
       # contour plot
       fig, ax = plt.subplots(1,2,figsize=(12,5))
       #ax.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
       #cntr2 = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")
       ax[0].tricontour(x, y, z[i]['eps'], linewidths=0.5, colors='k')
       cntr2 = ax[0].tricontourf(x, y, z[i]['eps'], cmap="RdBu_r")
       fig.colorbar(cntr2, ax=ax[0])
       #ax.plot(x, y, 'ko', ms=3)
       ax[0].set(xlim=(min(x), max(x)), ylim=(min(y), max(y)))
       #ax2.set_title('tricontour (%d points)' % npts)
       ax[1].tricontour(x, y, z[i]['sigma'], linewidths=0.5, colors='k')
       cntr2 = ax[1].tricontourf(x, y, z[i]['sigma'], cmap="RdBu_r")
       fig.colorbar(cntr2,ax=ax[1])
       ax[1].set(xlim=(min(x), max(x)), ylim=(min(y), max(y)))
       ax[0].set_ylabel("$\epsilon$ weight")
       ax[0].set_xlabel("$\sigma$ weight")
       ax[1].set_ylabel("$\epsilon$ weight")
       ax[1].set_xlabel("$\sigma$ weight")
       plot_name = "{}_sigma_eps".format(i)
       savefig(plot_name, dpi=300, bbox_inches='tight')
       close(fig)

    contour_xhi('contour_xhi2',contour_x,contour_y,contour_z)
    contour_xhi('iteration',contour_x,contour_y,iteration)
    
    quit()

def contour_xhi(plotname,x,y,z):

    # contour plot
    fig, ax = plt.subplots(1,figsize=(6,5))
    #ax.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
    #cntr2 = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")
    ax.tricontour(x, y, z, linewidths=0.5, colors='k')
    cntr2 = ax.tricontourf(x, y, z, cmap="RdBu_r")
    fig.colorbar(cntr2, ax=ax)
    #ax.plot(x, y, 'ko', ms=3)
    ax.set(xlim=(min(x), max(x)), ylim=(min(y), max(y)))
    #ax2.set_title('tricontour (%d points)' % npts)
    #ax.set_ylabel("angle weight")
    ax.set_ylabel("dihedral weight")
    ax.set_xlabel("dist weight")
    savefig(plotname, dpi=300, bbox_inches='tight')
    close(fig)

    plotname = plotname+'.txt'
    dist,angle,xhi  = zip(*sorted(zip(x, y,z)))
    with open(plotname,'w') as f:
      f.write("{:<10s} {:<10s} {:<10s}\n".format("dist","angle","xhi"))
      for count_i,i in enumerate(dist):
         f.write("{:<10.2f} {:<10.2f} {:<10.5e}\n".format(i,angle[count_i],xhi[count_i]))

    
def contour_vdw(x,y,z,title):

    for count_i,i in enumerate(z):
       # contour plot
       fig, ax = plt.subplots(1,2,figsize=(12,5))
       #ax.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
       #cntr2 = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")
       ax[0].tricontour(x, y, z[i]['eps'], linewidths=0.5, colors='k')
       cntr2 = ax[0].tricontourf(x, y, z[i]['eps'], cmap="RdBu_r")
       fig.colorbar(cntr2, ax=ax[0])
       #ax.plot(x, y, 'ko', ms=3)
       ax[0].set(xlim=(min(x), max(x)), ylim=(min(y), max(y)))
       #ax2.set_title('tricontour (%d points)' % npts)
       ax[1].tricontour(x, y, z[i]['sigma'], linewidths=0.5, colors='k')
       cntr2 = ax[1].tricontourf(x, y, z[i]['sigma'], cmap="RdBu_r")
       fig.colorbar(cntr2,ax=ax[1])
       fig.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
       ax[1].set(xlim=(min(x), max(x)), ylim=(min(y), max(y)))
       ax[0].set_ylabel("$\epsilon$ weight")
       ax[0].set_xlabel("$\sigma$ weight")
       ax[1].set_ylabel("$\epsilon$ weight")
       ax[1].set_xlabel("$\sigma$ weight")
       plot_name = "{}_sigma_eps".format(title[count_i])
       savefig('./'+plot_name, dpi=300, bbox_inches='tight')
       close(fig)
    return


if __name__ == "__main__":
    main(sys.argv[1:])
