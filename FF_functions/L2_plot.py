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
    for i in range(1,num_pair):
      z[i] = {}
      z[i]['eps'] = []
      z[i]['sigma'] = []
    with open('extract_vdw.log','r') as f:
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) == 0 and flag == 0: continue
         
            if flag == 0  and fields[0] == 'Parameter' and fields[1] == 'update':
               flag = 1
               readline = 0
               continue

            if flag == 0  and fields[0] == 'L2' and fields[1] == 'sigma' :
               contour_x.append(float(fields[3]))
               continue

            if flag == 0  and fields[0] == 'L2' and fields[1] == 'epsilon' :
               contour_y.append(float(fields[3]))
               continue

            if flag == 0  and fields[0] == 'iteration:' :
               iteration.append(float(fields[1]))
               continue
   
            if flag == 0  and fields[0] == 'Final' and fields[1] == 'xhi^2' :
               #if fields[3] == 'nan' :
               if fields[3] == 'nan' or float(fields[3])>0.6:
                  contour_x.pop()
                  contour_y.pop()
                  iteration.pop()
                  continue
               contour_z.append(float(fields[3]))
               continue
           
            if flag == 1 :
               if readline == 0:
                  readline += 1
                  continue
               if readline < num_pair:
                  if fields[-3] == 'nan':
                     readline += 1
                     continue
                  z[readline]['eps'].append(float(fields[-3])) 
                  z[readline]['sigma'].append(float(fields[-1]))
                  title.append(fields[1]+' '+fields[2])
                  readline += 1
                  continue
               if readline == num_pair:
                  flag = 0
                  continue
    #print(len(z[1]['eps']))
    #print(len(contour_x))
    contour_xhi('contour_xhi2_06',contour_x,contour_y,contour_z)
    contour_xhi('iteration',contour_x,contour_y,iteration)
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
    ax.set_ylabel("$\epsilon$ weight")
    ax.set_xlabel("$\sigma$ weight")
    savefig(plotname, dpi=300, bbox_inches='tight')
    close(fig)
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
