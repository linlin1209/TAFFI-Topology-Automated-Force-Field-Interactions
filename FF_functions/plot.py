import sys,argparse,os,re,fnmatch
from scipy.optimize import minimize, rosen, rosen_der
from scipy.optimize import OptimizeResult,approx_fprime
from matplotlib import pyplot as plt
import numpy as np
from copy import deepcopy
import random


    

def plot_convergence(x_fit):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []

    x_vals = np.arange(-5,5.1,0.1)
    y_vals = x_vals**2.0

    y_fit = [ i**2.0 for i in x_fit]


    plot_handle, = ax.plot(x_vals,y_vals,ls="-",linewidth=4.0,c=color_list[0]) 
    plot_handel, = ax.plot(x_fit,y_fit,marker='.',ls='--',markersize=20,color=color_list[1],markeredgewidth=0.0)
    # Set limits based on the largest range

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("f(x)",fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("x",fontsize=32,labelpad=10,fontweight='bold')
    plot_name = "toy.png"

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return

def plot_custom(x,y,name='toy_2',folder='.'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    for count_i,i in enumerate(x.keys()):
       plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)
       #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[0])

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("xhi2",fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("iteration",fontsize=32,labelpad=10,fontweight='bold')
    plot_name = folder+'/'+name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close(fig)


    # Save data to file
    dataname = folder+'/'+name+'.txt'
    niter = len(x[list(x.keys())[0]])
    #with open(dataname,'w') as  f:
    #     f.write('{:<15s} {}\n'.format('iter',' '.join([ "{:<15s}".format(i) for i in x ])))
    #     for i in range(niter):
    #        f.write('{:<15d} {}\n'.format(i,' '.join([ "{:< 15.8f}".format(y[j][i]) for j in x ])))
    return

def plot_custom_2(x,y,name='toy_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    for count_i,i in enumerate(x):
       plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("y",fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("x",fontsize=32,labelpad=10,fontweight='bold')
    plot_name = name+'.png'

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return
def plot_custom_3(x,y1,y2,y3,name='toy_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    plot_handel, = ax.plot(x,y1,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[0],label='DFT')
    plot_handel, = ax.plot(x,y2,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[1],label='TAFFI')
    plot_handel, = ax.plot(x,y3,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[2],label='Drude-TAFFI')
    #for count_i,i in enumerate(x):
    #   plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max

    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel(r'$\mu$$\mathrm{_{calc}}$',fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel(r'$\mu$$\mathrm{_{exp}}$',fontsize=32,labelpad=10,fontweight='bold')
    plot_name = name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False,fontsize=20)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5),fontsize=20)

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return

def plot_custom_5(x,y1,y2,name='toy_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    plot_handel, = ax.plot(x,y1,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[0],label='DFT')
    plot_handel, = ax.plot(x,y2,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[1],label='Drude-TAFFI')
    #for count_i,i in enumerate(x):
    #   plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max

    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel(r'$\alpha$$\mathrm{_{calc}(A^3)}$',fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel(r'$\alpha$$\mathrm{_{exp}(A^3)}$',fontsize=32,labelpad=10,fontweight='bold')
    plot_name = name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False,fontsize=20)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5),fontsize=20)

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return

def plot_custom_4(x,y1,name='toy_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    x1=[5.11]
    y2=[19.02]
    plot_handel, = ax.plot(x,y1,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[0],label='TAFFI')
    plot_handel, = ax.plot(x1,y2,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[1],label='Drude-TAFFI')
    #for count_i,i in enumerate(x):
    #   plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max

    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel(r'$\epsilon$$\mathrm{_{calc}}$',fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel(r'$\epsilon$$\mathrm{_{exp}}$',fontsize=32,labelpad=10,fontweight='bold')
    plot_name = name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False,fontsize=32)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5),fontsize=20)

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return

def plot_custom_6(x,y1,y2,y3,y4,name='toy_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    plot_handel, = ax.plot(x,y1,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[0],label='DFT')
    plot_handel, = ax.plot(x,y2,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[1],label='TAFFI')
    plot_handel, = ax.plot(x,y3,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[2],label='TAFFI(avg500ps)')
    plot_handel, = ax.plot(x,y4,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[3],label='Drude-TAFFI')
    #for count_i,i in enumerate(x):
    #   plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max

    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel(r'$\mu$$\mathrm{_{calc}}$(debye)',fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel(r'$\mu$$\mathrm{_{exp}}$(debye)',fontsize=32,labelpad=10,fontweight='bold')
    plot_name = name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False,fontsize=20)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5),fontsize=20)

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return
def plot_custom_7(x,y,name='toy_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    # unit:  atom -> mol        mol->g  A^3->cm^3
    Const = (1/(6.02*10**23)) * 98.96 * (10**24)
    N_mol = 625 
    y = [ i*2 for i in y]
    y2 = [ N_mol/(i**3)*Const for i in y] 
    x2 = [ int(i/1000) for i in x]
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    plot_handle, = ax.plot(x2,y,ls="-",linewidth=4.0,c=color_list[0],label='box length') 
    ax2 = ax.twinx()
    plot_handle, = ax2.plot(x2,y2,ls="-",linewidth=4.0,c=color_list[1],label='density') 
    #for count_i,i in enumerate(x):
    #   plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max

    #ax.set_xlim([x_min,x_max])
    #ax.set_ylim([y_min,y_max])


    # Set Labels and Save Name
    ax.set_xlabel('time(ps)',fontsize=32,labelpad=10,fontweight='bold')
    ax.set_ylabel('box length(A)',fontsize=32,labelpad=10,fontweight='bold')
    ax2.set_ylabel(r'$\rho$(g/cm^3)',fontsize=32,labelpad=10,fontweight='bold')
    plot_name = name+'.png'

    # Generate Legend
    ax.legend(loc=(0.45,0.5),frameon=False,fontsize=20)
    ax2.legend(loc=(0.45,0.4),frameon=False,fontsize=20)
    #handles, labels = ax.get_legend_handles_labels()
    #lgd = ax.legend(handles, labels, loc='center left',fontsize=15)
    #handles, labels = ax2.get_legend_handles_labels()
    #lgd = ax2.legend(handles, labels, loc='center left',fontsize=15)

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=20,pad=10,direction='out',width=4,length=4,color=color_list[0])
    ax.tick_params(axis='both', which='minor',labelsize=20,pad=10,direction='out',width=4,length=3)
    ax2.tick_params(axis='both', which='major',labelsize=20,pad=10,direction='out',width=4,length=4,color=color_list[1])
    ax2.tick_params(axis='both', which='minor',labelsize=20,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return 

def plot_custom_8(exp_data,result_data,lmp_data,opls_data,key,name='TOY_2'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    x = [ exp_data[i][key] for i in exp_data if exp_data[i]['T'] == 298.15 or exp_data[i]['T'] == 285.45] 
    y1 = [ result_data[i][key] for i in result_data if result_data[i]['T'] == 298.15 or result_data[i]['T'] == 285.45] 
    y2 = [ lmp_data[i][key] for i in lmp_data if lmp_data[i]['T'] == 298.15 or lmp_data[i]['T'] == 285.45] 
    y3 = [ opls_data[i][key] for i in opls_data if opls_data[i]['T'] == 298.15 or opls_data[i]['T'] == 285.45] 
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    plot_handel, = ax.plot(x,y1,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[0],label='Spoel setting')
    plot_handel, = ax.plot(x,y2,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[1],label='original lammps')
    plot_handel, = ax.plot(x,y3,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[2],label='opls')
    #plot_handel, = ax.plot(x,y4,marker='.',linestyle='None',markeredgewidth=0.0,markersize=20,color=color_list[3],label='Drude-TAFFI')
    #for count_i,i in enumerate(x):
    #   plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    if x_min < y_min: y_min = x_min;
    else: x_min = y_min
    if x_max > y_max: y_max = x_max;
    else: x_max = y_max

    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Plot the diagonal line
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0")
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("{}$\mathrm{{_{{calc}}}}$".format(key),fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("{}$\mathrm{{_{{exp}}}}$".format(key),fontsize=32,labelpad=10,fontweight='bold')
    plot_name = key+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False,fontsize=20)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5),fontsize=20)

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_inches='tight')
    plt.close(fig)

    return
def read_data(name):
   exp_dipole = []
   DFT_dipole = []
   TAFFI_dipole = []
   Drude_dipole = []
   exp_polar = []
   DFT_polar = []
   Drude_polar = [] 
   with open(name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            DFT_polar.append(float(fields[2]))
            DFT_dipole.append(float(fields[3]))
            Drude_polar.append((fields[4]))
            Drude_dipole.append(float(fields[5]))
            TAFFI_dipole.append(float(fields[6]))
            exp_polar.append(float(fields[7]))
            exp_dipole.append(float(fields[8]))
   return exp_dipole,DFT_dipole,TAFFI_dipole,Drude_dipole,exp_polar,DFT_polar,Drude_polar

def read_data_2(name):
   exp_eps = []
   TAFFI_eps = []
   with open(name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            TAFFI_eps.append(float(fields[2]))
            exp_eps.append(float(fields[3]))
   return TAFFI_eps,exp_eps

def read_data_3(name):
   exp_alpha = []
   DFT_alpha = []
   Drude_alpha = []
   with open(name,'r') as f:
        for lc,lines in enumerate(f):
            if lc == 0: continue
            fields = lines.split()
            DFT_alpha.append(float(fields[2]))
            Drude_alpha.append(float(fields[5]))
            exp_alpha.append(float(fields[8]))
   return DFT_alpha,Drude_alpha,exp_alpha
def read_data_4(name):
   exp_dipole = []
   DFT_dipole = []
   TAFFI_dipole = []
   TAFFI_avg_dipole = []
   Drude_dipole = []
   with open(name,'r') as f:
        for lc,lines in enumerate(f):
            if lc == 0: continue
            fields = lines.split()
            exp_dipole.append(float(fields[1]))
            DFT_dipole.append(float(fields[2]))
            TAFFI_dipole.append(float(fields[3]))
            TAFFI_avg_dipole.append(float(fields[4]))
            Drude_dipole.append(float(fields[5]))
   return exp_dipole,DFT_dipole,TAFFI_dipole,TAFFI_avg_dipole,Drude_dipole
def read_data_5(name):
   time = []
   box = []
   with open(name,'r') as f:
      for lc,lines in enumerate(f):
            fields = lines.split()
            time.append(float(fields[0]))
            box.append(float(fields[1]))
            #if fields[0] == '510000':break
   return time,box

def read_data_6(name):
   data = {}
   with open(name,'r') as f:
      for lc,lines in enumerate(f):
            if lc == 0: continue # skip title
            fields = lines.split()
            data[fields[0]] = {}
            data[fields[0]]['T'] = float(fields[1])
            data[fields[0]]['rho'] = float(fields[2])
            data[fields[0]]['H_vap'] = float(fields[4])
            data[fields[0]]['eps'] = float(fields[6])
            data[fields[0]]['iso'] = float(fields[8])
            data[fields[0]]['alpha'] = float(fields[10])
            #data[fields[0]]['tension'] = float(fields[12])
   return data


def main(argv):
   exp_data = read_data_6('exp_result.txt') 
   lmp_data = read_data_6('lammps_result.txt') 
   opls_data = read_data_6('opls_result.txt')
   result_data = read_data_6('cale_result.txt')
   for i in exp_data['methanol']:
      if i == 'T':continue
      plot_custom_8(exp_data,result_data,lmp_data,opls_data,i)
   quit()
   #exp_dipole,DFT_dipole,TAFFI_dipole,Drude_dipole,exp_polar,DFT_polar,Drude_polar = read_data('data.txt')
   #TAFFI_eps,exp_eps = read_data_2('eps.txt')
   #DFT_alpha,Drude_alpha,exp_alpha = read_data_3('Polarizability.txt')
   #plot_custom_5(exp_alpha,DFT_alpha,Drude_alpha,name='Polarizability')
   #exp_dipole,DFT_dipole,TAFFI_dipole,TAFFI_avg_dipole,Drude_dipole = read_data_4('dipole.txt')
   #plot_custom_6(exp_dipole,DFT_dipole,TAFFI_dipole,TAFFI_avg_dipole,Drude_dipole,name='Dipole')
   #plot_custom_3(exp_dipole,DFT_dipole,TAFFI_dipole,Drude_dipole,name='dipole')
   #plot_custom_4(exp_eps,TAFFI_eps,name='eps')

if __name__ == "__main__":
   main(sys.argv[1:])
