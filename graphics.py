#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:53:47 2022

Graphics helpers
@author: sarahdennis
"""
import numpy as np
from matplotlib import pyplot as pp
from matplotlib import colors
from matplotlib import patches

#------------------------------------------------------------------------------
# COLOURINGS
#------------------------------------------------------------------------------

# -------- 2D field ----------------velocity-----------------------------------

stream_cmap = pp.cm.Spectral_r(np.arange(pp.cm.Spectral_r.N))
colour_map_stream = colors.ListedColormap(stream_cmap)

# -------- 2D value ------------pressure/incompressiblity/error----------------

mesh_cmap = pp.cm.RdYlBu_r(np.arange(pp.cm.RdYlBu_r.N))
colour_map_mesh = colors.ListedColormap(mesh_cmap)

#---------LEGEND---------------------------------------------------------------

# colour_bar_scale=0.015    # for very long figures, H=1.25, L=4
colour_bar_scale=0.022      # for long figures like H=2, L=4
# colour_bar_scale=0.5      # for almost square figures like H=2.75, L=4

colour_bar_pad = 0.02

#------------------------------------------------------------------------------
# RESOLUTION
#------------------------------------------------------------------------------
dpi=1200

n_contours = 20
contour_width = 0.25
stream_width = 1
line_width = 1.5

linthresh = 1e-8

# font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 16

pp.rc('font', size=SMALL_SIZE)          # controls default text sizes
pp.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
pp.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
pp.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pp.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pp.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
pp.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#------------------------------------------------------------------------------
# LINE PLOTS
#------------------------------------------------------------------------------
def plot_2D(fs, xs, title, axis_labels, color='darkmagenta', marker='o'):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi

    pp.plot(xs, fs, color=color, linewidth=.8, marker=marker)

    pp.title(title, fontweight="bold")
    
    pp.xlabel(axis_labels[0])
    
    pp.ylabel(axis_labels[1])
    # pp.ylim(0, 1.25*max(fs))
    pp.minorticks_on()
   
    return fig

def plot_2D_multi(fs, xs, title, fun_labels, ax_labels, loc='upper', colors='pri'):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi
    ax = fig.add_subplot()
    if colors== 'pri':
        cs = ['r','b', 'forestgreen','darkmagenta', 'darkorange']
    else: 
        cs=['forestgreen', 'darkmagenta', 'darkorgange']
    markers = ['D', 'o', 's', '*', 'H', 'X']
    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=fun_labels[i], color=cs[i], linewidth=0.8, marker=markers[i])
    
    # ax.set_xlim([None, 18])
    # ax.set_ylim([0.25*np.min(fs),min(1.1*np.max(fs),25)])

    # ax.set_ylim([2,min(1.1*np.max(fs),5)])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    pp.title(title,  fontweight ="bold")
    pp.minorticks_on()
    
    if loc== 'upper':
        fig.legend(bbox_to_anchor=(0.9, 0.875))
    elif loc=='left':
        fig.legend(bbox_to_anchor=(0.325, 0.875))
    elif loc=='lower':
        
        fig.legend(bbox_to_anchor=(0.325, 0.45))
    elif loc=='center':
        fig.legend(bbox_to_anchor=(0.5, 0.875))
    else: 
        fig.legend(bbox_to_anchor=(0.9, 0.375))
    
    return fig

#------------------------------------------------------------------------------   
# LOG LINE PLOTS
#------------------------------------------------------------------------------

def plot_log(fs, xs, title, ax_labels):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi
    pp.loglog(xs, fs, color='b')   

    pp.title(title)
    
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    
    return fig


def plot_log_multi(fs, xs, title, f_labels, ax_labels, linthresh=linthresh, bigO_on=False, loc='left', log_x=False, log_y=True, colors='pri'):
    pp.rcParams['figure.dpi'] = dpi
    fig = pp.figure()
    
    
    ax = fig.add_subplot() 
    if colors== 'pri':
         cs = ['r','b', 'forestgreen','darkmagenta', 'darkorange']
    else: 
         cs=['forestgreen', 'darkmagenta', 'darkorgange']
    markers = ['D', 'o', 's', '*', 'H', 'X']

    pp.rcParams["lines.linewidth"] =line_width
    for i in range(len(fs)):
        ax.plot(xs, fs[i], label=f_labels[i], color=cs[i], linewidth=0.8, marker=markers[i], markevery=1)
    
    
    # reference lines
    if bigO_on:
        O1 = 1
        O2 = 1
        ax.plot(xs, [O1*x**-1 for x in xs], label="$\mathcal{O}(\Delta x)$", color='darkgrey')    
        ax.plot(xs, [O2*x**-2 for x in xs], label="$\mathcal{O}(\Delta x^2)$", color='k')
    
    if log_x:
        ax.set_xscale('log')
        
    if log_y:
        ax.set_yscale('log')
        # ax.set_yscale('symlog', linthresh=linthresh)

    # ax.set_ylim(min(1,0.5*np.min(fs)),2*np.max(fs))
    # ax.set_ylim(1,2*np.max(fs))

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    
    # ax.minorticks_on()
    # ax.yaxis.set_minor_locator(MultipleLocator(2.5))
    
    pp.title(title)
    
    if loc== 'upper': #upper-right
        fig.legend(bbox_to_anchor=(0.9, 0.875))
    elif loc=='left': #upper-left
        fig.legend(bbox_to_anchor=(0.3, 0.89))
    elif loc=='lower': #lower-right
        fig.legend(bbox_to_anchor=(0.3, 0.275))
    elif loc=='center': #lower-right
        fig.legend(bbox_to_anchor=(0.5, 0.875))
    else: #lower left
        fig.legend(bbox_to_anchor=(0.9, 0.275))  
        
    return fig

#------------------------------------------------------------------------------
# STREAMLINE PLOT <U(X,Y),V(X,Y)>
#------------------------------------------------------------------------------
        
def plot_stream_heat(vx, vy, xs, ys, color_map, title, ax_labels, vmin, vmax, log_cmap=False, linthresh=linthresh):

    pp.rcParams['figure.dpi'] = dpi
    
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[2,1]
    # stream_density=[xs.shape[0]/(2*ys.shape[0]),1]
  
    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=vmax, clip=False)
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=stream_width, color=color_map, cmap=colour_map_stream, norm=norm_symLog)
    else:
        no_norm = colors.CenteredNorm(vcenter=vmin + vmax/2, halfrange=vmax/2, clip=True)
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=stream_width, color=color_map, cmap=colour_map_stream, norm=no_norm)

    cb=pp.colorbar(stream_plot.lines, label=ax_labels[0], fraction=colour_bar_scale, pad=colour_bar_pad)
    ticks = np.linspace(vmin, vmax, num=5)
    cb.set_ticks(ticks)

    # remove contours arrows
    ax = pp.gca()
    for art in ax.get_children():
        if not isinstance(art, patches.FancyArrowPatch):
            continue
        art.remove()        


    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])
    ax.set_aspect('equal')
    pp.minorticks_on()
    pp.show()
    
#------------------------------------------------------------------------------       
# VALUE PLOT Z = F(X,Y)
#------------------------------------------------------------------------------

def plot_contour_mesh(zs, xs, ys, title, labels, vmin, vmax, vscale=None,log_cmap=False, linthresh=linthresh, n_contours=n_contours):
    pp.rcParams['figure.dpi'] = dpi
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    if vscale is not None:
        zs/= vscale
        vmin/= vscale
        vmax/= vscale
       
    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=0)#vmax, clip=False)
        color_plot = pp.pcolor(X, Y, zs, cmap=colour_map_mesh, norm=norm_symLog)
    else:
        color_plot = pp.pcolor(X, Y, zs, cmap=colour_map_mesh, vmin=vmin, vmax=vmax)
    
    cb = pp.colorbar(color_plot, label=labels[0], fraction=colour_bar_scale, pad=colour_bar_pad)
    
    if vmin is not None and vmax is not None:
        ticks = np.linspace(vmin, vmax, num=5)
        cb.set_ticks(ticks)

    pp.rcParams["lines.linewidth"] = contour_width
    pp.contour(X, Y, zs, n_contours, colors='black', negative_linestyles='solid')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])

    pp.minorticks_on()
    ax = pp.gca()
    ax.set_aspect('equal')    
    # ax.set_ylim(y_lim)

    pp.show()    

#------------------------------------------------------------------------------------
def grid_zoom_2D(grid, ex, x_start, x_stop, y_start, y_stop):
    i_start = int((x_start - ex.x0)/ex.dx)
    i_stop = int((x_stop - ex.x0)/ex.dx)
    j_start = int((y_start - ex.y0)/ex.dy)
    j_stop = int((y_stop - ex.y0)/ex.dy)
    i_max = int((ex.xf - ex.x0)/ex.dx)
    j_max = int((ex.yf - ex.y0)/ex.dy)
    if i_start < 0 or j_start < 0 or i_stop > i_max or j_stop > j_max or i_start > i_stop or j_start > j_stop:
        raise Exception('graphics zoom window out of bounds')
    return grid[j_start:j_stop,i_start:i_stop]

def grid_zoom_1D(grid_x, grid_y, ex, x_start, x_stop, y_start, y_stop):
    
    i_start = int((x_start - ex.x0)/ex.dx)
    i_stop = int((x_stop - ex.x0)/ex.dx)
    j_start = int((y_start - ex.y0)/ex.dy)
    j_stop = int((y_stop - ex.y0)/ex.dy)
    i_max = int((ex.xf - ex.x0)/ex.dx)
    j_max = int((ex.yf - ex.y0)/ex.dy)
    if i_start < 0 or j_start < 0 or i_stop > i_max or j_stop > j_max or i_start > i_stop or j_start > j_stop:
        raise Exception('graphics zoom window out of bounds')
    return grid_x[i_start:i_stop], grid_y[j_start:j_stop]
        
