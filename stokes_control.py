#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import numpy as np
from scipy.interpolate import interpn

import graphics
import readwrite as rw
import stokes_pressure as pressure

#------------------------------------------------------------------------------
from stokes_solver import run_spLU
#---------------------------------------------------------------------------
# graphics args

# zoom plot
lenx = 4
leny = 2
x_start = 6 
y_start = 0
x_stop= x_start + lenx
y_stop = y_start + leny

# plotting thresholds
vel_max = 5
p_min=60
p_max=110

# log plots
log_linthresh=1e-5
log_cmap_on = False

#---------------------------------------------------------------------------
class Stokes_Solver:
    def __init__(self, Example, args, BC, Re):

        self.Example = Example
        self.args = args
        
        self.BC = BC # BC <- U, Q

        self.Re = Re
        
        # iterative solution args
        self.max_iters = 500000
        self.write_mod = 500
        self.error_mod = 500
        self.err_tol = 1e-8
        

        
#------------------------------------------------------------------------------
    def new_run(self, N):
        ex = self.Example(self.args, self.BC, self.Re, N)
        
        u_init = np.zeros(ex.Nx * ex.Ny)
        v_init = np.zeros(ex.Nx * ex.Ny)
        psi_init = np.zeros(ex.Nx * ex.Ny)
        past_iters = 0
    
        u, v, psi = run_spLU(ex, u_init, v_init, psi_init, self.max_iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
    
        rw.write_stokes(ex, u, v, psi, self.max_iters)
    
    def new_run_many(self, N_0, dN, many):
        self.new_run(N_0)
        N_load = N_0
        for k in range (1, many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N)
            N_load = N

    def load_run_new_many(self, N_0, dN, many):
        N_load = N_0
        self.load_run(N_load)
        for k in range (many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N)
            N_load = N

    def load_run_many(self, N_0, dN, many):
        N = N_0
        for k in range (many): 
            self.load_run(N)
            N *= dN 
                                                                                                                                                                                                                                                                       
    def load_run(self, N):                                
        ex = self.Example(self.args, self.BC, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx*ex.Ny)
        u, v, psi = run_spLU(ex, u, v, psi, self.max_iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
        rw.write_stokes(ex, u, v, psi, self.max_iters+past_iters)
   
    def load_scale(self, N_load, N_scale):
        ex_load = self.Example(self.args, self.BC, self.Re, N_load)
        ex_scale = self.Example(self.args, self.BC, self.Re, N_scale)
        
        points_load = (ex_load.ys, ex_load.xs)
        
        u_load, v_load, psi_load, past_iters = rw.read_stokes(ex_load.filestr+".csv", ex_load.Ny*ex_load.Nx)
        u_load_2D = u_load.reshape((ex_load.Ny,ex_load.Nx))
        v_load_2D = v_load.reshape((ex_load.Ny,ex_load.Nx))
        psi_load_2D = psi_load.reshape((ex_load.Ny,ex_load.Nx))
    
        points_scale = np.meshgrid(ex_scale.ys, ex_scale.xs)
        
        u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale))
        v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale))
        psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale))
    
        u_scaled = u_scaled_2D.ravel(order='F')
        v_scaled = v_scaled_2D.ravel(order='F')
        psi_scaled = psi_scaled_2D.ravel(order='F')
    
        rw.write_stokes(ex_scale, u_scaled, v_scaled, psi_scaled, 0)
    
    def load(self,N):
        ex = self.Example(self.args, self.BC, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx*ex.Ny)
        p = pressure.pressure(ex, u, v)
        p = p.reshape((ex.Ny,ex.Nx))
        u = u.reshape((ex.Ny,ex.Nx))
        v = v.reshape((ex.Ny,ex.Nx))
        
        p = np.flip(p, axis=0)        
        u = np.flip(u, axis=0)
        v = -np.flip(v, axis=0) 
        
        return p, u, v, psi

#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
    def load_plot(self, N,zoom=False):
        ex = self.Example(self.args, self.BC, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)

    
    # Grid domain
        xs = ex.xs
        ys = ex.ys
        
    # Space plot (xi,yj: boundary/interior/exterior)
        # graphics.plot_contour_mesh(ex.space, xs, ys, 'space',['space', 'x', 'y'],log_cmap=False)

    # Zoom domain
        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(xs, ys, ex, x_start, x_stop, y_start, y_stop)

    # Pressure plot: 
        p = pressure.pressure(ex, u, v)
        dp = pressure.dp(ex, p) 
        
        p_2D = p.reshape((ex.Ny,ex.Nx))
        dp_str = ', $\Delta p =%.2f$'%(dp)

    
        ax_labels_p = ['$p$', '$x$', '$y$']
        title_p = 'Stokes\n' + ex.spacestr + dp_str
    
        p_ma = np.ma.masked_where(ex.space==-1, p_2D)
        p_ma = np.flip(p_ma, axis=0)
       
        if zoom:
            p_zoom = graphics.grid_zoom_2D(p_ma, ex, x_start, x_stop, y_start, y_stop)  
            graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, title_p, ax_labels_p, vmin=p_min, vmax=p_max, log_cmap=log_cmap_on, n_contours=100)#, y_lim = min(ex.y_peaks))
        else:
            graphics.plot_contour_mesh(p_ma, xs, ys, title_p, ax_labels_p,  vmin=p_min, vmax=p_max, log_cmap=log_cmap_on, n_contours=200)#, y_lim = np.min(ex.y_peaks))
    
    #  Velocity plot: 
    
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        title = 'Stokes\n' + ex.spacestr + dp_str
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        u_2D = u.reshape((ex.Ny,ex.Nx))
        v_2D = v.reshape((ex.Ny,ex.Nx)) 

        uv_mag = np.sqrt(u_2D**2 + v_2D**2)
        
        uv_mag = np.flip(uv_mag, axis=0)
        u_2D = np.flip(u_2D, axis=0)
        v_2D = -np.flip(v_2D, axis=0)
        
        if zoom:
            u_2D_zoom = graphics.grid_zoom_2D(u_2D, ex, x_start, x_stop, y_start, y_stop)
            v_2D_zoom = graphics.grid_zoom_2D(v_2D, ex, x_start, x_stop, y_start, y_stop)
            uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, ex, x_start, x_stop, y_start, y_stop)
            graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, title, ax_labels, vmin=0, vmax=vel_max, log_cmap=log_cmap_on)
        else:
            
            graphics.plot_stream_heat(u_2D, v_2D, xs, ys, uv_mag, title, ax_labels,  vmin=0, vmax=vel_max, log_cmap=log_cmap_on) 
        
    # Stream plot:
    
        # ax_labels = ['$\psi(x,y)$', '$x$', '$y$']
        # title = 'Stream $\psi(x,y)$ \n' + ex.spacestr + dp_str
        # stream_2D = psi.reshape((ex.Ny,ex.Nx))
        # stream_2D_ma = np.ma.masked_where(ex.space==-1, stream_2D)
        # stream_2D_ma = np.flip(stream_2D_ma, axis=0)
        # if zoom:
        #     stream_2D_zoom = graphics.grid_zoom_2D(stream_2D_ma, ex, x_start, x_stop, y_start, y_stop)
        #     graphics.plot_contour_mesh(stream_2D_zoom, xs_zoom, ys_zoom, title, ax_labels, log_cmap=log_cmap_on, n_contours=20, vmin=None, vmax=ex.BC.Q)
        # else:
        #     graphics.plot_contour_mesh(stream_2D_ma, xs, ys, title, ax_labels, log_cmap=log_cmap_on, n_contours=20, vmin=None, vmax=ex.BC.Q)
        return dp
