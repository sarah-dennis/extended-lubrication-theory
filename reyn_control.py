# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics

import domain as dm
import reyn_velocity as rv
import reyn_pressure as rp 
import reyn_perturbed as rpert

# zoom plot
lenx = 4
leny = 2
x_start = 6
y_start = 0
x_stop= x_start + lenx
y_stop = y_start + leny

# colorbar min max
vel_max = 5
p_min= 0
p_max = 40

# log plots
log_linthresh=1e-8  
log_cmap_on = False
        
class Reynolds_Solver: 
    def __init__(self, Example, BC, args=None):
        self.Example = Example #initialize height = Example(args) in the solver
        self.args = args
        
        self.BC = BC        


    def fd_solve(self, N, plot=True, scaled=False, zoom=False):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)

        reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)

        velocity = rv.Reyn_Velocity(height, self.BC, ps=reyn_pressure.ps_1D)
        if plot: 
            self.p_plot(height, reyn_pressure, velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, velocity, reyn_pressure, solver_title, scaled, zoom)

        return reyn_pressure, velocity

    def fd_VA_ELT_solve(self, N, plot=True, scaled=False, zoom=False):
        solver_title = "VA-ELT"
        
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        adj_pressure = rp.VA_ELT_Pressure(height, self.BC)

        adj_velocity = rv.VA_ELT_Velocity(height,self.BC, adj_pressure)

        if plot:
            self.p_plot(height, adj_pressure , adj_velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, adj_velocity, adj_pressure, solver_title, scaled, zoom)
       
        return adj_pressure, adj_velocity
    
    def fd_TG_ELT_solve(self, N, plot=True, scaled=False, zoom=False):
        solver_title = "T.G.-ELT"

        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        adj_pressure = rp.TG_ELT_Pressure(height, self.BC)
   
        adj_velocity = rv.Reyn_Velocity(height, self.BC, ps=adj_pressure.ps_1D)

        
        if plot:
            self.p_plot(height, adj_pressure , adj_velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, adj_velocity, adj_pressure, solver_title, scaled, zoom)
       
        return  adj_pressure, adj_velocity

    def fd_pert_solve(self, N, order,  plot=True, scaled=False, zoom=False, get_all = False):
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)

    
        reyn_velocity = rv.Reyn_Velocity(height, self.BC, reyn_pressure.ps_1D)                   
        
        pert = rpert.PerturbedReynSol(height, self.BC, order, reyn_pressure, reyn_velocity)

        
        if plot:
    
            if order > 3:
                solver_title ="$\epsilon^4$ PLT"
                self.p_plot(height, pert.pert4_pressure, pert.pert4_velocity.Q, solver_title, scaled, zoom)
                self.v_plot(self.BC, height, pert.pert4_velocity, pert.pert4_pressure, solver_title, scaled, zoom)
        
            # elif order > 1:
            solver_title = "$\epsilon^2$ PLT"
            self.p_plot(height, pert.pert2_pressure, pert.pert2_velocity.Q, solver_title, scaled, zoom)
            self.v_plot(self.BC, height, pert.pert2_velocity, pert.pert2_pressure, solver_title, scaled, zoom)
           
            
        if get_all:
            return pert
        else:
            if order < 3:
                return pert.pert2_pressure, pert.pert2_velocity
            
            else:
                return pert.pert4_pressure, pert.pert4_velocity

    
    def p_plot(self, height, pressure, flux, solver_title, scaled=False, zoom=False):

           
        # graphics.plot_2D(pressure.ps_1D, height.xs, p_title, p_labels)
        x_scale = height.xs[-1]-height.xs[0]
        y_scale = max(height.hs)
        u_scale = flux/y_scale
        p_scale = flux*x_scale/y_scale #*visc
        
        
        ps_2D = np.nan_to_num(pressure.ps_2D)

        dp_nondim = (sum(ps_2D[:,0])/height.hs[0] - sum(ps_2D[:,-1])/height.hs[-1])*height.dy/p_scale/y_scale
        
        dp_dim = (sum(ps_2D[:,0])/height.hs[0] - sum(ps_2D[:,-1])/height.hs[-1])*height.dy
        
        
        if scaled:
            x_scale = height.xs[-1]-height.xs[0]
            y_scale = max(height.hs)
            u_scale = flux/y_scale
            p_scale = flux*x_scale/y_scale #*visc
            
            paramstr = "$Q=%.2f$, $U=%.2f$, $\Delta P=%.2f$"%(1,self.BC.U/u_scale, dp_nondim)
            p_title = solver_title + '\n' + paramstr
            p_labels =  ['$P$','$X$', '$Y$']
            
            if zoom: 
                xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
                p_zoom = graphics.grid_zoom_2D(pressure.ps_2D, height, x_start, x_stop, y_start, y_stop)
                graphics.plot_contour_mesh(p_zoom/p_scale, xs_zoom/x_scale, ys_zoom/y_scale, p_title, p_labels, vmin=p_min, vmax=p_max, log_cmap=log_cmap_on,n_contours=100)
            else:
                graphics.plot_contour_mesh(pressure.ps_2D/p_scale, height.xs/x_scale, height.ys/y_scale, p_title, p_labels, vmin=p_min/p_scale, vmax=p_max/p_scale, log_cmap=log_cmap_on)
        
        else:
            paramstr = "$Q=%.2f$, $U=%.1f$, $\Delta p=%.2f$"%(flux, self.BC.U, dp_dim)
            p_title = solver_title +'\n' + paramstr
            p_labels = ["$p$", "$x$","$y$"]
           
            if zoom: 
                xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
                p_zoom = graphics.grid_zoom_2D(pressure.ps_2D, height, x_start, x_stop, y_start, y_stop)
                graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, p_title, p_labels, vmin=p_min, vmax=p_max, log_cmap=log_cmap_on,n_contours=50)
            else:
                graphics.plot_contour_mesh(pressure.ps_2D, height.xs, height.ys, p_title, p_labels, vmin=p_min, vmax=p_max, log_cmap=log_cmap_on)
        
    
    def v_plot(self, BC, height, velocity, pressure, solver_title, scaled=False, zoom=False):
                
        x_scale = height.xs[-1]-height.xs[0]
        y_scale = max(height.hs)
        u_scale = velocity.Q/y_scale
        v_scale = velocity.Q/x_scale
        p_scale = velocity.Q*x_scale/y_scale #*visc

        
        ps_2D = np.nan_to_num(pressure.ps_2D)
    
        if scaled:
            dp_nondim = (sum(ps_2D[:,0])/height.hs[0] - sum(ps_2D[:,-1])/height.hs[-1])*height.dy/p_scale/y_scale
            paramstr = "$Q=%.2f$, $U_b=%.2f$, $\Delta P=%.2f$"%(1, self.BC.U/u_scale, dp_nondim)
            v_title = solver_title + '\n' + paramstr
            v_ax_labels =  ['$|(U, V)|_2$','$X$', '$Y$']
            uv_mag = np.sqrt((velocity.u/u_scale)**2 + (velocity.v/v_scale)**2)
                
            if zoom:
                xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
                u_2D_zoom = graphics.grid_zoom_2D(velocity.u, height, x_start, x_stop, y_start, y_stop)
                v_2D_zoom = graphics.grid_zoom_2D(velocity.v, height, x_start, x_stop, y_start, y_stop)
                uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, height, x_start, x_stop, y_start, y_stop)
         
                graphics.plot_stream_heat(u_2D_zoom/u_scale, v_2D_zoom/v_scale, xs_zoom/x_scale, ys_zoom/y_scale, uv_mag_zoom, v_title, v_ax_labels, vmin=0, vmax=vel_max, log_cmap=log_cmap_on) 
            else:
                graphics.plot_stream_heat(velocity.u/u_scale, velocity.v/y_scale, height.xs/x_scale, height.ys/y_scale, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vel_max/velocity.Q, log_cmap=log_cmap_on)

        else:
            dp_dim = (sum(ps_2D[:,0])/height.hs[0] - sum(ps_2D[:,-1])/height.hs[-1])*height.dy
            paramstr = "$Q=%.2f$, $U=%.2f$, $\Delta p=%.2f$"%(velocity.Q, self.BC.U, dp_dim)
            v_title = solver_title + '\n' + paramstr
            v_ax_labels =  ['$|(  u,  v)|_2$','$x$', '$y$']  
            uv_mag = np.sqrt((velocity.u)**2 + (velocity.v)**2)
            
            if zoom:
                xs_zoom, ys_zoom = graphics.grid_zoom_1D(height.xs, height.ys, height, x_start, x_stop, y_start, y_stop)
                u_2D_zoom = graphics.grid_zoom_2D(velocity.u, height, x_start, x_stop, y_start, y_stop)
                v_2D_zoom = graphics.grid_zoom_2D(velocity.v, height, x_start, x_stop, y_start, y_stop)
                uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, height, x_start, x_stop, y_start, y_stop)
         
                graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, v_title, v_ax_labels, vmin=0, vmax=vel_max, log_cmap=log_cmap_on) 
            else:
            
               graphics.plot_stream_heat(velocity.u, velocity.v, height.xs, height.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vel_max, log_cmap=log_cmap_on)

      