#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:30:05 2025

@author: sarahdennis
"""
import numpy as np
import domain as dm
import graphics
import boundary as bc

from reyn_pressure import Pressure
from reyn_velocity import Velocity
from reyn_heights import PWL_Height
# H=1
# slope_k = 1/10
 
class PerturbedReynSol:
    def __init__(self, height, BC, order, reyn_pressure, reyn_velocity):
        self.order = order
        if order < 0 or order > 4:
            return Exception(f"order {order} not in range [0,4]")
        
        self.x_scale = (height.xf - height.x0)/2 
        self.y_scale = height.yf - height.y0
        
        self.Q_scale = reyn_velocity.Q

        self.P_scale = self.Q_scale * self.x_scale * (self.y_scale**-3)
        
        self.U_scale = self.Q_scale/self.y_scale
        self.V_scale = self.Q_scale/self.x_scale                      

        self.u0s = reyn_velocity.u /self.U_scale
        self.v0s = reyn_velocity.v /self.V_scale
        
        reyn_pressure.make_2D_ps(height)
        self.p0s = reyn_pressure.ps_2D /self.P_scale
        
        delta = self.y_scale/self.x_scale

        self.is_pwl_height = isinstance(height, PWL_Height)

        if order > 1:

            # make self.u2s, self.v2s, self.p2s, etc... 
            self.perturb_second(height)
            pert2_ps_2D = (self.p0s + (delta**2) * self.p2s) *self.P_scale
            pert2_us_2D = (self.u0s + (delta**2) * self.u2s) *self.U_scale
            pert2_vs_2D = (self.v0s + (delta**2) * self.v2s) *self.V_scale
            self.pert2_pressure = Pressure(height, BC, ps_1D = reyn_pressure.ps_1D, ps_2D=pert2_ps_2D)
            self.pert2_velocity = Velocity(reyn_velocity.Q, pert2_us_2D, pert2_vs_2D)

        if order > 2: 

            # make self.u4s, self.v4s, self.p4s, etc... 
            self.perturb_fourth(height)
            pert4_ps_2D = pert2_ps_2D + (delta**4) * self.p4s *self.P_scale
            pert4_us_2D = pert2_us_2D + (delta**4) * self.u4s *self.U_scale
            pert4_vs_2D = pert2_vs_2D + (delta**4) * self.v4s *self.V_scale
            self.pert4_pressure = Pressure(height, BC, ps_1D = reyn_pressure.ps_1D, ps_2D=pert4_ps_2D)
            self.pert4_velocity = Velocity(reyn_velocity.Q, pert4_us_2D, pert4_vs_2D)
    
        
    def perturb_second(self, height):

        ys = height.ys/self.y_scale
        hs = height.hs/self.y_scale
        
        h_xs = np.zeros(height.Nx)     # d/dx[h] @ xi
        h2_2xs = np.zeros(height.Nx)   # d^2/dx^2 [h^-2] @ xi
        h3_2xs = np.zeros(height.Nx)   # d^2/dx^2 [h^-3] @ xi
        h2_3xs = np.zeros(height.Nx)   # d^3/dx^3 [h^-2] @ xi
        h3_3xs = np.zeros(height.Nx)   # d^3/dx^3 [h^-3] @ xi
                
        c3_2xs = np.zeros(height.Nx)                # d^2/dx^2 [c3(x)]
        c3_xs = np.zeros(height.Nx)                 # d/dx [c3(x)]
        c3s = np.zeros(height.Nx)                   # intS d/dx[c3(x)] dx
        u0_xs = np.zeros((height.Ny, height.Nx))    # d/dx [reyn_velx] 
        
        u2s = np.zeros((height.Ny, height.Nx)) #
        v2s = np.zeros((height.Ny, height.Nx)) #

        dx = height.dx/(self.x_scale)
        
        # find h_x, h2_2x, h2_3x, h3_2x, h3_3x
        for i in range(height.Nx):
            h = hs[i]            
            
            if i < 2: #inlet
                h2_3x = dm.right_third(dx, [h**-2 for h in hs[i : i+5]])
                h3_3x = dm.right_third(dx,  [h**-3 for h in hs[i : i+5]])
                
                if i == 0: #
                    h_x = dm.right_first(dx, hs[i : i+3])
                    h2_2x = dm.right_second(dx, [h**-2 for h in hs[i : i+4]])
                    h3_2x = dm.right_second(dx, [h**-3 for h in hs[i : i+4]])
                    
                else: # i == 1
                    h_x = dm.center_first(dx, hs[i-1 : i+2])
                    h2_2x = dm.center_second(dx, [h**-2 for h in hs[i-1 : i+2]])
                    h3_2x = dm.center_second(dx, [h**-3 for h in hs[i-1 : i+2]])

            elif i > height.Nx-3: #outlet
                h2_3x = dm.left_third(dx, [h**-2 for h in hs[i-4 : i+1]])
                h3_3x = dm.left_third(dx, [h**-3 for h in hs[i-4 : i+1]])
                
                if i == height.Nx-1:
                    h_x = dm.left_first(dx, hs[i-2 : i+1])
                    h2_2x = dm.left_second(dx, [h**-2 for h in hs[i-3 : i+1]])
                    h3_2x = dm.left_second(dx, [h**-3 for h in hs[i-3 : i+1]])
                    
                else: # i == Nx-2
                    h_x = dm.center_first(dx, hs[i-1 : i+2])
                    h2_2x = dm.center_second(dx, [h**-2 for h in hs[i-1 : i+2]])
                    h3_2x = dm.center_second(dx, [h**-3 for h in hs[i-1 : i+2]])
                
            else: #interior  
                h_x = dm.center_first(dx, hs[i-1 : i+2])
                h2_2x = dm.center_second(dx, [h**-2 for h in hs[i-1 : i+2]])
                h3_2x = dm.center_second(dx, [h**-3 for h in hs[i-1 : i+2]])
                h2_3x = dm.center_third(dx, [h**-2 for h in hs[i-2 : i+3]])
                h3_3x = dm.center_third(dx, [h**-3 for h in hs[i-2 : i+3]])
            
            h_xs[i]   = h_x
            h2_2xs[i] = h2_2x
            h3_2xs[i] = h3_2x
            h2_3xs[i] = h2_3x
            h3_3xs[i] = h3_3x
            
            # make c3_x, c3_2x
            c3_xs[i]  = 6 * h2_2x * h - 18/5 * h3_2x * (h**2)
            c3_2xs[i] = 6 * (h2_2x * h_x +  h2_3x * h) - 18/5* (2*h*h_x * h3_2x + (h**2) * h3_3x)
        

        # correct finite differences at discontinuities
        for i in height.i_peaks[1:-1]:
            h_xs[i-1 : i+2] = dm.avg_x(h_xs[i-2 : i+3])
            h2_2xs[i-1 : i+2] = dm.avg_2x(h2_2xs[i-2 : i+3])
            h3_2xs[i-1 : i+2] = dm.avg_2x(h3_2xs[i-2 : i+3])
            h2_3xs[i-2 : i+3] = dm.avg_3x(h2_3xs[i-3 : i+4])
            h3_3xs[i-2 : i+3] = dm.avg_3x(h3_3xs[i-3 : i+4])
        
            c3_xs[i-1 : i+2] = dm.avg_x(c3_xs[i-2 : i+3])
            c3_2xs[i-2 : i+3] = dm.avg_3x(c3_2xs[i-3 : i+4])
            
        # find u2, v2   
        for i in range(height.Nx):
            h = hs[i]    
            h_x = h_xs[i]
            h2_2x = h2_2xs[i]
            h2_3x = h2_3xs[i]
            h3_2x = h3_2xs[i]
            h3_3x = h3_3xs[i]
            c3_x = c3_xs[i]
            c3_2x = c3_2xs[i]

            
            # make u2_v2
            for j in range (height.Ny):
                y = ys[j]
    
                if y <= h:
                    u0_xs[j,i] = 6*h_x * (-2*y*(h**-3) + 3*(y**2)* (h**-4))
    
                    u2_A = -2 * h2_2x * (y**3 - (h**2) * y) 
                    u2_B = h3_2x * (y**4 - (h**3) * y)
                    u2_C = (1/2) * c3_x * (y**2 - h * y)
                    u2s[j,i] = (u2_A + u2_B +  u2_C)
                    
                    v2_A = h2_3x * ((1/2) * (y**4) - (h**2) * (y**2))
                    v2_B = -2 * h2_2x * h_x * h * (y**2)
                    v2_C = -h3_3x * ((1/5) * (y**5) - (1/2) * (h**3) * (y**2))
                    v2_D = (3/2) * h3_2x * h_x * (h**2) * (y**2)
                    v2_E = -(1/2)*c3_2x* ((1/3) * (y**3) - (1/2)* h *(y**2))
                    v2_F = (1/4) * c3_x * h_x * (y**2)
                    v2s[j,i] = (v2_A + v2_B + v2_C + v2_D + v2_E + v2_F)

                    
        # p2s[j,0] = 0 (no correction at inlet)
        
        # integrate x0 -> xL
        c3s[0] = np.mean(u0_xs[:,0]) #u0_xs[0]
        c3s[1] = c3s[0] + c3_xs[1]*dx
        for i in range(2, height.Nx):
            c3s[i] =(4*c3s[i-1] -c3s[i-2] + 2*dx*c3_xs[i])/3
    
        c3s -= c3s[-1]
        p2s = -u0_xs + c3s 
        
        self.p2s = p2s
        self.u2s = u2s
        self.v2s = v2s
        
        self.c3s = c3s     # intS [x0, xL] d/dx [c3] dx 
        self.c3_xs = c3_xs # d/dx [c3]
        self.c3_2xs = c3_2xs # d^2/dx^2 

        self.h_xs = h_xs       # d/dx   [h]    
        self.h2_2xs = h2_2xs   # d^2/dx^2 [h^-2] 
        self.h3_2xs = h3_2xs   # d^2/dx^2 [h^-3] 
        self.h2_3xs = h2_3xs   # d^3/dx^3 [h^-2] 
        self.h3_3xs = h3_3xs   # d^3/dx^3 [h^-3] 

    def perturb_fourth(self, height):

        ys = height.ys/self.y_scale
        hs = height.hs/self.y_scale

        # find 8 intermediate derivatives
        h2_5xs = np.zeros(height.Nx)
        h3_5xs = np.zeros(height.Nx)
        h2_4xs = np.zeros(height.Nx)
        h3_4xs = np.zeros(height.Nx)
        h_3xs = np.zeros(height.Nx)
        h_2xs = np.zeros(height.Nx)
        c3_4xs = np.zeros(height.Nx)
        c3_3xs = np.zeros(height.Nx)
        
        # integration constant 
        f1_2xs = np.zeros(height.Nx)
        f1_3xs = np.zeros(height.Nx)
        f2_2xs = np.zeros(height.Nx)
        f2_3xs = np.zeros(height.Nx)
        f3_2xs = np.zeros(height.Nx)
        f3_3xs = np.zeros(height.Nx)
        c5_2xs = np.zeros(height.Nx)        # d/dx [c5(x)]
        c5_xs = np.zeros(height.Nx)         # d/dx [c5(x)]
        c5s = np.zeros(height.Nx)           # intS d/dx[c5(x)] dx
        
        # p4s = -u2_xs + dxx int_0^y [v0] dy + c5s
        u2_xs = np.zeros((height.Ny, height.Nx))        # d/dx [u2] 
        v0_Sys = np.zeros((height.Ny, height.Nx))       # int_0^y [v0] dy
        v0_Sy_2xs = np.zeros((height.Ny, height.Nx))    # d^2/dx^2 int_0^y [v0] dy
        u4s = np.zeros((height.Ny, height.Nx))          #
        v4s = np.zeros((height.Ny, height.Nx))          #
        
        dy = height.dy/(self.y_scale)
        dx = height.dx/(self.x_scale)

         # find h_2x, h_3x, h2_4x, h2_5x, h3_4x, h3_5x
        for i in range(height.Nx):
   
            h = hs[i]
            h_x = self.h_xs[i]
            h2_2x = self.h2_2xs[i]
            h3_2x = self.h3_2xs[i]
            h2_3x = self.h2_3xs[i]
            h3_3x = self.h3_3xs[i]
            c3_x = self.c3_xs[i] 
            c3_2x = self.c3_2xs[i]
            
            if i < 3: #inlet
                h2_5x = dm.right_fifth(dx,  [h**-2 for h in hs[i : i+7]])
                h3_5x = dm.right_fifth(dx,  [h**-3 for h in hs[i : i+7]])
                
                if i < 2: 
                    h2_4x = dm.right_fourth(dx, [h**-2 for h in hs[i : i+6]])
                    h3_4x = dm.right_fourth(dx, [h**-3 for h in hs[i : i+6]])
                    h_3x = dm.right_third(dx,  hs[i : i+5])          

                    if i == 0:
                        h_2x = dm.right_second(dx, hs[i : i+4])
                    else: #i == 1
                        h_2x = dm.center_second(dx, hs[i-1 : i+2])
                        
                else: #i == 2
                    h2_4x = dm.center_fourth(dx, [h**-2 for h in hs[i-2 : i+3]])
                    h3_4x = dm.center_fourth(dx, [h**-3 for h in hs[i-2 : i+3]])
                    h_3x = dm.center_third(dx, hs[i-2 : i+3])
                    h_2x = dm.center_second(dx, hs[i-1 : i+2])
                
            elif i > height.Nx-4: #outlet
                h2_5x = dm.left_fifth(dx,  [h**-2 for h in hs[i-6 : i+1]])
                h3_5x = dm.left_fifth(dx,  [h**-3 for h in hs[i-6 : i+1]])
                
                if i > height.Nx-3:
                    h2_4x = dm.left_fourth(dx, [h**-2 for h in hs[i-5 : i+1]])
                    h3_4x = dm.left_fourth(dx, [h**-3 for h in hs[i-5 : i+1]])
                    h_3x = dm.left_third(dx,  hs[i-4 : i+1])
                    
                    if i == height.Nx-1:
                        h_2x = dm.left_second(dx, hs[i-3 : i+1])
                    else: # i == Nx-2
                        h_2x = dm.center_second(dx, hs[i-1 : i+2])
                        
                else: #i == Nx-3
                    h_2x = dm.center_second(dx, hs[i-1 : i+2])
                    h_3x = dm.center_third(dx, hs[i-2 : i+3])
                    h2_4x = dm.center_fourth(dx, [h**-2 for h in hs[i-2 : i+3]])
                    h3_4x = dm.center_fourth(dx, [h**-3 for h in hs[i-2 : i+3]])
               
            else:  #interior 
                h_2x = dm.center_second(dx, hs[i-1 : i+2])
                h_3x = dm.center_third(dx, hs[i-2 : i+3])
                h2_4x = dm.center_fourth(dx, [h**-2 for h in hs[i-2 : i+3]])
                h3_4x = dm.center_fourth(dx, [h**-3 for h in hs[i-2 : i+3]])
                h2_5x = dm.center_fifth(dx, [h**-2 for h in hs[i-3 : i+4]])
                h3_5x = dm.center_fifth(dx, [h**-3 for h in hs[i-3 : i+4]])
               
            h_2xs[i] = h_2x
            h_3xs[i] = h_3x
            h2_4xs[i] = h2_4x
            h3_4xs[i] = h3_4x
            h2_5xs[i] = h2_5x
            h3_5xs[i] = h3_5x
            
            # make c5_x, c5_2x, c3_3x, c3_4x
        
            c3_3x_A = 6*(2*h2_3x*h_x + h2_2x*h_2x + h2_4x*h)
            c3_3x_B = (-18/5)*(2*(h_x**2)*h3_2x + 2*h*h_2x*h3_2x + 4*h*h_x*h3_3x + (h**2)*h3_4x)
            c3_3x = c3_3x_A + c3_3x_B
            c3_3xs[i] = c3_3x 
            
            c3_4x_A = 6*(3*h2_4x*h_x + 3*h2_3x*h_2x + h2_2x*h_3x + h2_5x*h)
            c3_4x_B = (-18/5)*(6*(h_x*h_2x*h3_2x+ (h_x**2)*h3_3x + h*h_2x*h3_3x + h*h_x*h3_4x) + 2*h*h_3x*h3_2x + (h**2)*h3_5x)
            c3_4x = c3_4x_A + c3_4x_B
            c3_4xs[i] = c3_4x
            
            f1_2x = (6*h*(h_x**2) + 3*(h**2)*h_2x)*h3_2x + 6*(h**2)*h_x*h3_3x + (h**3)*h3_4x
            f1_2xs[i] = f1_2x 
            
            f1_3x_A = (6*(h_x**3) + 18*h*h_x*h_2x + 3*(h**2)*h_3x) *h3_2x
            f1_3x_B = (9*(h**2)*h_2x + 18*h*(h_x**2)) *h3_3x + 9*(h**2)*h_x*h3_4x + (h**3)*h3_5x
            f1_3x = f1_3x_A + f1_3x_B
            f1_3xs[i] = f1_3x 
            
            f2_2x = 2*(h_x**2 + h*h_2x)*h2_2x + 4*h*h_x*h2_3x + (h**2)*h2_4x
            f2_2xs[i] = f2_2x 

            f2_3x_A = (6*h_x*h_2x + 2*h*h_3x)*h2_2x
            f2_3x_B = 6*(h*h_2x + h_x**2)*h2_3x + 6*h*h_x*h2_4x + (h**2)*h2_5x
            f2_3x = f2_3x_A + f2_3x_B 
            f2_3xs[i] = f2_3x 
    
            f3_2x = h_2x*c3_x + 2*h_x*c3_2x + h*c3_3x
            f3_2xs[i] = f3_2x 
            
            f3_3x = h_3x*c3_x + 3*h_2x*c3_2x + 3*h_x*c3_3x + h*c3_4x
            f3_3xs[i] = f3_3x 
    
            c5_x_A = (3/14)*h3_4x*(h**4) - (3/5)*h2_4x*(h**3)+ (3/10)*c3_3x*(h**2)
            c5_x_B = -(f1_2x - 2*f2_2x)*h  -(1/2)*f3_2x*h
            c5_x =  c5_x_A + c5_x_B
            c5_xs[i] = c5_x 
            
            c5_2x_A = (3/14)*(h3_5x*(h**4) + 4*h3_4x*h_x*(h**3)) 
            c5_2x_B = -(3/5)*(h2_5x*(h**3) + 3*h2_4x*h_x*(h**2))
            c5_2x_C = -((f1_3x - 2*f2_3x)*h + (f1_2x - 2*f2_2x)*h_x) 
            c5_2x_D = (3/10)*(c3_4x*(h**2) + 2*c3_3x*h_x*h) 
            c5_2x_E = -(1/2)*(f3_3x*h + f3_2x*h_x) 
            c5_2x = c5_2x_A + c5_2x_B + c5_2x_C + c5_2x_D + c5_2x_E 
            c5_2xs[i] = c5_2x 
            

        # correct finite differences at discontinuities
        for i in height.i_peaks[1:-1]:
            h_2xs[i-1 : i+2] = dm.avg_2x(h_2xs[i-2 : i+3])
            h_3xs[i-2 : i+3] = dm.avg_3x(h_3xs[i-3 : i+4])
            h2_4xs[i-3 : i+4] = dm.avg_4x(h2_4xs[i-4 : i+5])
            h3_4xs[i-3 : i+4] = dm.avg_4x(h3_4xs[i-4 : i+5])
            h2_5xs[i-4 : i+5] = dm.avg_5x(h2_5xs[i-5 : i+6])
            h3_5xs[i-4 : i+5] = dm.avg_5x(h3_5xs[i-5 : i+6])
             
            c3_3xs[i-3 : i+4] = dm.avg_4x(c3_3xs[i-4 : i+5])
            c3_4xs[i-4 : i+5] = dm.avg_5x(c3_4xs[i-5 : i+6])
            c5_xs[i-3 : i+4] = dm.avg_4x(c5_xs[i-4 : i+5])
            c5_2xs[i-4 : i+5] = dm.avg_5x(c5_2xs[i-5 : i+6])
            
            f1_2xs[i-3 : i+4] = dm.avg_4x(f1_2xs[i-4 : i+5])
            f1_3xs[i-4 : i+5] = dm.avg_5x(f1_3xs[i-5 : i+6])
            f2_2xs[i-3 : i+4] = dm.avg_4x(f2_2xs[i-5 : i+5])
            f2_3xs[i-4 : i+5] = dm.avg_5x(f2_3xs[i-5 : i+6])
            f3_2xs[i-3 : i+4] = dm.avg_4x(f3_2xs[i-4 : i+5])
            f3_3xs[i-4 : i+5] = dm.avg_5x(f3_3xs[i-5 : i+6])                        

        
        # find u4, v4
        for i in range(height.Nx):
            v0_Sy_i = 0
            v0_i_js = self.v0s[:,i]
            
            h = hs[i]
            
            h_x = self.h_xs[i]
            h2_2x = self.h2_2xs[i]
            h3_2x = self.h3_2xs[i]
            h2_3x = self.h2_3xs[i]
            h3_3x = self.h3_3xs[i]
            c3_x = self.c3_xs[i] 
            c3_2x = self.c3_2xs[i]
            
            h2_5x = h2_5xs[i] 
            h3_5x = h3_5xs[i]
            h2_4x = h2_4xs[i]
            h3_4x = h3_4xs[i]
            h_3x = h_3xs[i]
            h_2x = h_2xs[i]
            
            c3_3x = c3_3xs[i]
            c3_4x = c3_4xs[i]
            c5_x = c5_xs[i] 
            c5_2x = c5_2xs[i]
            
            f1_2x = f1_2xs[i]
            f1_3x = f1_3xs[i]
            f2_2x = f2_2xs[i]
            f2_3x = f2_3xs[i]
            f3_2x = f3_2xs[i]
            f3_3x = f3_3xs[i]
            
            # make u4, v4
            for j in range (height.Ny):
                y = ys[j]
                v0_Sy_i += v0_i_js[j] *  dy
                    
                v0_Sys[j,i] = v0_Sy_i

                
                if y <= h:

                    u2x_A = -2*h2_3x*((y**3)-(h**2)*y) + 4*h2_2x*h_x*h*y
                    u2x_B = h3_3x*((y**4)-(h**3)*y) -3*h3_2x*h_x*(h**2)*y
                    u2x_C = (1/2)*(c3_2x *((y**2)-h*y) - h_x*y) 
                    u2_xs[j,i] = u2x_A + u2x_B + u2x_C
                    
                    u4_A = (-1/20)*h3_4x*((y**6) - (h**5)*y) + (3/20)*h2_4x*((y**5) - (h**4)*y)
                    u4_B = ((1/3)*f1_2x - (2/3)*f2_2x+ (1/6)*f3_2x)*((y**3) - (h**2)*y)
                    u4_C = (-1/12)*c3_3x*((y**4) - (h**3)*y) + (1/2)*c5_x*((y**2)-h*y)
                    u4s[j,i] = (u4_A + u4_B + u4_C)
 
                    v4_A = (-1/20)*h3_5x*((1/7)*(y**7) - (1/2)*(h**5)*(y**2)) + (1/8)*h3_4x*h_x*(h**4)*(y**2)
                    v4_B = (3/20)*h2_5x*((1/6)*(y**6) - (1/2)*(h**4)*(y**2)) - (3/10)*h2_4x*h_x*(h**3)*(y**2)
                    v4_C = (1/3)*((f1_3x - 2*f2_3x)*((1/4)*(y**4) - (1/2)*(h**2)*(y**2)) - (f1_2x - 2*f2_2x)*h_x*h*(y**2))
                    v4_D = (1/6)*(f3_3x*((1/4)*(y**4) - (1/2)*(h**2)*(y**2)) - f3_2x*h_x*h*(y**2))
                    v4_E = (-1/12)*(c3_4x)*((1/5)*(y**5) - (1/2)*(h**3)*(y**2)) + (1/8)*c3_3x*h_x*(h**2)*(y**2)
                    v4_F = (1/2)*(c5_2x)*((1/3)*(y**3) - (1/2)*h*(y**2)) - (1/4)*c5_x*h_x*(y**2)
                    v4s[j,i] = -(v4_A + v4_B + v4_C + v4_D + v4_E + v4_F)
        
        # make p4s
         
        for j in range(height.Ny): 
            v0_Sy_2xs[j] = dm.center_second_diff(v0_Sys[j], height.Nx, dx)
            
            for i in height.i_peaks[1:-1]:
                v0_Sy_2xs[j, i-2 : i+3] = dm.avg_3x(v0_Sy_2xs[j, i-3 : i+4])
                
        p4s = np.zeros((height.Ny, height.Nx))

        c5s[0] = 0
        c5s[1] = c5s[0] + c5_xs[1]*dx
        
        for i in range(2, height.Nx):
            c5s[i] =(4*c5s[i-1] -c5s[i-2] + 2*dx*c5_xs[i])/3
            
        c5s -=c5s[-1]


        p4s = -u2_xs + v0_Sy_2xs + c5s 

        self.p4s = p4s
        self.u4s = u4s
        self.v4s = v4s

