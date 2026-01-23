# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain
import boundary as bc
import reyn_velocity_ELT as elt_v

class Velocity:
    def __init__(self, Q, u, v):

        self.u = u
        self.v = v
        self.Q = Q
        
            
class Reyn_Velocity(Velocity):
    def __init__(self, height, BC, ps=None) :
        if isinstance(BC, bc.Fixed):
            h0=height.hs[0]
            px0 = domain.center_first(height.dx, ps[0:3])
            Q = (BC.U*h0)/2 - (px0*(h0**3))/12 #/visc
            
        elif isinstance(BC, bc.Mixed):
            Q = BC.Q
            
        u, v = self.make_velocity(height, BC.U, Q)
        super().__init__(Q, u, v)
        
    # 2D velocity field from 1D pressure 
    def make_velocity(self, height, U, Q):

        u = np.zeros((height.Ny, height.Nx))
        v = np.zeros((height.Ny, height.Nx))

        for i in range(height.Nx):

            h = height.hs[i]
            hx = height.hxs[i]

            for j in range(height.Ny):
                y = height.ys[j]
                
                if y <= height.hs[i]:
                    u[j,i] = (h-y)*(U*(h-3*y)/h**2 + 6*Q*y/h**3)
                    v[j,i] = -2*hx * y**2 * (h-y) *(U/h**3 - 3*Q/h**4)
                    
                else:
                    u[j,i] = 0
                    v[j,i] = 0
                    
        return u, v


class VA_ELT_Velocity(Velocity):
    
    def __init__(self, height, BC, adj_pressure):
        u, v = elt_v.make_ELT_velocity(height, BC, adj_pressure)
        if isinstance(BC, bc.Mixed):
            Q = BC.Q
        else:
            h0=height.hs[0]
            px0 = domain.right_first(height.dx, adj_pressure.ps_2D[0,0:3])
            Q = (BC.U*h0)/2 - (px0*(h0**3))/12 #/visc
        super().__init__(Q, u, v)
        
