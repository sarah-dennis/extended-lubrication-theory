# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
import math
from domain import Space

#------------------------------------------------------------------------------
class PWLinear(Space):
    # peak_xs = [x0, ..., xi,..., xf] : x0 < xi < xf
    # peak_ys = [(yf,h_in),...,(hi_left, hi_right),...,(h_out,yf)]
    
    def __init__(self, x0, xf, y0, yf, N, BC, Re, namestr, x_peaks, y_peaks):
        super().__init__(x0, xf, y0, yf, N, BC, Re, namestr)
        
        # peaks must fall on the grid
        self.x_peaks = x_peaks
        self.y_peaks = y_peaks
        self.N_regions = len(self.x_peaks)-1
        self.make_space()

                              
        self.H_in = yf - y_peaks[0][1] 
        self.H_out = yf - y_peaks[-1][0]
        
        self.spacestr = "$Q=%.2f$, $U=%.1f$"%(BC.Q,BC.U)  # for plot title
        
        if self.H_in == 0: # closed cavity --> Q=0, dp=0
            self.dp_in = 0
        else: # gap entry --> dp ~ Q
            self.dp_in = (BC.Q - 0.5*BC.U*self.H_in) * (-12 / self.H_in**3)
    
    # 0 : boundary, -1: exterior, 1: interior
    def make_space(self):

        slopes = np.zeros(self.N_regions)
        
        for k in range(self.N_regions):
            dh = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            slopes[k] = dh/dx
        hs = np.zeros((self.Nx,2))
        grid = np.zeros((self.Ny, self.Nx))
    
        reg = 0        
        i_ref = 0
        for i in range(self.Nx):
            if math.isclose(self.xs[i], self.x_peaks[reg]):
                h_left = self.y_peaks[reg][0]
                h_right = self.y_peaks[reg][1] 
                i_ref = i
                reg +=1
                hs[i] = [h_left,h_right]
            else:
                h = slopes[reg-1]*(i - i_ref)*self.dx + self.y_peaks[reg-1][1]
                hs[i] = [h,h]
                
            for j in range(self.Ny):
                y = self.ys[j]
                if j == self.Ny-1: #upper boundary
                    grid[j,i] = 0
                    
                elif i == 0: #inlet boundary
                    if y >= self.y_peaks[0][1]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i] = -1
                    
                elif i == self.Nx-1:# outlet boundary
                    if y >= self.y_peaks[-1][0]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i] = -1
                
                elif i == i_ref: # pwl region change              
                    if math.isclose(y, h_left) or math.isclose(y, h_right): # true boundary point at region change
                        grid[j,i] = 0
                        
                    elif h_left < h_right:
                        if h_left < y and y < h_right: # x=h(y) vertical boundary
                            grid[j,i] = 0
                        elif y > h_right: # above vert boundary (interior)
                            grid[j,i] = 1
                        else:             # below vert boundary(exterior)
                            grid[j,i] = -1
                
                    else:
                        if h_left > y and y > h_right: # x=h(y) vertical boundary
                            grid[j,i] = 0
                        elif y > h_left: # above vert boundary (interior)
                            grid[j,i] = 1
                        else:              # below vert boundary (exterior)
                            grid[j,i] = -1

                else:
                    if math.isclose(y,h): # true boundary point not at region change (from dx | slope)
                        grid[j,i] = 0
                    elif y > h:            # above boundary (interior)

                        grid[j,i] = 1
                    else:                   # below boundary (exterior)
                        grid[j,i] = -1
        
        self.space = grid
        # self.slopes = slopes
        self.hs = hs

#------------------------------------------------------------------------------
# Boundary conditions on stream and velocity
#------------------------------------------------------------------------------
   
    def streamInlet(self, j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:
      
            u_term = self.BC.U*(0.5*(y**2 - self.yf**2) + (self.H_in-self.yf)*(y-self.yf))/self.H_in
        
            dp_term = -0.5*self.dp_in*((-1/3)*(y**3 -self.yf**3) + 0.5*(2*self.yf -self.H_in)*(y**2-self.yf**2) + self.yf*(self.H_in-self.yf)*(y-self.yf))

            psi = u_term + dp_term + self.BC.Q
            return psi
        else:
            return 0
    
    def velInlet(self, j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:

            u = -(self.BC.U/self.H_in - 0.5*self.dp_in * (self.yf-y)) * (self.yf-self.H_in - y) 
            
            return  u
        else: 
            return 0
        
#------------------------------------------------------------------------------
# Boundary interpolation
#------------------------------------------------------------------------------
    def interp(self, scale, v_opp, v_bdry=0, p=False):
        v_nbr = v_bdry + (v_bdry - v_opp)*scale

        if p:
            print(v_nbr, v_opp)

        return v_nbr

    def scale_S(self, i,j):
        y_N = self.ys[j+1]
        y_nbr = self.ys[j-1]
        y_bdry = self.hs[i][0] # arbitrary 
    
        l1 = np.abs(y_nbr-y_bdry)
        l2 = np.abs(y_N-y_bdry) 
        
        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l1/l2        
    
        return scale
    


    def scale_E(self, i,j):
        x_W = self.xs[i-1]
        h_W = self.hs[i-1][1]
        x_nbr = self.xs[i+1]
        h_nbr = self.hs[i+1][0]
        
        y_bdry = self.ys[j]
        x_bdry = x_W + (x_nbr-x_W) * (y_bdry-h_W)/(h_nbr-h_W)

        l1 = np.abs(x_nbr-x_bdry)
        l2 = np.abs(x_W-x_bdry)
        
        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l1/l2      
        return scale
    
    def scale_W(self, i,j):
        x_E = self.xs[i+1]
        h_E = self.hs[i+1][0]
        x_nbr = self.xs[i-1]
        h_nbr = self.hs[i-1][1]
        
        y_bdry = self.ys[j]
        x_bdry = x_E + (x_nbr-x_E) * (y_bdry-h_E)/(h_nbr-h_E)

        l1 = np.abs(x_nbr-x_bdry)
        l2 = np.abs(x_E-x_bdry)
        
        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l1/l2 

        return scale


    
    def scale_NE(self, i,j): 

        x_SW = self.xs[i-1]
        y_SW = self.ys[j-1]
        h_SW = self.hs[i-1][1]
        
        x_nbr = self.xs[i+1]
        y_nbr = self.ys[j+1]
        h_nbr = self.hs[i+1][0]

        slope = (h_nbr-h_SW)/(x_nbr-x_SW)
        x_bdry = (y_SW-h_SW)/(slope-1) + x_SW
        y_bdry = (x_bdry-x_SW) + y_SW
        
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = np.sqrt((x_SW-x_bdry)**2 + (y_SW-y_bdry)**2)
        
        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l1/l2 

        return scale
        
    def scale_SW(self, i,j): 

        x_NE = self.xs[i+1]
        y_NE = self.ys[j+1]
        h_NE = self.hs[i+1][0]
        
        x_nbr = self.xs[i-1]
        y_nbr = self.ys[j-1]
        h_nbr = self.hs[i-1][1]

        slope = (h_nbr-h_NE)/(x_nbr-x_NE)
        x_bdry = (y_NE-h_NE)/(slope-1) + x_NE
        y_bdry = (x_bdry-x_NE) + y_NE
        
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = np.sqrt((x_NE-x_bdry)**2 + (y_NE-y_bdry)**2)
        
        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l1/l2  
        return scale
    
    def scale_NW(self, i,j): 

        x_SE = self.xs[i+1]
        y_SE = self.ys[j-1]
        h_SE = self.hs[i+1][0]
        
        x_nbr = self.xs[i-1]
        y_nbr = self.ys[j+1]
        h_nbr = self.hs[i-1][1]

        slope = (h_nbr-h_SE)/(x_nbr-x_SE)
        x_bdry = (y_SE-h_SE)/(slope+1) + x_SE
        y_bdry = -(x_bdry-x_SE) + y_SE
        
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = np.sqrt((x_SE-x_bdry)**2 + (y_SE-y_bdry)**2)

        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l1/l2   
        return scale
    
    def scale_SE(self, i,j): 

        x_NW = self.xs[i+1]
        y_NW = self.ys[j-1]
        h_NW = self.hs[i+1][0]
        
        x_nbr = self.xs[i-1]
        y_nbr = self.ys[j+1]
        h_nbr = self.hs[i-1][1]

        slope = (h_nbr-h_NW)/(x_nbr-x_NW)
        x_bdry = (y_NW-h_NW)/(slope+1) + x_NW
        y_bdry = -(x_bdry-x_NW) + y_NW
        
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = np.sqrt((x_NW-x_bdry)**2 + (y_NW-y_bdry)**2)
        
        if np.isclose(l2,0):
            scale = 0
        else:
            scale = l2/l1  

        return scale
#------------------------------------------------------------------------------















