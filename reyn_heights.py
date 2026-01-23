 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import numpy as np
from domain import Height

#------------------------------------------------------------------------------

class PWL_Height(Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions # = len(h_peaks)-1
        
        hs, self.slopes, self.widths, i_peaks = self.make_hs(x0, xf, N, x_peaks, h_peaks)
        
        y0 = 0
        yf = max(hs)  
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks)
        
    def make_hs(self, x0, xf, N, x_peaks, h_peaks):
        slopes = np.zeros(self.N_regions)
        widths = np.zeros(self.N_regions)
        i_peaks = np.zeros(self.N_regions+1, dtype=int)
        for r in range(self.N_regions):            
            slopes[r] = (h_peaks[r+1,0] - h_peaks[r,1])/(x_peaks[r+1] - x_peaks[r])
            
        Nx = int((xf-x0)*N + 1)
        hs = np.zeros(Nx)
        dx = 1/N
        
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
            
            i_peaks[r+1] = i    
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1] + slopes[r] * (xi - x_peaks[r])

        return  hs, slopes, widths, i_peaks
    
#------------------------------------------------------------------------------

class PWC_Height(Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks):

        #solver requires minimum 3 regions
        while N_regions <3 :
            x_peak_mid = x_peaks[-1] - x_peaks[-2]/2
            x_peak_end = x_peaks[-1]
            x_peaks[-1] = x_peak_mid
            x_peaks = np.append(x_peaks, x_peak_end)
            h_peaks = np.append(h_peaks, h_peaks[-1])
            h_peaks = np.reshape(h_peaks, (N_regions+2,2))
            N_regions +=1 
            
        self.N_regions = N_regions
        self.x_peaks=x_peaks
        self.h_peaks=h_peaks
        self.h_steps = h_peaks[:,1][:-1] 
        self.slopes =np.zeros(self.N_regions)
        
        hs, self.widths, i_peaks = self.make_hs(x0, xf, N, x_peaks, h_peaks)
        y0 = 0
        yf = max(hs)  
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks)
        
    def make_hs(self, x0, xf, N, x_peaks, h_peaks):
        widths = np.zeros(self.N_regions)
        i_peaks = np.zeros(self.N_regions+1, dtype=int)
       
        Nx = int((xf-x0)*N + 1)
        hs = np.zeros(Nx)
        dx = 1/N
        
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
            i_peaks[r+1] = i    
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1]
        return  hs, widths, i_peaks
    
#------------------------------------------------------------------------------

class LogisticHeight(Height):
    def __init__(self, x0, xf, N, H, h, center, delta):
       
        Nx = (xf-x0)*N + 1
        dx = 1/N
        self.h = h
        self.H = H
        self.delta = delta #slope = delta*(H-h)/4
        self.center = center
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])  
        y0 = 0
        yf = np.max(hs)
        i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks)

    def h_fun(self, x):
        
        return (self.h + (self.H-self.h) / ( 1 + np.exp(self.delta*(x-self.center))))  
