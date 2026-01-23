# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 18:10:34 2025

@author: sarah
"""
import numpy as np
from stokes_heights import PWLinear


class BFS(PWLinear):
    def __init__ (self, args, BC, Re, N):
        h_in, h_out, l_in, l_out = args
        x0 = 0
        xf = l_in+l_out
  
        x_peaks = [0, l_in, l_in+l_out]
        y0 = 0
        yf = max(h_in,h_out) 
        y_peaks=[[yf,yf-h_in],[yf-h_in,h_out],[h_out,yf]]
        namestr= f'BFS_hin{h_in}hout{h_out}lin{l_in}lout{l_out}_U{BC.U}_Q{BC.Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, BC, Re,namestr, x_peaks, y_peaks)
        
class BFS_pwl(PWLinear):
    def __init__ (self, args, BC, Re, N):
        H, L, delta = args
        h=1
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        l = L/2
        x_peaks = [x0, l-delta, l, l+delta, xf]
        
        y_peaks=[[yf,0],[0,0],[yf-h-(H-h)/2,yf-h-(H-h)/2],[yf-h,yf-h],[yf-h,yf]]
        
        namestr= f'BFS_pwl_H{H}L{L}d{delta}_U{BC.U}_Q{BC.Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, BC, Re, namestr, x_peaks, y_peaks)
     
class Logistic(PWLinear):
    def __init__(self, args, BC, Re, N):
        H, h, L, delta = args
        self.delta = delta
        self.h = h
        self.H = H
        self.L = L
        self.x0 = 0
        self.xf = L
        self.center = L//2
        self.h_in = self.h + (self.H-self.h) / ( 1 + np.exp(-self.delta*(self.center-self.x0)))
        
        self.h_out = self.h + (self.H-self.h) / ( 1 + np.exp(-self.delta*(self.center-self.xf)))
        
        self.h_max = max(self.h_in, self.h_out)
        
        x_peaks = [self.x0 + i/N for i in range (int(1 + self.L*N))]
               
        y_peaks_L = [self.h+self.H] + [self.h_fun(x) for x in x_peaks[:-1]]
        y_peaks_R =  [self.h_fun(x) for x in x_peaks[1:]] + [self.h+self.H]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 
        y0 = 0
        yf = self.h_max
        
        namestr= f'logistic_H{H}L{L}d{delta}_U{BC.U}_Q{BC.Q}_Re{Re}'
        super().__init__(self.x0, self.xf, y0, yf, N, BC, Re,namestr, x_peaks, y_peaks)   

    def h_fun(self, x):
        return self.h_max - (self.h + (self.H-self.h) / ( 1 + np.exp(self.delta*(x-self.center))))
     
class Sinusoid(PWLinear):
    def __init__(self, args, BC, Re, N):
        H, delta, k,  l, L = args

        self.H=H
        self.k = k
        self.delta = delta
        self.l = l
        self.L = L
        self.period = np.pi * k / l
        
        self.x0 = -self.L
        self.xf = self.L
        
        y0 = 0
        yf = self.H*(1+self.delta)
        self.yf = yf
        
        
        x_peaks = [self.x0 + i/N for i in range (int(1 + 2*self.L*N))]
        h_peaks = [self.h_fun(x) for x in x_peaks]
        y_peaks_L = [0] + h_peaks[:-1]
        y_peaks_R = h_peaks[1:] + [0]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 

        namestr= f'sinusoid_H{H}L{L}d{delta}k{k}_U{BC.U}_Q{BC.Q}_Re{Re}'
        super().__init__(self.x0, self.xf, y0, yf, N, BC, Re,namestr, x_peaks, y_peaks)


    def h_fun(self,x):
        if x < -self.l or x > self.l:
            if self.k%2 == 0:
                return 0 #yf - hin
            else:
                return self.yf - self.H*(1-self.delta)
            
        else:
            return self.yf -self.H*(1+self.delta * np.cos(self.period * x))
            
        
class TriCavity(PWLinear):
    def __init__ (self, args, BC, Re, N):
        H, L = args
        x0 = 0
        xf = 2
        l=xf//2
        y0 = 0
        yf = H
        x_peaks = [x0, x0+l, xf]
        y_peaks=[[yf,yf],[0,0],[yf,yf]]
        namestr = f"TriCavity_H{H}L{L}_Re{Re}_Q{BC.Q}_U{BC.U}"
        super().__init__(x0, xf, y0, yf, N, BC, Re,namestr, x_peaks, y_peaks)
        
        
class TriSlider(PWLinear):
    def __init__ (self, args, BC, Re, N):

        Hin, H,  Hout, Lin, La, Lb, Lout = args

        x0=0
        xf = Lin+La+Lb+Lout
        x_peaks = [x0,Lin,Lin+La,Lin+La+Lb,xf]
        
        y0 = 0
        yf = max(Hin, Hout, H)
    
        
        y_peaks=[[yf,yf-Hin],[yf-Hin,yf-Hin],[yf-H,yf-H],[yf-Hout,yf-Hout], [yf-Hout, yf]]
        namestr = f"TriCavity_Hin{Hin}H{H}Hout{Hout}La{La}Lb{Lb}Lin{Lin}Lout{Lout}_Re{Re}_Q{BC.Q}_U{BC.U}"
        super().__init__(x0, xf, y0, yf, N, BC, Re,namestr, x_peaks, y_peaks)