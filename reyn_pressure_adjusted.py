# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:16:58 2025

@author: sarah
"""
import numpy as np
import domain as dm
import boundary as bc
import graphics
import reyn_pressure_finDiff as fd

def make_adj_ps(height, BC, reyn_ps, TG=False):
    ps_adj = np.zeros((height.Ny, height.Nx))

    hs = height.hs
    hxs = height.hxs
    pxs = dm.center_diff(reyn_ps, height.Nx, height.dx)
    p2xs = dm.center_second_diff(reyn_ps, height.Nx, height.dx)
    p3xs = dm.center_third_diff(reyn_ps, height.Nx, height.dx)
    p4xs = dm.center_fourth_diff(reyn_ps, height.Nx, height.dx)
    
    if TG:
        sigmas, sigma_xs, sigma_2xs = 0, 0, 0 
    else: #VA
        sigmas, sigma_xs, sigma_2xs = make_sigmas(height,BC, pxs,p2xs,p3xs,p4xs)

    #---------------------------------------------------------------------------

    for i in range(height.Nx): 
        h = hs[i]
        hx = hxs[i]
        px = pxs[i]
        pxx = p2xs[i]

        phi1x = -(pxx*h + px*hx)/2 + BC.U/(h**2)*hx #*visc

        for j in range(height.Ny):
            y = height.ys[j]
            if y > h:
                ps_adj[j,i] = None
                continue

            else:  
                adj = -pxx*(y**2)/2 - phi1x*y 
                if TG: 
                    ps_adj[j,i] = reyn_ps[i] + adj
                    
                else:
                    ps_adj[j,i] = reyn_ps[i] + adj + sigmas[i] #*visc 


    sigma_derivs = [sigmas, sigma_xs, sigma_2xs]
    reyn_derivs = [pxs, p2xs, p3xs, p4xs]
    
    return ps_adj, reyn_derivs, sigma_derivs

def adj_rhs(height, BC, pxs, p2xs, p3xs, p4xs):
    vs = np.zeros(height.Nx)
   
    for i in range(1,height.Nx-1):
        h = height.hs[i]
        hx = height.hxs[i]
        h2x = height.h2xs[i]
        h3x = height.h3xs[i]
        px = pxs[i]
        p2x = p2xs[i]
        p3x = p3xs[i]
        p4x = p4xs[i]

        v_a = (h**5)*p4x + 5*(h**4)*p3x*hx #/visc
        v_b = (h*p4x+ 3*hx*p3x + h3x*px + 3*h2x*p2x)*(h**4) #/visc
        v_c = 4*(h*p3x + 2*hx*p2x + h2x*px)*(h**3)*hx#/visc
        v_d = ((h**2)*h3x-2*(hx**3)-2*h*hx*h2x)
    
        vs[i] = 3*v_a/20 - (v_b + v_c)/4 + v_d*BC.U/2 

    vs[0] = 0
    vs[-1] = 0
    return vs
 
def make_sigmas(height, BC, pxs, p2xs, p3xs, p4xs):
    s = np.zeros( height.Nx)
    sx = np.zeros( height.Nx)
    sxx = np.zeros( height.Nx)
    
    if isinstance(BC, bc.Fixed): # match reyn dP
        M = fd.make_mat(height, BC) 
        rhs = adj_rhs(height, BC, pxs, p2xs, p3xs, p4xs) 
        
        s  = np.linalg.solve(M, rhs)
        
        sx = dm.center_diff(s, height.Nx, height.dx)
        sxx = dm.center_second_diff(s, height.Nx, height.dx)
                    
    elif isinstance(BC, bc.Mixed): #match reyn Flux
      
        for i in range(height.Nx):
            h = height.hs[i]
            hx = height.hxs[i]
            h2x = height.h2xs[i]
            h3x = height.h3xs[i]
            p4x = p4xs[i]
            px = pxs[i]
            p2x = p2xs[i]
            p3x = p3xs[i]

            sx_A = (3/20)*(h**2) * p3x - (1/4)*h * (h*p3x + 2*p2x*hx + px*h2x)
            sx_B = BC.U/2 * (-2*(hx**2)/(h**2)+ h2x/h)
            sx[i] = sx_A + sx_B 
            
            s2x_A = 3/20*(2*p3x*h*hx+p4x*(h**2))-1/4*(h*p4x+3*p2x*h2x+3*p3x*hx+px*h3x)*h 
            s2x_B = -1/4*(h*p3x+2*p2x*hx+px*h2x)*hx  
            s2x_C = BC.U/2*(4/(h**3)*(hx**3)-5/(h**2)*hx*h2x+1/h*h3x)
            sxx[i] = s2x_A + s2x_B +s2x_C
    

        for i in range(height.Nx):
            if i > 1:
                s[i] = (4*s[i-1] -s[i-2] + 2*height.dx*sx[i])/3
            elif i > 0:
                s[i] = s[i-1] + sx[i]*height.dx
                
    s -= s[-1]  #b.c. p(xf)=0 --> s(xf)=0
          
    return s, sx, sxx   
