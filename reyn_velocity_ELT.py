# -*- coding: utf-8 -*-
"""
Created on Thu May  8 14:11:23 2025

@author: sarah
"""
import numpy as np
import domain as dm

def make_ELT_velocity(height, BC, adj_pressure):

    # derivatives of the reynolds pressure
    pxs = adj_pressure.reyn_pxs
    p2xs = adj_pressure.reyn_p2xs
    p3xs = adj_pressure.reyn_p3xs
    p4xs = adj_pressure.reyn_p4xs
    
    # sigmas = adj_pressure.sigmas
    sigma_xs = adj_pressure.sigma_xs
    sigma_2xs = adj_pressure.sigma_2xs

    us = make_us(height, BC.U, pxs, p2xs, p3xs, sigma_xs)
    vs = make_vs(height,  BC.U, pxs, p2xs, p3xs, p4xs, sigma_xs, sigma_2xs)
    
    for i in height.i_peaks[1:-1]:
        for j in range(height.Ny):
            us[j,i-1: i+2] = dm.avg_x(us[j,i-2:i+3])
            vs[j,i-1: i+2] = dm.avg_x(vs[j,i-2:i+3])
  
    return us, vs


def make_us(height, U, pxs, p2xs, p3xs, sigmaxs):
    us = np.zeros((height.Ny, height.Nx))

    visc = 1#height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        h2x = height.h2xs[i]
        
        px = pxs[i]
        p2x = p2xs[i]
        p3x = p3xs[i]
        sigmax = sigmaxs[i]
        
        h2 = h**2
        h3 = h**3
        hx2 = hx**2
        for j in range(height.Ny):
            y = height.ys[j]
            y2 = y**2
            y3 = y**3
            y4 = y**4
            
            if y > h:
                continue
            else:
              # reynolds u
                uA = 1/(2*visc)  * px                   * (y2 - h*y) + U*(1 - y/h)
                
              # adjusted u
                uB = -1/(24*visc)* p3x                  * (y4 - h3*y)
                uC = 1/(12*visc) * (h*p3x + 2*p2x*hx + px*h2x)  * (y3 - h2*y)
                uD = -U/6        * (-2*hx2/h3 + h2x/h2) * (y3 - h2*y)
                uE = (1/2)       * sigmax               * (y2 - h*y)
                
                us[j,i] = (uA+uB+uC+uD+uE)

    return us


def make_vs(height, U, pxs, p2xs, p3xs, p4xs, sigmaxs, sigma2xs):
    vs = np.zeros((height.Ny, height.Nx))
    visc = 1# height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        h2x = height.h2xs[i]
        h3x = height.h3xs[i]
        px = pxs[i]
        p2x = p2xs[i]
        p3x = p3xs[i]
        p4x =  p4xs[i]
        sigmax = sigmaxs[i]
        sigma2x = sigma2xs[i]
        
        h2 = h**2
        h3 = h**3
        h4 = h**4
        hx3 = hx**3
        for j in range(height.Ny):
            y = height.ys[j]
            y2 = y**2
            y3 = y**3
            y4 = y**4
            y5 = y**5
            
            if y > h:
                continue
            
            else:
                
                # v(x,y) = integral -du/dx  dy
                vA = -1/(2*visc)* (p2x*(y3/3-h*y2/2) - px*hx*y2/2) - U*hx*y2/(2*h2) 
                vB = 1/(24*visc)*(p4x*(y5/5-h3*y2/2) - 3*h2*y2*p3x*hx/2 )
                vCa = (h*p4x + 3*p3x*hx + 3*p2x*h2x + px*h3x)*(y4/4 - h2*y2/2)
                vCb = -(h*p3x + 2*p2x*hx + px*h2x)*h*hx*y2
                vC = -1/(12*visc)*(vCa+vCb)
                vDa = (y4-2*h2*y2)/(4*h2)*h3x + (3*y4-2*h2*y2)/(2*h4)*hx3 + (4*h2*y2-3*y4)/(2*h3)*hx*h2x

                vD = U/6 * vDa
                vE = -1/2 * (sigma2x*(y3/3 -h*y2/2) - sigmax*hx*y2/2)
                vs[j,i] = vA+vB+vC+vD+vE

    return vs

    