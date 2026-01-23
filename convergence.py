# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import readwrite as rw
import stokes_pressure as pressure
#-------------------------------------------------------------------------------
def stokes_cnvg_self(Ex, args, BC, Re, N_min, Ns, N_max, p_err=True):
    ex_min = Ex(args, BC, Re, N_min)
    ex_max = Ex(args, BC, Re, N_max)
    
    # Load max grid for 'true' values
    u_max, v_max, psi_max, past_iters = rw.read_stokes(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
    psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))

    if p_err:
        p_max = pressure.pressure(ex_max, u_max, v_max)
        dp_max, res_max = pressure.resistance(ex_max, p_max) 
        p_max = p_max.reshape((ex_max.Ny,ex_max.Nx))
        
    mult_max = int(N_max/N_min)
    
    err = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_inf = np.zeros(len(Ns)+1)
    err_l1 = np.zeros(len(Ns)+1)
    err_l2 = np.zeros(len(Ns)+1)
    

    if p_err:
        err_p = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
        err_p_inf = np.zeros(len(Ns)+1)
        err_p_l1 = np.zeros(len(Ns)+1)
        err_p_l2 = np.zeros(len(Ns)+1)
        err_dp = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            ex_n = Ex(args, BC, Re, N_min)
            mult = 1
        else:
            ex_n = Ex(args, BC, Re, Ns[n-1])
            mult = int(Ns[n-1]/N_min)
            
            
        u_n, v_n, psi_n, past_iters = rw.read_stokes(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
        psi_n=psi_n.reshape((ex_n.Ny, ex_n.Nx))

        if p_err:
            p_n = pressure.pressure(ex_n, u_n, v_n)
            dp_n, res_n = pressure.resistance(ex_n, p_n)
            p_n=p_n.reshape((ex_n.Ny,ex_n.Nx))
            
        for k_min in range(ex_min.Ny*ex_min.Nx):
        # all indices (i,j) on grid N_min
            i = k_min % ex_min.Nx
            j = (k_min // ex_min.Nx)
            
            # -> indices on grid N || N_min
            i_n = mult * i
            j_n = mult * j

            # -> indices on grid N_max || N_min
            i_max = mult_max * i
            j_max = mult_max * j
            
            err[n,j,i] = abs(psi_max[j_max,i_max] - psi_n[j_n,i_n])
            err_l1[n] += err[n,j,i]
            err_l2[n] += err[n,j,i]**2
            
            if p_err:
                err_p[n,j,i] = abs(p_max[j_max,i_max] - p_n[j_n,i_n])
                err_p_l1[n]+= err_p[n,j,i] 
                err_p_l2[n]+= (err_p[n,j,i]**2)
                
        size_n= ex_n.Nx*ex_n.Ny   
        
        err_inf[n] = np.max(err[n])
        err_l1[n] /= size_n
        err_l2[n] = np.sqrt(err_l2[n] / size_n)
        
        if p_err:
            err_p_inf[n] = np.max(err_p[n])
            err_p_l2[n] = np.sqrt(err_p_l2[n])
            err_dp[n]=np.abs(dp_max-dp_n)
            
    l1_rate = convg_rate(err_l1)
    l2_rate = convg_rate(err_l2)
    inf_rate = convg_rate(err_inf)
    cnvg_rates = np.stack([l1_rate, l2_rate, inf_rate], axis=0)
    
    if p_err:
        p_l1_rate = convg_rate(err_p_l1)
        p_l2_rate = convg_rate(err_p_l2)
        p_inf_rate = convg_rate(err_p_inf)
        dp_rate = convg_rate(err_dp)
        p_cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate, dp_rate], axis=0)
        
    print("stream cnvg rates")
    print("l1: " + np.array2string(cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(cnvg_rates[1], precision=2))
    print("linf" + np.array2string(cnvg_rates[2], precision=2))
    
    if p_err:
        print("pressure cnvg rates")
        print("l1: " + np.array2string(p_cnvg_rates[0], precision=2))
        print("l2: " + np.array2string(p_cnvg_rates[1], precision=2))
        print("linf" + np.array2string(p_cnvg_rates[2], precision=2))
        print("dp: " + np.array2string(p_cnvg_rates[3], precision=2))

    return err_l1, err_l2, err_inf, cnvg_rates, ex_min

#------------------------------------------------------------------------------
    
def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates

#------------------------------------------------------------------------------


