#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 15:41:51 2025

@author: sarahdennis
"""
import reyn_control
import reyn_examples
import boundary as bc
import stokes_control
import stokes_examples

import graphics

import numpy as np
#------------------------------------------------------------------------------
# Errors 
def linf(ax, ay, bx, by):
    return np.max((np.max(np.abs(ax-bx)), np.max(np.abs(ay-by))))

def linf_(ax, ay, bx, by):
    norm_x = np.max(np.abs(ax-bx))
    norm_y = np.max(np.abs(ay-by))
    return np.max((norm_x, norm_y))

def l1(ax,ay,bx,by):
    return np.sum(np.abs(ax-bx)) + np.sum(np.abs(ay-by))

def l2(ax,ay,bx,by):
    return np.sum((ax-bx)**2 + (ay-by)**2) **(1/2)

#----------------
plots_on =  not True # plots p(x,y) contour-mesh and (u,v) streamlines
zoom_on = False # plot a zoomed frame (set params in control)
scaled_on=False # plot on scaled axis (set params in control)

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =0


Q=1


BC = bc.Mixed(U, Q)

Re=0


N = 80 # grid size |1|= N

#------------------------------------------------------------------------------
Reyn_Example = reyn_examples.Logistic
Stokes_Example= stokes_examples.Logistic
h_in = 2   # inlet height
h_out= 1   # outlet height
l = 16   #  length

tests = [2, 3, 4, 6, 8, 16, 32]
test_args = [[h_in , h_out, l, lam] for lam in tests]
exstr = 'Logistic Step'
label = '$\lambda$'
#------------------------------------------------------------------------------
# Reyn_Example = reyn_examples.TriSlider
# Stokes_Example = stokes_examples.TriSlider

# h_in=1  # inlet height
# h0=1/4   # apex height 
# h_out = 1  #oulet height
# l_in = 7  # inlet length
# l_out = 7  #outlet length
# l_a = 1.25  # base length A  
# l_b = 0.75  # base length B 

# #test h0

# tests = [1/16, 1/8, 1/4, 1/2, 3/4, 5/4, 3/2, 7/4, 2]
# test_args = [[h_in, h0, h_out, l_in, l_a, l_b, l_out] for h0 in tests]

# exstr = 'Triangular Slider'
# label = '$H_v$'
# #------------------------------------------------------------------------------
# Reyn_Example = reyn_examples.BFS
# Stokes_Example = stokes_examples.BFS

# h_in=1  # inlet height
# # h_out=2   # outlet height 
# l_in = 2  # inlet length
# l_out = 2  #outlet length

# #test h_out
# tests = [2]
# test_args = [[h_in, h_out, l_in, l_in+l_out] for h_out in tests]

# exstr = 'BFS'
# label = '$H$'
 #------------------------------------------------------------------------------
k = 0
num_tests=len(test_args)

fun_labels= ['Reyn',  '$\epsilon^2$-PLT', '$\epsilon^4$-PLT','VA-ELT', 'TG-ELT']
num_models = 5 #reyn, VA-TG adj, e2 pert, e4 pert

l1_V_errs= np.zeros((num_models,num_tests))
linf_V_errs = np.zeros((num_models,num_tests))
l2_V_errs = np.zeros((num_models,num_tests))


l1_P_errs = np.zeros((num_models,num_tests))
linf_P_errs = np.zeros((num_models,num_tests))
l2_P_errs = np.zeros((num_models,num_tests))

for args in test_args:
#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
    reyn_solver = reyn_control.Reynolds_Solver(Reyn_Example, BC, args)
    reyn_P, reyn_V = reyn_solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)
    reyn_ps,reyn_us, reyn_vs = np.nan_to_num(reyn_P.ps_2D),reyn_V.u,reyn_V.v
    
    adj_P, adj_V = reyn_solver.fd_VA_ELT_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)
    adj_ps,adj_us, adj_vs = np.nan_to_num(adj_P.ps_2D),adj_V.u,adj_V.v
    
    adj_TG_P, adj_TG_V = reyn_solver.fd_TG_ELT_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)
    adj_TG_ps,adj_TG_us,adj_TG_vs = np.nan_to_num(adj_TG_P.ps_2D),adj_TG_V.u,adj_TG_V.v

    pert = reyn_solver.fd_pert_solve(N, order=4, plot=plots_on, scaled=scaled_on, zoom=zoom_on,get_all=True)
    e2_ps, e2_us, e2_vs =  np.nan_to_num(pert.pert2_pressure.ps_2D), pert.pert2_velocity.u, pert.pert2_velocity.v

    e4_ps, e4_us, e4_vs =  np.nan_to_num(pert.pert4_pressure.ps_2D), pert.pert4_velocity.u, pert.pert4_velocity.v

#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
    
    stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, BC, Re)
    stokes_ps, stokes_us, stokes_vs = stokes_solver.load(N)
    stokes_ps = np.nan_to_num(stokes_ps)

    if plots_on:
        stokes_solver.load_plot(N)   
    #------------------------------------------------------------------------------
    l1_stokes_V = l1(stokes_us, stokes_vs, 0, 0)
    linf_stokes_V = linf(stokes_us, stokes_vs, 0, 0)
    l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)

    l1_stokes_P = l1(stokes_ps, 0, 0, 0)
    linf_stokes_P = linf(stokes_ps, 0, 0, 0)
    l2_stokes_P = l2(stokes_ps, 0, 0, 0)
    #------------------------------------------------------------------------------
    # fun_labels= ['Reyn',  '$\epsilon^2$-PLT', '$\epsilon^4$-PLT','VA-ELT', 'TG-ELT']
    test_us = [[reyn_us, reyn_vs], [e2_us,e2_vs], [e4_us,e4_vs], [adj_us, adj_vs],  [adj_TG_us, adj_TG_vs]]
    test_ps = [reyn_ps, e2_ps, e4_ps, adj_ps, adj_TG_ps]
    
    for i in range(len(test_us)):
        # l1_V_errs[i,k] = l1(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/l1_stokes_V *100
        l2_V_errs[i,k] = l2(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/l2_stokes_V *100
        # linf_V_errs[i,k] = linf(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/linf_stokes_V *100
    
    for i in range(len(test_ps)):
    
        # l1_P_errs[i,k] = l1(stokes_ps, 0, test_ps[i], 0)/l1_stokes_P *100
        l2_P_errs[i,k] = l2(stokes_ps, 0, test_ps[i], 0)/l2_stokes_P *100
        # linf_P_errs[i,k] = linf(stokes_ps, 0, test_ps[i], 0)/linf_stokes_P *100
 
    k+=1
    
    
# graphics.plot_log_multi(l1_V_errs[L:-1], tests, f'$L_1$ rel. %-error Velocity, {exstr} $Q=${Q:.1f}', fun_labels, [label, '$L_1$ rel. %-error'],loc='left')
# graphics.plot_log_multi(l2_V_errs[:-1], tests, f'Velocity $L_2$ rel. %-error, {exstr}', fun_labels, [label, 'Velocity $L_2$ rel. %-error'],loc='center')#,loc='left')
graphics.plot_2D_multi(l2_V_errs[:-1], tests, f'Velocity $L_2$ rel. %-error, {exstr}', fun_labels, [label, 'Velocity $L_2$ rel. %-error'],loc='right')#,loc='right')
# graphics.plot_log_multi(linf_V_errs, tests, f'$L_\infty$ rel. %-error Velocity, {exstr} $Q=${Q:.1f}',  fun_labels,  [label, '$L_\infty$ rel. %-error'],loc='left')


# graphics.plot_log_multi(l1_P_errs, tests, f'$L_1$ rel. %-error Pressure, {exstr} $Q=${Q:.1f}', fun_labels, [label, '$L_1$ rel. %-error'],loc='left')
# graphics.plot_log_multi(l2_P_errs, tests, f'Pressure $L_2$ rel. %-error, {exstr}',  fun_labels, [label, 'Presure $L_2$ rel. %-error '],loc='right')#,loc='left')
graphics.plot_2D_multi(l2_P_errs, tests, f'Pressure $L_2$ rel. %-error, {exstr}',  fun_labels, [label, 'Pressure $L_2$ rel. %-error '],loc='left')#,loc='lower')

# graphics.plot_log_multi(linf_P_errs, tests, f'$L_\infty$ rel. %-error Pressure, {exstr} $Q=${Q:.1f}',  fun_labels,  [label, '$L_\infty$ rel. %-error'],loc='left')


