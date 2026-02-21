# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import boundary as bc
import reyn_examples as examples

#-------------------plotting---------------------------------------------------
plots_on = True
zoom_on =  not False    # plot a zoomed-in window, set window position in reyn_control.py

#------------------------------------------------------------------------------
# Examples
#------------------------------------------------------------------------------

# Example = examples.BFS
# h_inlet =2
# h_outlet =1
# l_inlet =8
# L_total =16
# args =  [h_inlet, h_outlet, l_inlet, L_total]

#-----------------------------------------------------------------------------
# Example = examples.TriSlider
# h_in=1
# h_vertex=2
# h_out = h_in
# l_in = 7
# l_out = l_in
# l_a = 1.25
# l_b = 0.75
# args =  [h_in, h_vertex, h_out, l_in, l_a, l_b, l_out]

#-----------------------------------------------------------------------------
Example = examples.Logistic
lam = 8     # max slope: lam*(h_in-hout)/4
h_in = 2     
h_out = 1    
L_total = 16     
args = [ h_in, h_out, L_total, lam]

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------
# {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

# {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

Solver = control.Reynolds_Solver(Example, BC, args)

N = 160

#------------------------------------------------------------------------------

Solver.fd_solve(N, plot=plots_on, zoom=zoom_on)

Solver.fd_TG_ELT_solve(N, plot=plots_on, zoom=zoom_on)

Solver.fd_VA_ELT_solve(N, plot=plots_on, zoom=zoom_on)

Solver.fd_pert_solve(N, order=2,  plot=plots_on, zoom=zoom_on)

Solver.fd_pert_solve(N, order=4,  plot=plots_on, zoom=zoom_on)
