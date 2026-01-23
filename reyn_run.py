# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import boundary as bc
import reyn_examples as examples
import graphics

#-------------------plotting---------------------------------------------------

plots_on = True
uv_on =  not True # plot u(x,y) & v(x,y)
zoom_on =   False    # plot a zoomed-in window, set location in reyn_control.py
scaled_on = False  # plot in scaled variables x/X, y/Y etc.

#------------------------------------------------------------------------------
# Examples
#------------------------------------------------------------------------------
# Example = examples.BFS
# H=1
# h=2
# l=8
# L=16
# args =  [h, H, l, L]
#-----------------------------------------------------------------------------
Example = examples.TriSlider
h_in=1
h=2
h_out = h_in
l_in = 1
l_out = 1
l_a = 1.25
l_b = 0.75
args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]
#-----------------------------------------------------------------------------
# Example = examples.Logistic
# delta = 8 # max slope: delta*(H-h)/4
# H = 2     # outlet height
# h = 1       # inlet height
# L = 16     # total length
# args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 0
# BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1

# U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = control.Reynolds_Solver(Example, BC, args)

N = 320

#------------------------------------------------------------------------------

solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)

solver.fd_TG_ELT_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)

solver.fd_VA_ELT_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)

solver.fd_pert_solve(N, order=4,  plot=plots_on, scaled=scaled_on, zoom=zoom_on)
