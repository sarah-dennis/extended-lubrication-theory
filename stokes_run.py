# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples
import boundary as bc

zoom_on = False   

#------------------------------------------------------------------------------
Example = examples.BFS
h_inlet = 2
h_outlet = 1
l_in = 8
l_out=8
args = [h_inlet, h_outlet, l_in, l_out]

#------------------------------------------------------------------------------
# Example = examples.TriSlider

# l_in = 7
# l_out = l_in
# h_in = 1 
# h_out = h_in
# h_vertex = 2
# l_a = 1.25
# l_b = 0.75
# args =  [h_in, h_vertex, h_out, l_in, l_a, l_b, l_out]


#------------------------------------------------------------------------------
# Example = examples.Logistic

# h_in = 2
# h_out = 1
# L_total = 16
# delta = 8         #slope: delta*(h_in-h_out)/4
# args = [h_in, h_out, L_total, delta]

#------------------------------------------------------------------------------

U=0
Q=1
Re=0

BC = bc.Mixed(U,Q)

#------------------------------------------------------------------------------

Solver = control.Stokes_Solver(Example, args, BC, Re)                

N=80

# Solver.new_run(N) 

# Solver.load_run(N)
# Solver.load_scale(N,2*N) 

# k = 4


# Solver.new_run_many(N, 2, k)  
# Solver.load_run_new_many(N, 2, k)

# Solver.load_run_many(N, 2, k)

Solver.load_plot(N, zoom=zoom_on)







