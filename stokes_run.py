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

max_iters = 500000 
#------------------------------------------------------------------------------

h_in = 2
h_out = 1
l_in = 8
l_out=8
args = [h_in, h_out, l_in, l_out]
Example = examples.BFS

#------------------------------------------------------------------------------

# H=2
# L=16
# delta=1 
# args = [H, L, delta]
# Example = examples.BFS_pwl
#------------------------------------------------------------------------------

# H=2
# h=1
# L=16
# delta = 8 #slope: -delta*(H-h)/4
# args = [H, h, L, delta]
# Example = examples.Logistic
#------------------------------------------------------------------------------

# H = 1
# delta = 1/4
# k = 2 #int 
# l = 1
# L = 8
# args = [H, delta, k,  l, L]
# Example = examples.Sinusoid

#------------------------------------------------------------------------------
# H = 2
# L = 4
# args = [H, L]
# Example = examples.TriCavity

#------------------------------------------------------------------------------

l=7
h=1
H = 2
args =  [h, H, h, l, 1.25, 0.75, l]
Example = examples.TriSlider

#------------------------------------------------------------------------------

U=0
Q=1
Re=0

BC = bc.Mixed(U,Q)

#------------------------------------------------------------------------------

solver = control.Stokes_Solver(Example, args, BC, Re, max_iters=max_iters)                

N=80

       

# solver.new_run(10) 


# solver.load_scale(20,40) 
solver.load_run(40)

# solver.load_run_many(20, 2, 4)

# solver.new_run_many(N, 2, 3)  
# solver.load_run_new_many(N, 2, 3)

solver.load_plot(40, zoom=zoom_on)

# ------------------------------------------------------------------------------
# solver.compare(10,[20,40,80],160)







