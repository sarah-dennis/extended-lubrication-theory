# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWC_Height,  PWL_Height, LogisticHeight
# -----------------------------------------------------------------------------------------------------------------------------------

class BFS(PWC_Height):
    def __init__(self, args, N):
        h, H, l, L = args
        x0 = 0
        xf = L
        N_regions = 2
        x_peaks = np.asarray([0, l, L], float)
        h_peaks = np.asarray([[h, h], [h, H], [H, H]], float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks)

# -----------------------------------------------------------------------------------------------------------------------------------

class TriSlider(PWL_Height):
    def __init__(self, args, N):
        h_in, h, h_out, l_in, l_a, l_b, l_out=args
        N_regions = 4

        x0 = 0
        xf = x0 + l_in + l_a + l_b + l_out
        x_peaks = np.asarray([x0, x0+l_in, x0+l_in+l_a, x0+l_in+l_a+l_b, xf], float)

        h_peaks = np.asarray(([[0, h_in], [h_in, h_in], [h, h], [h_out, h_out], [h_out, 0]]), float)
       
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks)

#-----------------------------------------------------------------------------------------------------------------------------------

class Logistic(LogisticHeight):
    def __init__(self, args, N):
        H_in, H_out, L, delta = args
        x0 = 0
        xf = L
        center = L//2
       
        super().__init__(x0, xf, N, H_in, H_out, center, delta)


