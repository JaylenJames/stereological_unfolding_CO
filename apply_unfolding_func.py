#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 14:49:33 2020

@author: JaylenJames


This function will apply Cruz-Orive's unfolding algorthm to a matrix
"""
import numpy as np
from stereogram_function import stereogram



def apply_unfolding_func(two_dim_frequencies, zero_negatives = False, bins_per_var = 10):
    """
    Return a matrix describing the frequency of 3D prolate spheroids given the
    frequencies of 2D measurments.
    
    Parameters
    ----------
    two_dim_frequencies : n by n matrix where n is the number of rows and columns.
        The values in the matrix are the frequencies associated 
        with the m and y^2 measurments from a section of a population of 
        spheroidal particles.
    dtype : float
        
    zero_negatives : bool, optional
    If this is set to True, any frequency values calculated in the unfolded
    distribution that are negative will be replaced with 0. 
        
        
    
    Returns
    -------
    out : ndarray
        Array of ones with the same shape and type as `a`.
        
    
    Examples
    --------
   
    
   
    
    """
    B = 1
    s = bins_per_var
    k = 10
    delta = B/s
    P, Q = stereogram(bins_per_var = s) #Obtain size and shape factor values. These are the inverses.
    
    
    H_bar = 90/(80*np.pi) #estimate for original data. It should be adjusted for your data set so that g_ij_sum = 1
                          #     See Cruz-Orive's paper Particle size-shape distributions: the general spheroid prolem Part I:
                          #     Mathematical Model  
    
    
    g = np.zeros([s,s])           #initialize blank g matrix
    val_mat = np.zeros([s,s])
        
    i_index_row = 0
    
    for i in range(1, s+1):  #cycle through row elements of g matrix
        j_index_col = 0
        
        for j in range(1, s+1):   #cycle through column elements of g matrix
            alpha_index_row = 0
            
            for alpha in range(1,s+1): 
                beta_index_col = 0
                
                for beta in range(1,s+1):
                    
                    val_mat[alpha_index_row][beta_index_col] = P[i_index_row][alpha_index_row]*two_dim_frequencies[alpha_index_row][beta_index_col]*Q[beta_index_col][j_index_col]
                
                    beta_index_col += 1
                
                alpha_index_row += 1
                
            val_mat_sum = np.sum(val_mat)
            mltd_val_mat = (H_bar/delta)*val_mat_sum
            g[i_index_row][j_index_col] = mltd_val_mat
            
            j_index_col += 1
            
        i_index_row += 1
    
    
    g_ij_sum = np.sum(g)
    
    
    
    if zero_negatives == True:
        g = g.clip(min=0)
    
    
    return g

if __name__ == '__main__':
    apply_unfolding_func()