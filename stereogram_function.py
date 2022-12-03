#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:06:57 2019

@author: JaylenJames

Function used to calculate the size and shape correction factors as described
    in Cruz-Orive's paper: Particle size-shape distributions: the general 
    spheroid problem Part II. Stochastic model and practical guide.
    

    
"""

import numpy as np


def stereogram(bins_per_var=10):
    
    bpv = bins_per_var
    k = np.float64(bpv) #
    
    
    
    P = np.zeros([bpv,bpv])           #initialize blank matrix to P Matrix
    
    
    index_row = 0
    
    for alpha in range(1, bpv+1):  #cycle through alpha, size domain
        index_col = 0
        
        for i in range(1, bpv+1):   #cycle through beta, shape domain
            
            if alpha == i:
                p_ai = np.sqrt(i-(3/4))
                
            elif alpha > i:
                p_ai = 0
                
            else:
                p_ai = np.sqrt((i-0.5)**2 - (alpha-1)**2) - np.sqrt((i - 0.5)**2 - alpha**2)
                
            P[index_row][index_col] = p_ai
            index_col += 1
            
        index_row += 1
    
    
    
    Q = np.zeros([bpv,bpv])  #Initialize blank prolate spheroid matrix
    
    
    index_row = 0
    
    #for j in np.linspace(0.1, 1.0, 10)-0.05:  #cycle through alpha, size domain #this is a float64
    for j in np.linspace(1, bpv, bpv):
        index_col = 0
        
        #for beta in np.linspace(0.1, 1.0, 10)-0.05:   #cycle through beta, shape domain
        for beta in np.linspace(1, bpv, bpv):    
            
            
            def f(t):
                val = t/((t**2.0 - 1.0) )  + ((0.5)*np.log((t+1.0)/(t-1.0)))
                return val
    
            
            t_1 = ((((2.0*k) - (2.0*1.0) + 2.0))/((2.0*j) - (2.0*1.0) + 1.0))**(0.5)
            
            if j == beta:
                
                t_j = ((((2.0*k) - (2.0*j) + 2.0))/((2.0*j) - (2.0*beta) + 1.0))**(0.5)
                
                f_t_j = f(t_j)
                
                
                p_jb = np.sqrt(t_1**2.0 - 1.0)*f_t_j   
                
            elif beta > j: 
                p_jb = 0.0
                
            else:
                
                t_beta = ((2.0*k-(2.0*beta)+2.0)/(2.0*j - (2.0*beta)+1.0))**(0.5)
                t_beta_pone = ((2.0*k-(2.0*(beta+1))+2.0)/(2.0*j - (2.0*(beta+1))+1.0))**(0.5)
                
                
                p_jb = np.sqrt((t_1**2.0) - 1.0)*(f(t_beta) - f(t_beta_pone))
                
                
            Q[index_row][index_col] = p_jb
        
            index_col += 1
            
        index_row += 1
    
    P_inv = np.linalg.inv(P) 
    Q_inv = np.linalg.inv(Q)

    
    #Check if matrix inverses are correct
    P_ID_1 = np.dot(P,P_inv)
    P_ID_1[np.abs(P_ID_1)<=1E-13] = 0 #Change very small numbers to 0
    
    Q_ID_1 = np.dot(Q,Q_inv)
    Q_ID_1[np.abs(Q_ID_1)<=1E-13] = 0 #Change very small numbers to 0
    
    
    return P_inv, Q_inv


if __name__ == '__main__':
    P_inverse, Q_inverse = stereogram(10)      

    

    