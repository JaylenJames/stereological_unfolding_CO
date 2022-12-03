#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:06:57 2019

@author: JaylenJames

Function used to create frequency distribution stereogram from 2D measurement
    data.
"""

import numpy as np


def stereogram():
    k = np.float64(10.0)
    s = np.float64(10.0)
    
    P = np.zeros([10,10])           #initialize blank matrix to P Matrix
    
    
    index_row = 0
    
    for alpha in range(1, 11):  #cycle through alpha, size domain
        index_col = 0
        
        for i in range(1, 11):   #cycle through beta, shape domain
            
            if alpha == i:
                p_ai = np.sqrt(i-(3/4))
                
            elif alpha > i:
                p_ai = 0
                
            else:
                p_ai = np.sqrt((i-0.5)**2 - (alpha-1)**2) - np.sqrt((i - 0.5)**2 - alpha**2)
                
            P[index_row][index_col] = p_ai
            index_col += 1
            
        index_row += 1
    
    
    
    
    
    #t = 2   #find appropriate values and replace
    #t_1 = 1 #find appropriate values and replace
    
    
    Q = np.zeros([10,10])  #Init blank prolate spheroid matrix
    
    index_row = 0
    
    for j in np.linspace(0.1, 1,10)-0.05:  #cycle through alpha, size domain
        index_col = 0
        
        for beta in np.linspace(0.1, 1,10)-0.05:   #cycle through beta, shape domain
            
            #t_beta = (((2*j)-(2*beta)+1)/((2*k)-(2*j)+1))**0.5
            #t_beta_plus_one = (((2*j)-(2*(beta+1))+1)/((2*k)-(2*j)+1))**0.5
            
            
            def f(t):
                val = t/((t**2.0 - 1.0)  +  np.arctanh(t))  #Update equation correcting for arg tanh
                return val
    
#            def t(b):
#                value = ((((2.0*k) - (2.0*b) + 2.0))/((2.0*j) - (2.0*b) + 2.0))**(0.5)
#                return value
            
            t_1 = ((((2.0*k) - (2.0*1.0) + 2.0))/((2.0*j) - (2.0*beta) + 2.0))**(0.5)
            
            if j == beta:
                
                t_j = ((((2.0*k) - (2.0*j) + 2.0))/((2.0*j) - (2.0*beta) + 2.0))**(0.5)
                
                f_t_j = f(t_j) #t_j/(t_j**2.0 - 1.0) +  np.arctanh(j)
                
                
                p_jb = np.sqrt(t_1**2.0) - f_t_j   #replace t with t_j here
                
            elif j > beta:
                p_jb = 0.0
                
            else:
                
                #f_t_beta = beta/(beta**2.0 - 1.0)  +  np.arctanh(beta)
                #f_t_beta_pone = (beta+1.0)/((beta+1.0)**2.0 - 1.0) + np.arctanh((beta+1.0))
                t_beta = ((2.0*k-(2.0*beta)+2.0)/(2.0*j - (2.0*beta)+1))**(0.5)
                t_beta_pone = ((2.0*k-(2.0*(beta+1.0))+2.0)/(2.0*j - (2.0*(beta+1.0))+1))**(1.0/2.0)
                
                #print(f_t_beta_pone)
                p_jb = np.sqrt(t_1**2 - 1.0)*(f(t_beta) - f(t_beta_pone))
                #p_jb_test = np.sqrt(t_1**2 - 1.0*(f(t_1)) - f(f_t_beta_pone))
                
            Q[index_row][index_col] = p_jb
            index_col += 1
            
        index_row += 1
    
    
    P_inv = np.linalg.inv(P) 
    Q_inv = np.linalg.inv(Q)

    return P_inv, Q_inv


if __name__ == '__main__':
    stereogram()      

#np.arange(1, 11, dtype=np.float64) #way to create a list of values of a
    # certaiin data type
    
    
    
    
    

    