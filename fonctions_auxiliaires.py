#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:12:42 2021

@author: willy
"""


#______________________________________________________________________________
#                       fonctions annexes: debut 
#______________________________________________________________________________
def diff_positive(op1, op2):
    """
    difference positive between op1 and op2.
    op1, op2 : float
    return diff 
    mean that :
        diff = 0 if op1 < op2
        diff = op1 - op2 if op1 >= op2
    """
    
    return 0 if op1 < op2 else op1 - op2 

#______________________________________________________________________________
#                       fonctions annexes: fin 
#______________________________________________________________________________