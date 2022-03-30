#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:12:42 2021

@author: willy
"""
import constances as csts

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

#______________________________________________________________________________
#       compute q_plus_i, q_minus_i, phi_EPO_plus, phi_EPO_minus : debut
#______________________________________________________________________________
def compute_qi_plus_minus(Pi, Ci, Si, Si_max):
    """
    compute q_plus et q_minus for player_i

    """
    q_minus_k_i = diff_positive(Ci, Pi) - diff_positive(Pi, Ci+Si_max-Si)
    q_plus_k_i = diff_positive(Pi, Ci) - diff_positive(Ci, Pi+Si)
                    
    return q_plus_k_i, q_minus_k_i
        
def compute_phi_EPO_plus(q_plus_k, quantity, a):
    """
    compute the benefit of energy sold by SG to EPO
    """
    return quantity * pow(q_plus_k, a)
        
def compute_phi_EPO_minus(q_minus_k, quantity, b):
    """
    compute the cost of energy bought by SG to EPO
    """
    return quantity * pow(q_minus_k, b)
    
#______________________________________________________________________________
#       compute q_plus_i, q_minus_i, phi_EPO_plus, phi_EPO_minus : fin
#______________________________________________________________________________

#______________________________________________________________________________
#          compute 
#           beta_{sg,t-1}^{+}, beta_{sg,t-1}^{-}, 
#           pi_0_plus, pi_0_minus: 
#           Debut
#______________________________________________________________________________
def compute_beta_sg(t_period, dico_chosen_strats_t, 
                    beta_sg_0_minus, beta_sg_0_plus,
                    quantity_a, a, quantity_b, b):
    """
    dico_chosen_t_strats = {"t_0":{"player_0":{ ... } }}
    """
    if t_period == csts.PERIOD_ROOT + str(t_period.split("_")[1]):
        return beta_sg_0_minus, beta_sg_0_plus
    else:
        sum_T_diff_plus, sum_T_diff_minus = 0, 0
        sum_T_prod, sum_T_cons = 0, 0
        for t in range(0, int(t_period.split("_")[1] )):
            t_ = csts.PERIOD_ROOT + str(t)
            sum_N_cons_i = 0; sum_N_prod_i = 0
            for player_i, dico_values in dico_chosen_strats_t[t_].items():
                prod_i = dico_chosen_strats_t[t_][player_i]["prod_i"]
                cons_i = dico_chosen_strats_t[t_][player_i]["cons_i"]
                
                sum_N_cons_i += cons_i
                sum_N_prod_i += prod_i
               
            # 
            diff_minus_t = 0 if sum_N_cons_i < sum_N_prod_i \
                               else sum_N_cons_i - sum_N_prod_i
            diff_plus_t = 0 if sum_N_cons_i > sum_N_prod_i \
                               else sum_N_prod_i - sum_N_cons_i
             
            #
            sum_T_diff_plus += diff_plus_t
            sum_T_diff_minus += diff_minus_t
           
            ##
            sum_T_prod += sum_N_prod_i;  sum_T_cons += sum_N_cons_i
           
        
        beta_sg_t_minus_1_plus =  compute_phi_EPO_plus(
                                    q_plus_k=sum_T_diff_plus, 
                                    quantity=quantity_a, a=a)
        beta_sg_t_minus_1_minus = compute_phi_EPO_minus(
                                    q_minus_k=sum_T_diff_minus, 
                                    quantity=quantity_b, b=b)
        
        return beta_sg_t_minus_1_plus, beta_sg_t_minus_1_minus
    
def compute_pi_0(t_period, beta_sg_t_minus_1_plus, beta_sg_t_minus_1_minus,
                 pi_EPO_plus, pi_EPO_minus):
    pi_0_minus = beta_sg_t_minus_1_minus
    pi_0_plus = beta_sg_t_minus_1_minus * pi_EPO_plus / pi_EPO_minus
    pi_0_plus = round(pi_0_plus, csts.ARRONDI)
    
    return pi_0_plus, pi_0_minus
    
def compute_b0_c0(pi_0_plus, pi_0_minus, Out_sg_k, In_sg_k, 
                  quantity_a, a, quantity_b, b):
    """
    
    """
    b0, c0 = None, None
    if In_sg_k >= Out_sg_k:
        c0 = pi_0_minus
        phi_EPO_plus = compute_phi_EPO_plus(q_plus_k= In_sg_k-Out_sg_k, 
                                            quantity=quantity_a, a=a)
        b0 = (Out_sg_k * pi_0_plus + phi_EPO_plus )/ In_sg_k
    else:
        b0 = pi_0_plus
        phi_EPO_minus = compute_phi_EPO_minus(q_minus_k= Out_sg_k-In_sg_k, 
                                              quantity=quantity_b, 
                                              b=b)
        c0 = (phi_EPO_minus + In_sg_k * pi_0_minus)/Out_sg_k
        
    b0 = round(b0, csts.ARRONDI)
    c0 = round(c0, csts.ARRONDI)
    return b0, c0
#______________________________________________________________________________
#          compute 
#           beta_{sg,t-1}^{+}, beta_{sg,t-1}^{-}, 
#           pi_0_plus, pi_0_minus,
#           b0, c0 
#           fin
#______________________________________________________________________________


