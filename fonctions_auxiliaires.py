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

#______________________________________________________________________________
#          compute 
#           S_minus_i, S_plus_i, I_m, I_M, O_m, O_M 
#           ben_i, cst_i 
#           debut
#______________________________________________________________________________
def compute_Si_plusminus_ImM_OmM(state_i, Si_max, Si, Ci, Pi):
    """
    compute S_plus_i, S_minus_i and also Ii_m, Ii_M, Oi_m, Oi_M 
    
    """
    S_plus_i, S_minus_i = None, None
    Ii_m, Ii_M, Oi_m, Oi_M = None, None, None, None
    if state_i == csts.STATES[0]:
        # Deficit
        S_plus_i = Si; S_minus_i = 0
        Oi_m = Ci - (Pi + Si)
        Oi_M = Ci - Pi if Oi_M is None else Oi_M + (Ci - Pi)
    elif state_i == csts.STATES[1]:
        # Self
        S_plus_i = Si; S_minus_i = Si - (Ci - Pi);
        Oi_M = Ci - Pi if Oi_M is None else Oi_M + (Ci - Pi)
    elif state_i == csts.STATES[2]:
        # Surplus
        S_minus_i = Si; S_plus_i = max(Si_max, Si + (Ci - Pi));
        Ii_m = diff_positive(op1=Pi, op2=Ci+(Si_max-Si))
        Ii_M = Pi - Ci
    else:
        print(f" ---> state={state_i} DON'T EXIST <---")
        
    return S_plus_i, S_minus_i, Ii_m, Ii_M, Oi_m, Oi_M

#______________________________________________________________________________
#          compute 
#           S_minus_i, S_plus_i, I_m, I_M, O_m, O_M 
#           ben_i, cst_i 
#           fin
#______________________________________________________________________________


#______________________________________________________________________________
#          compute u_i
#           debut
#______________________________________________________________________________  
def compute_utility_fonction_ui(dico_chosen_strats_k, 
                                t_period):
    """
    compute utility fonction u_i
    dico_chosen_strats_k[player_i] = {"strategy_name_i":, 
                                      "q_plus_k_i":, "q_minus_k_i":,
                                      "P_t_plus_1_i":, "C_t_plus_1_i": ,
                                      "prod_i":, "cons_i":, "r_i":, 
                                      "S_t_plus_1_i": , 
                                      "Pi": , 
                                      "pp_t_i": , 
                                      "ben_i": , "cst_i":,
                                      "c0": , "b0":, 
                                      "pi_0_plus": , "pi_0_minus": , 
                                      "beta_sg_t_minus_1_minus":,
                                      "beta_sg_t_minus_1_plus":,
                                      "bg_i":, 'bg_min_i':, 'bg_max_i':, 
                                    }
    """
    bool_is_One_Ui_None = False 
    for player_i in dico_chosen_strats_k.keys():
        bg_min_i = dico_chosen_strats_k[player_i]["bg_min_i"]
        bg_max_i = dico_chosen_strats_k[player_i]["bg_max_i"]
        bg_i = dico_chosen_strats_k[player_i]["bg_i"]
        
        u_i = 1 - (bg_max_i - bg_i)/ (bg_max_i - bg_min_i) \
            if bg_max_i != bg_min_i \
            else None
            
        bool_is_One_Ui_None = False \
            if bg_max_i != bg_min_i and bool_is_One_Ui_None == False \
            else True
            
        dico_chosen_strats_k[player_i]["u_i"] = u_i
        
    return dico_chosen_strats_k, bool_is_One_Ui_None

#______________________________________________________________________________
#          compute u_i
#           fin
#______________________________________________________________________________ 