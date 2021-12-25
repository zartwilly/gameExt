#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Sun Nov 28 13:11:51 2021

@author: willy

goal: play game_ext1
"""

import os
import math
import json
import h5py
import time
import itertools as it
import numpy as np

import constances as csts
import generate_strategies as gene_strats
import fonctions_auxiliaires as fct_aux

from pathlib import Path

#______________________________________________________________________________
#                  compute prod, cons, r_i, remain_Si: debut
#______________________________________________________________________________
def compute_ProdConsRiSi(dico_tperiods_players, dico_PiStateiModei,
                         t_period, player_i, strategy_name_i):
    """
    STATE1_STRATS = ("CONS+", "CONS-")                                             # strategies possibles pour l'etat 1 de a_i
    STATE2_STRATS = ("DIS", "CONS-")                                               # strategies possibles pour l'etat 2 de a_i
    STATE3_STRATS = ("DIS", "PROD") 
    """
    Ci = dico_tperiods_players[t_period][player_i]["Ci"]
    Si = dico_tperiods_players[t_period][player_i]["Si"]
    Si_max = dico_tperiods_players[t_period][player_i]["Si_max"]
    Pi = dico_PiStateiModei["Pi"] 
    state_i = dico_PiStateiModei["state_i"]
    mode_i = dico_PiStateiModei["mode_i"] 
    strategy_name_i = strategy_name_i
    prod_i, cons_i, S_t_plus_1_i = None, None, None
    if state_i == csts.STATES[0]:
        prod_i = 0
        cons_i = Ci - (Pi + Si) if mode_i == csts.STATE1_STRATS[0] else Ci - Pi
        S_t_plus_1_i = 0 if mode_i == csts.STATE1_STRATS[0] else Si
    elif state_i == csts.STATES[1]:
        prod_i = 0
        cons_i = 0 if mode_i == csts.STATE2_STRATS[0] else Ci - Pi
        S_t_plus_1_i = Si - (Ci - Pi) if mode_i == csts.STATE2_STRATS[0] else Si
    elif state_i == csts.STATES[2]:
        cons_i = 0
        prod_i = fct_aux.diff_positive(Pi-Ci, Si_max-Si) \
            if mode_i == csts.STATE3_STRATS[0] else Pi - Ci
        S_t_plus_1_i = min(Si_max, Si+(Pi-Ci)) \
            if mode_i == csts.STATE3_STRATS[0] else Si
    else:
        print("ERROR state_i={state_i} UNKNOW")
        
    r_i = None
    if mode_i == csts.STATE1_STRATS[0]:
        r_i = 0
    elif mode_i == csts.STATE1_STRATS[1]:
        r_i = Si
    elif mode_i == csts.STATE2_STRATS[0] and state_i == csts.STATES[1]:
        r_i = 0
    elif mode_i == csts.STATE2_STRATS[0] and state_i == csts.STATES[2]:
        r_i = min(Si_max - Si, Pi - Ci)
    elif mode_i == csts.STATE3_STRATS[1]:
        r_i = 0
        
    return prod_i, cons_i, r_i, S_t_plus_1_i
    
#______________________________________________________________________________
#              compute prod, cons, r_i, remain_Si or Si_t_plus_1: fin
#______________________________________________________________________________

#______________________________________________________________________________
#       compute q_plus_i, q_minus_i, phi_EPO_plus, phi_EPO_minus : debut
#______________________________________________________________________________
def compute_qi_plus_minus(Pi, Ci, Si, Si_max):
    """
    compute q_plus et q_minus for player_i

    """
    q_minus_k_i = fct_aux.diff_positive(Ci, Pi) \
                    - fct_aux.diff_positive(Pi, Ci+Si_max-Si)
    q_plus_k_i = fct_aux.diff_positive(Pi, Ci) \
                    - fct_aux.diff_positive(Ci, Pi+Si)
                    
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
#                  compute gamma_i, Si_plus, minus: debut
#______________________________________________________________________________
def get_Pi_Ci_t_plus_1(dico_tperiods_players, t_periods, 
                       t_period, strategy_name_i, player_i):
    """
    """
    t = None
    if int(t_period.split('_')[1]) < t_periods-1:
        t = int(t_period.split('_')[1]) +1
    else:
        t = int(t_period.split('_')[1])
    P_t_plus_1_i = dico_tperiods_players\
                    [csts.PERIOD_ROOT+str(t)]\
                    [player_i]\
                    ["strategies"]\
                    [strategy_name_i]\
                    ["Pi"]
                    
    C_t_plus_1_i = dico_tperiods_players\
                    [csts.PERIOD_ROOT+str(t)]\
                    [player_i]\
                    ["Ci"]
    return P_t_plus_1_i, C_t_plus_1_i
                    
# def get_Si_plus_minus(state_i, Pi, dico_playeri_data):
#     """
#     compute maximun Si and minimun Si for one strategy
#     """
#     Si_plus, Si_minus = None, None
#     if state_i == csts.STATES[0]:
#         Si_minus = 0; Si_plus = dico_playeri_data["Si"]
#     elif state_i == csts.STATES[1]:
#         Si_minus = dico_playeri_data["Si"] - (dico_playeri_data["Ci"], Pi); 
#         Si_plus = dico_playeri_data["Si"]
#     elif state_i == csts.STATES[2]:
#         Si_minus = dico_playeri_data["Si"]
#         Si_plus = max(dico_playeri_data["Si_max"], 
#                       dico_playeri_data["Si"] + (Pi - dico_playeri_data["Ci"]))
#     else:
#         print(f"ERROR_state: {state_i} is UNKNOWN")
#     return Si_plus, Si_minus
    
# def compute_gamma(dico_players_tperiods, t_period, player_i, strategy_name_i,
#                   Pi_t_plus_1, Ci_t_plus_1,
#                   pi_0_plus, pi_0_minus, ppi_t):
#     """
#     compute gamma for player_i at the t_period
    
#     """
    
#     dico_strategy = dico_players_tperiods[t_period]\
#                                          [player_i]\
#                                          ["strategies"]\
#                                          [strategy_name_i]
#     Si_plus, Si_minus = None, None
#     Si_plus, Si_minus = get_Si_plus_minus(
#                             state_i = dico_strategy["state_i"], 
#                             Pi = dico_strategy["Pi"], 
#                             dico_playeri_data = dico_players_tperiods\
#                                                     [t_period][player_i])
#     print(f"Si_minus={Si_minus} <= Si_plus={Si_plus} OK") \
#         if Si_minus <= Si_plus \
#         else print(f"Si_minus={Si_minus} > Si_plus={Si_plus} NOK")
        
#     dico_strategy["Si_plus"] = Si_plus
#     dico_strategy["Si_minus"] = Si_minus
    
#     # compute pp_i
#     if ppi_t is None:
#         ppi_t = math.sqrt( min (1, 
#                                 fct_aux.diff_positive(
#                                     fct_aux.diff_positive(Ci_t_plus_1,
#                                                           Pi_t_plus_1),
#                                     Si_minus) / \
#                                     Si_plus - Si_minus
#                                 )
#                           )
#     # compute gamma
#     rd_draw = np.random.uniform(low=0.0, high=1.0, size=None)
#     rho_i_t = 1 if rd_draw < ppi_t else 0
#     X_gamV5 = None
#     X_gamV5 = pi_0_minus if dico_strategy["state_i"] == csts.STATES[0] \
#                             or dico_strategy["state_i"] == csts.STATES[1] \
#                          else pi_0_plus
#     gamma_i = rho_i_t * (X_gamV5 + 1)
#     dico_strategy["gamma_i"] = gamma_i
        
        
#     return 

#______________________________________________________________________________
#                  compute gamma_i, Si_plus, minus: fin
#______________________________________________________________________________

#______________________________________________________________________________
#                       choose strategy for all players: debut
#______________________________________________________________________________
def choose_strategy_4_all_players(dico_tperiods_players, t_period):
    """
    choose strategy for all players by following distribution 
    p_0s, p_1s, p_2s, p_3s, p_4s, p_5s
    
    NB:
        dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]}
    """
    q_plus_k, q_minus_k = 0, 0
    Out_sg_k, In_sg_k = 0, 0
    dico_chosen_strats_k = dict()     # key : player_X, values = dict()
    for player_i in dico_tperiods_players[t_period].keys():
        is_playing = dico_tperiods_players[t_period][player_i]["is_playing"][-1]
        if is_playing:
            keys_probas \
                = list(
                    filter(lambda x: x.startswith(csts.PROB_NAME), 
                           list(dico_tperiods_players[t_period][player_i].keys()) 
                          )
                    )
            vals_probas \
                = list( 
                    map(lambda x: dico_tperiods_players\
                                    [t_period][player_i][x][-1],
                        keys_probas) 
                    )
            
            strategy_name_i = np.random.choice(keys_probas, p=vals_probas)
            dico_PiStateiModei = dico_tperiods_players\
                                    [t_period][player_i]\
                                    ["strategies"][strategy_name_i]
                    
            prod_i, cons_i, r_i, S_t_plus_1_i \
                = compute_ProdConsRiSi(
                    dico_tperiods_players=dico_tperiods_players,
                    dico_PiStateiModei=dico_PiStateiModei,
                    t_period=t_period,
                    player_i=player_i, 
                    strategy_name_i=strategy_name_i
                    )
            Out_sg_k += cons_i
            In_sg_k += prod_i
            
            q_plus_k_i, q_minus_k_i \
                = compute_qi_plus_minus(
                    Pi=dico_tperiods_players[t_period][player_i]\
                                            ["strategies"][strategy_name_i]["Pi"],
                    Ci=dico_tperiods_players[t_period][player_i]["Ci"], 
                    Si=dico_tperiods_players[t_period][player_i]["S_t_plus_1_is"][-1],
                    Si_max=dico_tperiods_players[t_period][player_i]["Si_max"]
                    )
            q_plus_k += q_plus_k_i
            q_minus_k += q_minus_k_i
                
            P_t_plus_1_i, C_t_plus_1_i \
                = get_Pi_Ci_t_plus_1(
                    dico_tperiods_players=dico_tperiods_players, 
                    t_periods=t_periods, t_period=t_period, 
                    strategy_name_i=strategy_name_i, 
                    player_i=player_i
                    )
                
            dico_chosen_strats_k[player_i] \
                = {"strategy_name_i":strategy_name_i, 
                   "mode_i": dico_PiStateiModei['mode_i'], 
                   "state_i": dico_PiStateiModei['state_i'],
                   "Pi": dico_PiStateiModei['Pi'],
                   "q_plus_k_i":q_plus_k_i, "q_minus_k_i": q_minus_k_i,
                   "P_t_plus_1_i":P_t_plus_1_i, "C_t_plus_1_i": C_t_plus_1_i,
                   "prod_i":prod_i, "cons_i":cons_i, "r_i":r_i, 
                   "S_t_plus_1_i": S_t_plus_1_i
                   }
            
    return q_plus_k, q_minus_k, \
            Out_sg_k, In_sg_k, \
            dico_chosen_strats_k 
#______________________________________________________________________________
#                       choose strategy for all players: fin
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
    if t_period == 0:
        return beta_sg_0_minus, beta_sg_0_plus
    else:
        sum_T_diff_plus, sum_T_diff_minus = 0, 0
        sum_T_prod, sum_T_cons = 0, 0
        for t in range(0, int(t_period.split("_")[1] )):
            sum_N_cons_i = 0; sum_N_prod_i = 0
            for player_i, dico_values in dico_chosen_strats_t[t].items():
                prod_i = dico_chosen_strats_t[t][player_i]["prod_i"]
                cons_i = dico_chosen_strats_t[t][player_i]["cons_i"]
                
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
        Ii_m = fct_aux.diff_positive(op1=Pi, op2=Ci+(Si_max-Si))
        Ii_M = Pi - Ci
    else:
        print(f" ---> state={state_i} DON'T EXIST <---")
        
    return S_plus_i, S_minus_i, Ii_m, Ii_M, Oi_m, Oi_M

def compute_gamma_Siplusminus_ImM_OmM_bencstis(dico_chosen_strats_k, 
                                               dico_tperiods_players, 
                                               t_period, pp_t_i, 
                                               pi_0_minus, pi_0_plus, 
                                               beta_sg_t_minus_1_plus, 
                                               beta_sg_t_minus_1_minus,
                                               b0, c0):
    """
    compute gamma_i, Si_plis, Si_minus for all players
    compute also I_m, I_M, O_m, O_M 
    compute also ben_is, cst_is
    
    remember 
    dico_chosen_strats_k[player_i] = {"strategy_name_i":, 
                                      "q_plus_k_i":, "q_minus_k_i":,
                                      "P_t_plus_1_i":, "C_t_plus_1_i": ,
                                      "prod_i":, "cons_i":, "r_i":, 
                                      "S_t_plus_1_i": , 
                                      "Pi":
                                    }
    dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                    'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]}
    """
    I_m, I_M, O_m, O_M = 0, 0, 0, 0
    for player_i, dico_strats in dico_chosen_strats_k.items():
        strategy_name_i = dico_strats["strategy_name_i"]
        state_i = dico_tperiods_players\
                    [t_period][player_i]["strategies"][strategy_name_i]\
                    ['state_i']
        Pi = dico_tperiods_players\
                    [t_period][player_i]["strategies"][strategy_name_i]\
                    ['Pi']
        Ci = dico_tperiods_players[t_period][player_i]['Ci']
        Si = dico_tperiods_players[t_period][player_i]['Si']
        Si_max = dico_tperiods_players[t_period][player_i]['Si_max']
        
        S_plus_i, S_minus_i, Ii_m, Ii_M, Oi_m, Oi_M \
            = compute_Si_plusminus_ImM_OmM(state_i=state_i, 
                                            Si_max=Si_max, 
                                            Si=Si, Ci=Ci, Pi=Pi)
        # compute I_m, I_M, O_m, O_M
        I_m = Ii_m if Ii_m is None else I_m + Ii_m
        I_M = Ii_M if Ii_M is None else I_M + Ii_M
        O_m = Oi_m if Oi_m is None else O_m + Oi_m
        O_M = Oi_M if Oi_M is None else O_M + Oi_M
        
        # verification condition 
        assert S_minus_i <= S_minus_i
        
        
        # compute ppi_t
        P_t_plus_1_i = dico_chosen_strats_k[player_i]["P_t_plus_1_i"]
        C_t_plus_1_i = dico_chosen_strats_k[player_i]["C_t_plus_1_i"]
        if pp_t_i is None:
            pp_t_i = math.sqrt( min (1, 
                                    fct_aux.diff_positive(
                                        fct_aux.diff_positive(C_t_plus_1_i,
                                                               P_t_plus_1_i),
                                        S_minus_i) / \
                                        S_plus_i - S_minus_i
                                    )
                              )
                
        # compute gamma
        rd_draw = np.random.uniform(low=0.0, high=1.0, size=None)
        rho_i_t = 1 if rd_draw < pp_t_i else 0
        X_gamV5 = None
        X_gamV5 = pi_0_minus \
                    if state_i == csts.STATES[0] or state_i == csts.STATES[1] \
                    else pi_0_plus
        gamma_i = rho_i_t * (X_gamV5 + 1)
        dico_chosen_strats_k[player_i]["gamma_i"] = gamma_i
        dico_chosen_strats_k[player_i]["S_plus_i"] = S_plus_i
        dico_chosen_strats_k[player_i]["S_minus_i"] = S_minus_i
        dico_chosen_strats_k[player_i]["pp_t_i"] = pp_t_i
        
        # compute ben_i, cst_i
        prod_i = dico_chosen_strats_k[player_i]["prod_i"]
        cons_i = dico_chosen_strats_k[player_i]["cons_i"]
        r_i = dico_chosen_strats_k[player_i]["r_i"]
        ben_i = b0 * prod_i + gamma_i * r_i
        cst_i = c0 * cons_i
        dico_chosen_strats_k[player_i]["ben_i"] = ben_i
        dico_chosen_strats_k[player_i]["cst_i"] = cst_i
        dico_chosen_strats_k[player_i]["c0"] = c0
        dico_chosen_strats_k[player_i]["b0"] = b0
        dico_chosen_strats_k[player_i]["pi_0_minus"] = pi_0_minus
        dico_chosen_strats_k[player_i]["pi_0_plus"] = pi_0_plus
        dico_chosen_strats_k[player_i]["beta_sg_t_minus_1_plus"] = beta_sg_t_minus_1_plus
        dico_chosen_strats_k[player_i]["beta_sg_t_minus_1_minus"] = beta_sg_t_minus_1_minus
        
    return dico_chosen_strats_k, I_m, I_M, O_m, O_M
        
#______________________________________________________________________________
#          compute 
#           S_minus_i, S_plus_i, I_m, I_M, O_m, O_M 
#           ben_i, cst_i 
#           fin
#______________________________________________________________________________    

#______________________________________________________________________________
#          compute bg
#           debut
#______________________________________________________________________________  
def compute_bg(dico_chosen_strats_k, dico_tperiods_players, t_period, c0_M):
    """
    determine bg, bg_i_min, bg_i_max
    dico_chosen_strats_k[player_i] = {"strategy_name_i":, 
                                      "q_plus_k_i":, "q_minus_k_i":,
                                      "P_t_plus_1_i":, "C_t_plus_1_i": ,
                                      "prod_i":, "cons_i":, "r_i":, 
                                      "S_t_plus_1_i": , 
                                      "S_minus_i": , "S_plus_i": 
                                      "Pi": , 
                                      "pp_t_i": , 
                                      "ben_i": , "cst_i":,
                                      "c0": , "b0:"
                                      "pi_0_plus": , "pi_0_minus": , 
                                      "beta_sg_t_minus_1_minus":,
                                      "beta_sg_t_minus_1_plus":
                                    }
    dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]}
    """
    for player_i, dico_strats_k in dico_chosen_strats_k.items():
        strategy_name_i = dico_strats_k["strategy_name_i"]
        Pi = dico_tperiods_players\
                    [t_period][player_i]["strategies"][strategy_name_i]\
                    ['Pi']
        Ci = dico_tperiods_players[t_period][player_i]['Ci']
        
        bg_i = dico_tperiods_players[t_period][player_i]['ben_i'] \
                - dico_tperiods_players[t_period][player_i]['cst_i'] \
                + c0_M * fct_aux.diff_positive(op1=Ci, op2=Pi)
                
        dico_chosen_strats_k[player_i]["bg_i"] = bg_i
        
        bg_min_i ,bg_max_i = None, None
        
        bg_min_i = min(dico_tperiods_players\
                           [t_period][player_i]\
                               .get(csts.BG_ALL_VALUES))
        bg_max_i = max(dico_tperiods_players\
                           [t_period][player_i]\
                               .get(csts.BG_ALL_VALUES))
        
        if bg_min_i is None or bg_min_i > bg_i:
            bg_min_i = bg_i
            
        if bg_max_i is None or bg_max_i < bg_i:
            bg_max_i = bg_i
        
        dico_chosen_strats_k[player_i]["bg_i_min"] = bg_min_i
        dico_chosen_strats_k[player_i]["bg_i_max"] = bg_max_i
        
    return dico_chosen_strats_k

#______________________________________________________________________________
#          compute bg
#           fin
#______________________________________________________________________________  

#______________________________________________________________________________
#          compute u_i
#           debut
#______________________________________________________________________________  
def compute_utility_fonction_ui(dico_chosen_strats_k, 
                                dico_tperiods_players, 
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
                                      "bg_i":, 'bg_i_min':, 'bg_i_max':, 
                                    }
    dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]}
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

#______________________________________________________________________________
#          update p_Xs
#           debut
#______________________________________________________________________________
def update_probabilities_pXs(dico_chosen_strats_k, 
                             dico_tperiods_players, 
                             t_period, 
                             learning_rate):
    """
    update probability for choosen strategy
    
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
                                      "bg_i":, 'bg_min_i':, 'bg_max_i':,
                                      "beta_sg_t_minus_1_plus": ,
                                      "beta_sg_t_minus_1_minus": ,
                                    }
    dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]}
    
    """
    
    for player_i in dico_chosen_strats_k.keys():
        u_i = dico_chosen_strats_k[player_i]["u_i"]
        strategy_name_i = dico_chosen_strats_k[player_i]["strategy_name_i"]
        p_X = dico_tperiods_players[t_period][player_i][strategy_name_i][-1]
        
        increasing_val = learning_rate * u_i * (1 - p_X)
        p_X_new = p_X + increasing_val
        dico_chosen_strats_k[player_i][strategy_name_i] = p_X_new
        
        is_playing = True if p_X_new < csts.MAX_STRATEGY_PROBA else False
        dico_chosen_strats_k[player_i][is_playing] = is_playing
        
        # reduce remain strategies to  p_X_new
        # strats = {1,2,3,4,5,6}, X in starts
        # if X=1 and p_X = 0.2 and p_X_new = 0.3 then
        #  increasing_val = 0.1
        #  Y = strats - {X}
        #  p_Y_new = p_Y -  increasing_val / len(Y)
        strategy_names \
                = set(
                    filter(lambda x: x.startswith(csts.PROB_NAME), 
                           list(dico_tperiods_players[t_period][player_i].keys()) 
                          )
                    )
        '''
        Y = strategy_names - set({strategy_name_i})
        for y in Y:
            # y is the name of not choosen strategy
            p_y = dico_players_tperiods[t_period][player_i][y][-1]
            p_y_new = p_y - increasing_val / len(Y)
            dico_players_tperiods\
                [t_period]\
                [player_i]\
                [y].append(p_y_new)
        '''
           
        Y = strategy_names - set({strategy_name_i})
        no_choose_strategy_names = strategy_names - set({strategy_name_i})
        for no_choose_strategy_name_i in no_choose_strategy_names:
            p_y = dico_tperiods_players\
                    [t_period][player_i]\
                    [no_choose_strategy_name_i][-1]
            p_y_new = p_y - increasing_val / len(Y)
            dico_chosen_strats_k[player_i][no_choose_strategy_name_i] = p_y_new
                
    return dico_chosen_strats_k

#______________________________________________________________________________
#          update p_Xs
#           fin
#______________________________________________________________________________  

#______________________________________________________________________________
#          update learning variables 
#           debut
#______________________________________________________________________________  
def update_learning_variables(dico_chosen_strats_k, 
                              dico_tperiods_players, 
                              t_period, is_repeated_kstep):
    """
    update learning variables ie add LEARNING_VARIABLES to dico_players_tperiods 
    
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
                                      "beta_sg_t_minus_1_minus": ,
                                      "beta_sg_t_minus_1_plus": ,
                                      "bg_i":, 'bg_min_i':, 'bg_max_i':, 
                                    }
    dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]}
    """
    for player_i in dico_chosen_strats_k.keys():
        mode_i = dico_chosen_strats_k[player_i]["mode_i"]
        state_i = dico_chosen_strats_k[player_i]["state_i"]
        Pi = dico_chosen_strats_k[player_i]["Pi"]
        S_t_plus_1_i = dico_chosen_strats_k[player_i]["S_t_plus_1_i"]
        k_step = dico_chosen_strats_k[player_i]["k_step"]
        S_minus_i = dico_chosen_strats_k[player_i]["S_minus_i"]
        S_plus_i = dico_chosen_strats_k[player_i]["S_plus_i"]
        pp_t_i = dico_chosen_strats_k[player_i]["pp_t_i"]
        gamma_i = dico_chosen_strats_k[player_i]["gamma_i"]
        q_minus_k_i = dico_chosen_strats_k[player_i]["q_minus_k_i"]
        q_plus_k_i = dico_chosen_strats_k[player_i]["q_plus_k_i"]
        prod_i = dico_chosen_strats_k[player_i]["prod_i"]
        cons_i = dico_chosen_strats_k[player_i]["cons_i"]
        r_i = dico_chosen_strats_k[player_i]["r_i"]
        ben_i = dico_chosen_strats_k[player_i]["ben_i"]
        cst_i = dico_chosen_strats_k[player_i]["cst_i"]
        V_i = ben_i - cst_i
        bg_i = dico_chosen_strats_k[player_i]["bg_i"]
        bg_min_i = dico_chosen_strats_k[player_i]["bg_min_i"]
        bg_max_i = dico_chosen_strats_k[player_i]["bg_max_i"]
        u_i = dico_chosen_strats_k[player_i]["u_i"]
        
        is_playing = dico_chosen_strats_k[player_i]["is_playing"] 
        
        c0 = dico_chosen_strats_k[player_i]["c0"]
        b0 = dico_chosen_strats_k[player_i]["b0"]
        pi_0_plus = dico_chosen_strats_k[player_i]["pi_0_plus"]
        pi_0_minus = dico_chosen_strats_k[player_i]["pi_0_minus"]
        beta_sg_t_minus_1_minus = dico_chosen_strats_k[player_i]["beta_sg_t_minus_1_minus"]
        beta_sg_t_minus_1_plus = dico_chosen_strats_k[player_i]["beta_sg_t_minus_1_plus"]
        
        vars_2_add = [("mode_is", mode_i), ("state_i", state_i),
                      ("P_is", Pi), ("S_t_plus_1_is", S_t_plus_1_i), 
                      ("k_steps", k_step),
                      ("S_minus_is", S_minus_i), ("S_plus_is", S_plus_i), 
                      ("pp_t_is", pp_t_i), ("gamma_is", gamma_i), 
                      ("q_minus_is", q_minus_k_i), ("q_plus_is", q_plus_k_i), 
                      ("prod_is", prod_i), ("cons_is", cons_i), ("r_is", r_i), 
                      ("ben_is", ben_i), ("cst_is", cst_i), ("V_is", V_i),
                      ("u_is", u_i), 
                      ("is_playing", is_playing),
                      ("bg_is", bg_i), ("bg_min_is", bg_min_i), 
                      ("bg_max_is", bg_max_i), (csts.BG_ALL_VALUES, bg_i), 
                      ("c0_s", c0), ("b0_s", b0), 
                      ("pi_0_plus_s", pi_0_plus), ("pi_0_minus_s", pi_0_minus),
                      ("beta_sg_t_minus_1_minus_is", beta_sg_t_minus_1_minus),
                      ("beta_sg_t_minus_1_plus_is", beta_sg_t_minus_1_plus)
                      ]
                      
        if is_repeated_kstep:
            var_i = csts.BG_ALL_VALUES
            val_i = bg_i
            dico_tperiods_players[t_period][player_i][var_i].add(val_i)
        else:
            for var_i, val_i in vars_2_add:
                if var_i != csts.BG_ALL_VALUES:
                    dico_tperiods_players[t_period][player_i][var_i].append(val_i)
                else:
                    dico_tperiods_players[t_period][player_i][var_i].add(val_i)
            
    return dico_tperiods_players
    
#______________________________________________________________________________
#          update learning variables 
#           fin
#______________________________________________________________________________  
    
#______________________________________________________________________________
#          execute one learning step
#           debut
#______________________________________________________________________________
def execute_one_learning_step_4_one_period(dico_tperiods_players, 
                                           dico_chosen_strats_t, 
                                           t_period, k_step,
                                           cpt_repeated_kstep, args ):
    """
    """
    q_plus_k, q_minus_k, \
    Out_sg_k, In_sg_k, \
    dico_chosen_strats_k \
        = choose_strategy_4_all_players(
            dico_tperiods_players=dico_tperiods_players, 
            t_period=t_period
            )
    
    # compute q_t_plus, q_t_minus
    q_plus_k = 0 if q_plus_k < 0 else q_plus_k
    q_minus_k = 0 if q_minus_k < 0 else q_minus_k
    phi_EPO_plus = compute_phi_EPO_plus(q_plus_k=q_plus_k, 
                                        quantity=args["quantity_a"], 
                                        a=args["a"])
    phi_EPO_minus = compute_phi_EPO_minus(q_minus_k=q_minus_k, 
                                          quantity=args["quantity_b"], 
                                          b=args["b"])
    pi_EPO_plus = round(phi_EPO_plus/q_plus_k, csts.ARRONDI) \
                    if q_plus_k > 0 \
                    else 0
    pi_EPO_minus = round(phi_EPO_minus/q_minus_k, csts.ARRONDI) \
                    if q_minus_k > 0 \
                    else 0
    
    beta_sg_t_minus_1_plus, beta_sg_t_minus_1_minus \
        = compute_beta_sg(t_period=t_period, 
                          dico_chosen_strats_t=dico_chosen_strats_t,
                          beta_sg_0_minus=pi_EPO_minus-1, 
                          beta_sg_0_plus=pi_EPO_plus-1,
                          quantity_a=args["quantity_a"],  
                          a=args["a"],  
                          quantity_b=args["quantity_b"],        
                          b=args["b"]
                          )

    pi_0_plus, pi_0_minus \
        = compute_pi_0(t_period=t_period, 
                       beta_sg_t_minus_1_plus=beta_sg_t_minus_1_plus, 
                       beta_sg_t_minus_1_minus=beta_sg_t_minus_1_minus,
                       pi_EPO_plus=pi_EPO_plus, 
                       pi_EPO_minus=pi_EPO_minus
                       )

    b0, c0 = compute_b0_c0(pi_0_plus=pi_0_plus, 
                       pi_0_minus=pi_0_minus, 
                       Out_sg_k=Out_sg_k, In_sg_k=In_sg_k, 
                       quantity_a=args["quantity_a"], a=args["a"], 
                       quantity_b=args["quantity_b"], b=args["b"])
    
    # compute gamma_i, Si_plus, Si_minus, pp_i for all players
    # compute also I_m, I_M, O_m, O_M 
    # compute also ben_is, cst_is
    dico_chosen_strats_k, I_m, I_M, O_m, O_M \
        = compute_gamma_Siplusminus_ImM_OmM_bencstis(
            dico_chosen_strats_k=dico_chosen_strats_k, 
            dico_tperiods_players=dico_tperiods_players, 
            t_period=t_period, 
            ppi_t=args["pp_t_i"], 
            pi_0_minus=pi_0_minus, pi_0_plus=pi_0_plus, 
            beta_sg_t_minus_1_plus=beta_sg_t_minus_1_plus, 
            beta_sg_t_minus_1_minus=beta_sg_t_minus_1_minus,
            b0=b0, c0=c0
            )
    
    assert In_sg_k >= I_m and In_sg_k <= I_m
    assert Out_sg_k >= O_m and Out_sg_k <= O_m

    # compute c0_M
    c0_M = min(
            round(((O_M - I_m)*pi_EPO_minus + I_M*pi_0_minus)/O_m, 
                      csts.ARRONDI), 
            pi_0_minus)
    assert c0 <= c0_M
    
    # compute bg for all players
    # add fields bg_min_i and bg_max_i in dico_chosen_strats_k
    dico_chosen_strats_k = compute_bg(
                            dico_chosen_strats_k=dico_chosen_strats_k, 
                            dico_tperiods_players=dico_tperiods_players, 
                            t_period=t_period)
    
    # compute utility fonction
    dico_chosen_strats_k, bool_is_One_Ui_None \
        = compute_utility_fonction_ui(
            dico_chosen_strats_k=dico_chosen_strats_k, 
            dico_tperiods_players=dico_tperiods_players, 
            t_period=t_period
            )

    # update probabilities in case of all u_i are different to None
    if bool_is_One_Ui_None == False \
        or cpt_repeated_kstep >= csts.MAX_REPEATED_KSTEP-1:
        dico_chosen_strats_k \
            = update_probabilities_pXs(
                dico_chosen_strats_k=dico_chosen_strats_k, 
                dico_tperiods_players=dico_tperiods_players, 
                t_period=t_period, 
                learning_rate=args['learning_rate'])
        
        # update startegy keys and other variables to dico_tperiods_players
        is_repeated_kstep = False
        dico_tperiods_players \
            = update_learning_variables(
                dico_chosen_strats_k=dico_chosen_strats_k, 
                dico_tperiods_players=dico_tperiods_players, 
                t_period=t_period, 
                is_repeated_kstep=is_repeated_kstep)
            
        k_step += 1
        cpt_repeated_kstep = csts.MAX_REPEATED_KSTEP
        
    else:
        """ else means bool_is_One_Ui_None == True \
                        and cpt_repeated_kstep < csts.MAX_REPEATED_KSTEP
        """
        
        # update startegy keys and other variables to dico_tperiods_players
        is_repeated_kstep = True
        dico_tperiods_players \
            = update_learning_variables(
                dico_chosen_strats_k=dico_chosen_strats_k, 
                dico_tperiods_players=dico_tperiods_players, 
                t_period=t_period, 
                is_repeated_kstep=is_repeated_kstep)
        
        k_step = k_step
        cpt_repeated_kstep += 1
        
    dico_chosen_strats_t[t_period] = dico_chosen_strats_k
    
    return dico_tperiods_players, dico_chosen_strats_t, \
            k_step, cpt_repeated_kstep
    
#______________________________________________________________________________
#          execute one learning step
#           fin
#______________________________________________________________________________

#______________________________________________________________________________
#          send residual stock to next period
#           debut
#______________________________________________________________________________
def send_residual_stock_2_next_period(dico_tperiods_players, t_period, args):
    """
    pass the remaining stock at t to the new stock at t+1 
    """
    t = int(t_period.split('_')[-1])
    if t < args["t_periods"]:
        next_period = csts.PERIOD_ROOT + str(t+1)
    else:
        next_period = t_period
        
    liste_players = list()
    for player_i in dico_tperiods_players[t_period].keys():
        dico_tperiods_players[next_period]\
                             [player_i]\
                             ["Si"] \
            = dico_tperiods_players\
                [t_period]\
                [player_i]\
                ["S_t_plus_1_is"][-1]
                
        liste_players.append(player_i)
                    
    return dico_tperiods_players, list(set(liste_players))
        
#______________________________________________________________________________
#          send residual stock to next period
#           fin
#______________________________________________________________________________

#______________________________________________________________________________
#                   compute IB, IC and Perft
#                           debut
#______________________________________________________________________________
def compute_IB_IC_Perft(dico_tperiods_players, dico_chosen_strats_t, 
                        liste_players):
    """
    compute IB, IC and Perft
    
    NB: dico_chosen_strats_t = {
        't_0': {player_0: dico_pl0, ..., player_M: dico_plM}, 
        ...,
        't_T': {player_0: dico_pl0, ..., player_M: dico_plM}, 
        }
     with dico_pl0 = {"strategy_name_i":, 
                        "q_plus_k_i":, "q_minus_k_i":,
                        "P_t_plus_1_i":, "C_t_plus_1_i": ,
                        "prod_i":, "cons_i":, "r_i":, 
                        "S_t_plus_1_i": , 
                        "Pi": , 
                        "pp_t_i": , 
                        "ben_i": , "cst_i":,
                        "c0": , "b0":, 
                        "pi_0_plus": , "pi_0_minus": , 
                        "beta_sg_t_minus_1_minus": ,
                        "beta_sg_t_minus_1_plus": ,
                        "bg_i":, 'bg_min_i':, 'bg_max_i':, 
                        }
     
     IB = sum_{t=1}^{T} b0_t*prod_i^t, IC = sum_{t=1}^{T} c0_t*cons_i^t
    """
    dico_IB_IC = dict()
    for player_i in liste_players:
        dico_IB_IC[player_i] = {"IB_i_t":[], "IC_i_t":[], 
                                "CONS_i_t":[], "PROD_i_t":[]}
        
    dico_Perft = dict()
    for t_period in dico_chosen_strats_t.keys():
        Perf_t = 0
        for player_i in dico_chosen_strats_t[t_period].keys():
            c0_t = dico_chosen_strats_t[t_period][player_i]["c0"]
            cons_i_t = dico_chosen_strats_t[t_period][player_i]["cons_i"]
            b0_t = dico_chosen_strats_t[t_period][player_i]["b0"]
            prod_i_t = dico_chosen_strats_t[t_period][player_i]["prod_i"]
            ben_i_t = dico_chosen_strats_t[t_period][player_i]["ben_i"]
            cst_i_t = dico_chosen_strats_t[t_period][player_i]["cst_i"]
            beta_sg_t_minus = dico_chosen_strats_t[t_period][player_i]\
                                        ["beta_sg_t_minus_1_minus"]
            beta_sg_t_plus = dico_chosen_strats_t[t_period][player_i]\
                                        ["beta_sg_t_minus_1_plus"]
            
            # beta_sg_t_minus = dico_chosen_strats_t[t_period][player_i]["beta_sg_t_minus"]
            # beta_sg_t_plus = dico_chosen_strats_t[t_period][player_i]["beta_sg_t_plus"]
            
            IB_i_t = b0_t * prod_i_t
            IC_i_t = c0_t * cons_i_t
            V_i_t = ben_i_t - cst_i_t
            Perf_t += V_i_t
            
            dico_IB_IC[player_i]["IB_i_t"].append(IB_i_t)
            dico_IB_IC[player_i]["IC_i_t"].append(IC_i_t)
            dico_IB_IC[player_i]["CONS_i_t"].append(cons_i_t)
            dico_IB_IC[player_i]["PROD_i_t"].append(prod_i_t)
            dico_IB_IC[player_i]["beta_sg_t_minus"] = beta_sg_t_minus
            dico_IB_IC[player_i]["beta_sg_t_plus"] = beta_sg_t_plus
            
        dico_Perft[t_period] = Perf_t
            
    APerf = sum(dico_Perft.values()) / len(list(dico_chosen_strats_t.keys()))
        
    VR = 0
    for player_i in dico_IB_IC.keys():
        dico_IB_IC[player_i]["IB_i"] \
            = sum(dico_IB_IC[player_i]["IB_i_t"])
        dico_IB_IC[player_i]["IC_i"] \
            = sum(dico_IB_IC[player_i]["IC_i_t"])
            
        VR += dico_IB_IC[player_i]["IB_i"] - dico_IB_IC[player_i]["IC_i"]
        
        dico_IB_IC[player_i]["CONS_i"] \
            = sum(dico_IB_IC[player_i]["CONS_i_t"])
        dico_IB_IC[player_i]["PROD_i"] \
            = sum(dico_IB_IC[player_i]["PROD_i_t"])
        dico_IB_IC[player_i]["BB_i"] \
            = dico_IB_IC[player_i]["beta_sg_t_plus"] \
                * sum(dico_IB_IC[player_i]["PROD_i_t"])
        dico_IB_IC[player_i]["CC_i"] \
            = dico_IB_IC[player_i]["beta_sg_t_minus"] \
                * sum(dico_IB_IC[player_i]["CONS_i_t"])
            
    dico_VR_APerf = {"VR":VR, "APerf":APerf}
    
    return dico_IB_IC, dico_Perft, dico_VR_APerf

#______________________________________________________________________________
#                   compute IB, IC and Perft
#                           fin
#______________________________________________________________________________

#______________________________________________________________________________
#                   save learning variables
#                           debut
#______________________________________________________________________________
def save_learning_variables(dico_tperiods_players, dico_chosen_strats_t,
                            dico_IB_IC_Perft, dico_Perft, dico_VR_APerf, 
                            path_2_save):
    """
    save learning variables to json format or csv


    """
    Path(path_2_save).mkdir(parents=True, exist_ok=True)
    # save to json dico_tperiods_players
    with open(os.path.join(csts.PATH_ROOT,'data_tperiods_players.json'), 'w') \
        as fp:
        json.dump(dico_tperiods_players, fp)
        
    # save to json dico_chosen_strats_t
    with open(os.path.join(csts.PATH_ROOT,'data_chosen_strategy_4_players.json'), 
              'w') \
        as fp:
        json.dump(dico_chosen_strats_t, fp)
        
    # save to json dico_IB_IC_Perft
    with open(os.path.join(csts.PATH_ROOT,'data_players_IB_IC_Perft.json'), 'w') \
        as fp:
        json.dump(dico_IB_IC_Perft, fp)
    
    # save to json dico_Perft
    with open(os.path.join(csts.PATH_ROOT,'data_tperiods_Perft.json'), 'w') \
        as fp:
        json.dump(dico_Perft, fp)
        
    # save to json dico_VR_APerf
    with open(os.path.join(csts.PATH_ROOT,'data_VR_Perft.json'), 'w') \
        as fp:
        json.dump(dico_VR_APerf, fp)
        
    
    # save all dicos to hdf5
    liste_tup_jsonDico \
        = [("data_tperiods_players", 
            json.dumps(dico_tperiods_players, indent=4)), 
           ("data_chosen_strategy_4_players", 
            json.dumps(dico_chosen_strats_t, indent=4)), 
           ("data_players_IB_IC_Perft",
            json.dumps(dico_IB_IC_Perft, indent=4)),
           ("data_tperiods_Perft", json.dumps(dico_Perft, indent=4)), 
           ("data_VR_Perft", json.dumps(dico_VR_APerf, indent=4))
          ]
    with h5py.File(
            os.path.join(csts.DATA_ROOT, 
                         f'data_{n_instances}_instances.hdf5'), 'w') as f:
        for (id_json, json_data) in liste_tup_jsonDico:
            f.create_dataset(id_json, data=json_data)
            
    print("----> Saving all dicos with SUCCES......")
    
    
#______________________________________________________________________________
#                   save learning variables
#                           fin
#______________________________________________________________________________


#______________________________________________________________________________
#                       game extension 1: debut
#______________________________________________________________________________
def game_ext1(dico_tperiods_players, args):
    
    """
    * add probabilities for all strategies : 
            already insert in the generate_strategies function
    *  for each t 
        * for each k_step
            * for each player
                * choose one strategy from probabilities list
                * compute cons_i, prod_i, ...
            * compute in_sg, out_sg, beta_sg_{+,-}
            * 
    NB : dico_tperiods_players is a variable associated to one instance
    """
    dico_chosen_strats_t = dict(); liste_players = list()
    for t_period in dico_tperiods_players.keys():
        
        k_step = 0; cpt_repeated_kstep = 0
        while k_step < args["k_steps"] \
                and cpt_repeated_kstep < csts.MAX_REPEATED_KSTEP:
                    
            dico_tperiods_players, dico_chosen_strats_t, \
            k_step, cpt_repeated_kstep \
                = execute_one_learning_step_4_one_period(
                    dico_tperiods_players=dico_tperiods_players, 
                    dico_chosen_strats_t=dico_chosen_strats_t,
                    t_period=t_period, 
                    k_step=k_step,
                    cpt_repeated_kstep=cpt_repeated_kstep,
                    args=args
                    )
                
        # pass the remaining stock at t to the new stock at t+1 
        dico_tperiods_players, liste_players \
            = send_residual_stock_2_next_period(
                dico_tperiods_players=dico_tperiods_players, 
                t_period=t_period, 
                args=args)
        # TODO : TO THINK- save each kstep or each t_period
        
    # compute IB, IC, EB_i, EC_i
    dico_IB_IC_Perft, dico_Perft, dico_VR_APerf \
        = compute_IB_IC_Perft(dico_tperiods_players=dico_tperiods_players, 
                              dico_chosen_strats_t=dico_chosen_strats_t, 
                              liste_players=liste_players)
        
    # save learning process to all periods
    save_learning_variables(dico_tperiods_players=dico_tperiods_players, 
                            dico_chosen_strats_t=dico_chosen_strats_t,
                            dico_IB_IC_Perft=dico_IB_IC_Perft, 
                            dico_Perft=dico_Perft, 
                            dico_VR_APerf=dico_VR_APerf,
                            path_2_save=args["path_2_save"])
    
    
#______________________________________________________________________________
#                       game extension 1: fin
#______________________________________________________________________________


#______________________________________________________________________________
#                       test : debut
#______________________________________________________________________________
def test_gameExt1(args):
    
    dico_insti_tperiods_players \
        = gene_strats.read_json_filename(json_filename=args['json_filename'])
        
    dico_insti_tperiods_players \
        = gene_strats.add_fields_players(dico_data=dico_insti_tperiods_players)
        
    # test generations strategies
    dico_insti_tperiods_players \
        = gene_strats.generation_strategies_tperiods(
            dico_data=dico_insti_tperiods_players, 
            t_periods=args['t_periods'])
        
    game_ext1(dico_tperiods_players=dico_insti_tperiods_players, 
              args=args)
    pass
#______________________________________________________________________________
#                       test : fin
#______________________________________________________________________________


if __name__ == "__main__":
    t_periods = 5
    k_steps = 1000
    n_instances = 10
    n_instance = np.random.choice(range(10))
    
    json_filename = f"data_instance_{n_instance}_tperiods_{t_periods}.json"
    
    quantity_a = 10; a = 1;
    quantity_b = 30; b = 1;
    pp_t_i = None
    learning_rate = 0.1
    path_2_save = os.path.join(".", csts.PATH_ROOT, 
                               csts.INSTANCE_ROOT+str(n_instance))
    
    args = {'json_filename': json_filename, 
            't_periods': t_periods,
            'n_instance': n_instance, 
            "k_steps": k_steps,
            "quantity_a": quantity_a, "a": a,
            "quantity_b": quantity_b, "b": b,
            "pp_t_i": pp_t_i, 
            "learning_rate": learning_rate, 
            "path_2_save": path_2_save
            }
    
    test_gameExt1(args)
    
    
    