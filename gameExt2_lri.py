#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 23:22:17 2022

@author: willy

production_game (f2) with initial_game (f1)
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
#                       choose f12-strategy for all players: debut
#______________________________________________________________________________
def choose_f12_strategy_4_all_players(dico_insti_tperiods_players, t_period, r_step):
    """
    choose f12_strategy for all players by following distribution 
    r_0_s, r_1_s, r_2s
    
    NB:
        dico_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and
        r_step \in {0, ..., 100} 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES],
                            r_0_s = [r_0_1, r_0_2, ..., r_0_100]
                         r_1_s = [r_1_1, r_1_2, ..., r_1_100]
                         r_2_s = [r_2_1, r_2_2, ..., r_2_100]
                         f12_strategies_names = ["r_0_s", "r_1_s", "r_2_s"]
                         f12_strategies_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]   
                            }
    """
    dico_chosen_strats_t_r = dict()     # key : player_X, values = dict()
    for player_i in dico_insti_tperiods_players[t_period].keys():
        is_playing = dico_insti_tperiods_players[t_period][player_i]["is_playing"][-1]
        if is_playing:
            keys_probas \
                = list(
                    filter(lambda x: x.startswith(csts.F21_NAME) \
                                       and x.endswith("_s"), 
                           list(dico_insti_tperiods_players[t_period][player_i].keys()) 
                          )
                    )
            vals_probas \
                = list( 
                    map(lambda x: dico_insti_tperiods_players\
                                    [t_period][player_i][x]["prob"],
                        keys_probas) 
                    )
            #print("proba = {}".format( list(zip(keys_probas,vals_probas)) ))
            
            f12_strategy_name_i = np.random.choice(keys_probas, p=vals_probas)
            dico_insti_tperiods_players\
                [t_period][player_i]\
                    ["f12_strategies_names_is"].append(f12_strategy_name_i)
                    
            # save variables values of player_i to dictionnary
            dico_prob_Pi_p_01s = dict()
            dico_prob_Pi_p_01s["Pi"] = dico_insti_tperiods_players\
                                        [t_period][player_i]\
                                        [f12_strategy_name_i]["di"]
            dico_prob_Pi_p_01s["p_0s"] = dico_insti_tperiods_players\
                                            [t_period][player_i]\
                                            [f12_strategy_name_i]["p_0s"][-1]
            dico_prob_Pi_p_01s["p_1s"] = dico_insti_tperiods_players\
                                            [t_period][player_i]\
                                            [f12_strategy_name_i]["p_1s"][-1]
            dico_prob_Pi_p_01s["prob"] = dico_insti_tperiods_players\
                                            [t_period][player_i]\
                                            [f12_strategy_name_i]["prob"]
            dico_prob_Pi_p_01s["f12_strategy_name_i"] = f12_strategy_name_i
            
                
            dico_chosen_strats_t_r[player_i] = dico_prob_Pi_p_01s
            
    return dico_chosen_strats_t_r
#______________________________________________________________________________
#                       choose f12-strategy for all players: fin
#______________________________________________________________________________  

#______________________________________________________________________________
#                       fonctions auxiliaires de 
#           f12_f22_execute_one_learning_step_4_one_period: debut
#______________________________________________________________________________  
def compute_StateModeProdConsRiSi(dico_insti_tperiods_players,
                                  t_period, player_i, Pi, probas_p_0s1s):
    """
    """
    Si = dico_insti_tperiods_players[t_period][player_i]["Si"]
    Si_max = dico_insti_tperiods_players[t_period][player_i]["Si_max"]
    Ci = dico_insti_tperiods_players[t_period][player_i]["Ci"]
    state_i, mode_i = None, None
    prod_i, cons_i, S_t_plus_1_i = 0, 0, 0
    if Pi + Si <= Ci:
        state_i = csts.STATES[0]                                                # Deficit
        mode_i = np.random.choice(csts.STATE1_STRATS, p=probas_p_0s1s)
        prod_i = 0
        cons_i = Ci - (Pi + Si) if mode_i == csts.STATE1_STRATS[0] else Ci - Pi
        S_t_plus_1_i = 0 if mode_i == csts.STATE1_STRATS[0] else Si
    elif Pi + Si > Ci and Pi <= Ci:
        state_i = csts.STATES[1]                                                # Self
        mode_i = np.random.choice(csts.STATE2_STRATS, p=probas_p_0s1s)
        prod_i = 0
        cons_i = 0 if mode_i == csts.STATE2_STRATS[0] else Ci - Pi
        #S_t_plus_1_i = Si - (Ci - Pi) if mode_i == csts.STATE2_STRATS[0] else Si
        S_t_plus_1_i = fct_aux.diff_positive(op1=Si, op2=Ci-Pi) \
            if mode_i == csts.STATE2_STRATS[0] else Si
    elif Pi >= Ci:
        state_i = csts.STATES[2]                                                # Surplus
        mode_i = np.random.choice(csts.STATE3_STRATS, p=probas_p_0s1s)
        cons_i = 0
        prod_i = fct_aux.diff_positive(Pi-Ci, Si_max-Si) \
            if mode_i == csts.STATE3_STRATS[0] else Pi - Ci
        S_t_plus_1_i = min(Si_max, Si+(Pi-Ci)) \
            if mode_i == csts.STATE3_STRATS[0] else Si
    
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
        
    return state_i, mode_i, prod_i, cons_i, r_i, S_t_plus_1_i

def get_Pi_Ci_t_plus_1(dico_insti_tperiods_players, t_periods, 
                       t_period, strategy_name_i, player_i):
    """
    """
    t = None
    if int(t_period.split('_')[1]) < t_periods-1:
        t = int(t_period.split('_')[1]) +1
    else:
        t = int(t_period.split('_')[1])
        
    # TODO IMPORTANT: problem to get P_t_plus_1_i, a voir avec Dominique
    P_t_plus_1_i = dico_insti_tperiods_players\
                    [csts.PERIOD_ROOT+str(t)]\
                    [player_i]\
                    ["strategies"]\
                    [strategy_name_i]\
                    ["Pi"]
                    
    C_t_plus_1_i = dico_insti_tperiods_players\
                    [csts.PERIOD_ROOT+str(t)]\
                    [player_i]\
                    ["Ci"]
    return P_t_plus_1_i, C_t_plus_1_i

#______________________________________________________________________________
#                       fonctions auxiliaires de 
#           f12_f22_execute_one_learning_step_4_one_period: fin
#______________________________________________________________________________  

#______________________________________________________________________________
#                       choose strategy for all players: debut
#______________________________________________________________________________
def choose_strategy_4_all_players(dico_insti_tperiods_players, 
                                    dico_chosen_strats_t_r, t_period):
    """
    """
    Out_sg_k, In_sg_k = 0, 0
    q_plus_k, q_minus_k = 0, 0
    for player_i in dico_insti_tperiods_players[t_period].keys():
        is_playing = dico_insti_tperiods_players[t_period][player_i]\
                                                ["is_playing"][-1]
        if is_playing:
            probas_p_0s1s = [dico_chosen_strats_t_r[player_i]["p_0s"],
                             dico_chosen_strats_t_r[player_i]["p_1s"]]
            state_i, mode_i, prod_i, cons_i, r_i, S_t_plus_1_i \
                = compute_StateModeProdConsRiSi(
                    dico_insti_tperiods_players=dico_insti_tperiods_players,
                    t_period=t_period,
                    player_i=player_i,
                    Pi=dico_chosen_strats_t_r[player_i]["Pi"],
                    probas_p_0s1s=probas_p_0s1s
                    )
            Out_sg_k += cons_i
            In_sg_k += prod_i
            
            q_plus_k_i, q_minus_k_i \
                = fct_aux.compute_qi_plus_minus(
                    Pi=dico_chosen_strats_t_r[player_i]["Pi"],
                    Ci=dico_insti_tperiods_players[t_period][player_i]["Ci"], 
                    Si=dico_insti_tperiods_players[t_period][player_i]["S_t_plus_1_is"][-1],
                    Si_max=dico_insti_tperiods_players[t_period][player_i]["Si_max"]
                    )
            q_plus_k += q_plus_k_i
            q_minus_k += q_minus_k_i
            
            P_t_plus_1_i, C_t_plus_1_i \
                = get_Pi_Ci_t_plus_1(
                    dico_insti_tperiods_players=dico_insti_tperiods_players, 
                    t_periods=t_periods, t_period=t_period, 
                    player_i=player_i)
                
            variables = [("mode_i", mode_i), ("state_i", state_i),
                         ("q_plus_k_i", q_plus_k_i), ("q_minus_k_i", q_minus_k_i),
                         ("P_t_plus_1_i", P_t_plus_1_i), ("C_t_plus_1_i", C_t_plus_1_i),
                         ("prod_i", prod_i), ("cons_i", cons_i), ("r_i", r_i), 
                         ("S_t_plus_1_i", S_t_plus_1_i)]
            for var, val in variables:
                dico_chosen_strats_t_r[player_i][var] = val
                
                
    return q_plus_k, q_minus_k, \
            Out_sg_k, In_sg_k, \
            dico_chosen_strats_k 
        
        
#______________________________________________________________________________
#                       choose strategy for all players: fin
#______________________________________________________________________________

#______________________________________________________________________________
#          compute 
#           S_minus_i, S_plus_i, I_m, I_M, O_m, O_M 
#           ben_i, cst_i 
#           debut
#______________________________________________________________________________

def compute_gamma_Siplusminus_ImM_OmM_bencstis(dico_chosen_strats_t_r, 
                                               dico_insti_tperiods_players, 
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
    dico_chosen_strats_t_r[player_i] = {"f12_strategy_name_i":, 
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
        dico_t_player_i = 
            {keys=csts.LEARNING_PROPERTIES, 
            values=[values' list of LEARNING_PROPERTIES],
            r_0_s = [r_0_1, r_0_2, ..., r_0_100]
            r_1_s = [r_1_1, r_1_2, ..., r_1_100]
            r_2_s = [r_2_1, r_2_2, ..., r_2_100]
            f12_strategies_names = ["r_0_s", "r_1_s", "r_2_s"]
            f12_strategies_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]}
            }
                                   
    NB : if  S_plus_i == S_minus_i then frac \in [0.01, 0.1]
    """
    I_m, I_M, O_m, O_M = 0, 0, 0, 0
    for player_i, dico_strats in dico_chosen_strats_t_r.items():
        state_i = dico_chosen_strats_t_r[player_i]['state_i']
        Pi = dico_chosen_strats_t_r[player_i]['Pi']
        Ci = dico_insti_tperiods_players[t_period][player_i]['Ci']
        Si = dico_insti_tperiods_players[t_period][player_i]['Si']
        Si_max = dico_insti_tperiods_players[t_period][player_i]['Si_max']
        
        S_plus_i, S_minus_i, Ii_m, Ii_M, Oi_m, Oi_M \
            = fct_aux.compute_Si_plusminus_ImM_OmM(state_i=state_i, 
                                            Si_max=Si_max, 
                                            Si=Si, Ci=Ci, Pi=Pi)
        # compute I_m, I_M, O_m, O_M
        I_m = I_m if Ii_m is None else I_m + Ii_m
        I_M = I_M if Ii_M is None else I_M + Ii_M
        O_m = O_m if Oi_m is None else O_m + Oi_m
        O_M = O_M if Oi_M is None else O_M + Oi_M
        print(f"{player_i}: O_m={O_m}, O_M={O_M}, Pi={Pi}, Ci={Ci}, Si={Si}, {state_i}")
        
        # verification condition 
        assert S_minus_i <= S_minus_i
        
        # compute ppi_t
        P_t_plus_1_i = dico_chosen_strats_t_r[player_i]["P_t_plus_1_i"]
        C_t_plus_1_i = dico_chosen_strats_t_r[player_i]["C_t_plus_1_i"]
        if pp_t_i is None:
            frac \
                = fct_aux.diff_positive(
                    fct_aux.diff_positive(C_t_plus_1_i, P_t_plus_1_i), 
                    S_minus_i) / (S_plus_i - S_minus_i) \
                if S_minus_i != S_minus_i \
                else np.random.uniform(low=0.01, high=0.1)
            pp_t_i = math.sqrt( min(1, frac) )
                
        # compute gamma
        rd_draw = np.random.uniform(low=0.0, high=1.0, size=None)
        rho_i_t = 1 if rd_draw < pp_t_i else 0
        X_gamV5 = None
        X_gamV5 = pi_0_minus \
                    if state_i == csts.STATES[0] or state_i == csts.STATES[1] \
                    else pi_0_plus
        gamma_i = rho_i_t * (X_gamV5 + 1)
        dico_chosen_strats_t_r[player_i]["gamma_i"] = gamma_i
        dico_chosen_strats_t_r[player_i]["S_plus_i"] = S_plus_i
        dico_chosen_strats_t_r[player_i]["S_minus_i"] = S_minus_i
        dico_chosen_strats_t_r[player_i]["pp_t_i"] = pp_t_i
        
        # compute ben_i, cst_i
        prod_i = dico_chosen_strats_t_r[player_i]["prod_i"]
        cons_i = dico_chosen_strats_t_r[player_i]["cons_i"]
        r_i = dico_chosen_strats_t_r[player_i]["r_i"]
        ben_i = b0 * prod_i + gamma_i * r_i
        cst_i = c0 * cons_i
        dico_chosen_strats_t_r[player_i]["ben_i"] = ben_i
        dico_chosen_strats_t_r[player_i]["cst_i"] = cst_i
        dico_chosen_strats_t_r[player_i]["c0"] = c0
        dico_chosen_strats_t_r[player_i]["b0"] = b0
        dico_chosen_strats_t_r[player_i]["pi_0_minus"] = pi_0_minus
        dico_chosen_strats_t_r[player_i]["pi_0_plus"] = pi_0_plus
        dico_chosen_strats_t_r[player_i]["beta_sg_t_minus_1_plus"] = beta_sg_t_minus_1_plus
        dico_chosen_strats_t_r[player_i]["beta_sg_t_minus_1_minus"] = beta_sg_t_minus_1_minus
        
    return dico_chosen_strats_t_r, I_m, I_M, O_m, O_M
        
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
def compute_bg(dico_chosen_strats_t_r, dico_insti_tperiods_players, t_period, c0_M):
    """
    determine bg, bg_min_i, bg_max_i
    dico_chosen_strats_t_r[player_i] = {"strategy_name_i":, 
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
    dico_insti_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES], 
                           r_0_s = [r_0_1, r_0_2, ..., r_0_100]
                            r_1_s = [r_1_1, r_1_2, ..., r_1_100]
                            r_2_s = [r_2_1, r_2_2, ..., r_2_100]
                            f12_strategies_names = ["r_0_s", "r_1_s", "r_2_s"]
                            f12_strategies_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]}
                            }
                                  
    """
    for player_i in dico_chosen_strats_t_r.keys():
        Pi = dico_chosen_strats_t_r[player_i]['Pi']
        Ci = dico_insti_tperiods_players[t_period][player_i]['Ci']
        
        bg_i = dico_chosen_strats_t_r[player_i]['ben_i'] \
                - dico_chosen_strats_k[player_i]['cst_i'] \
                + c0_M * fct_aux.diff_positive(op1=Ci, op2=Pi)
                
        dico_chosen_strats_t_r[player_i]["bg_i"] = bg_i
        
        bg_min_i, bg_max_i = None, None
        bg_values = dico_insti_tperiods_players[t_period][player_i]\
                                           .get(csts.BG_ALL_VALUES)
                                           
        # TODO A SUPPRIMER LE commentaire PRINT
        #print(f"{player_i} bg_values = {bg_values}")
        
        bg_min_i = min(bg_values) if len(bg_values) != 0 else None
        bg_max_i = max(bg_values) if len(bg_values) != 0 else None
        
        if bg_min_i is None or bg_min_i > bg_i:
            bg_min_i = bg_i
            
        if bg_max_i is None or bg_max_i < bg_i:
            bg_max_i = bg_i
        
        dico_chosen_strats_t_r[player_i]["bg_min_i"] = bg_min_i
        dico_chosen_strats_t_r[player_i]["bg_max_i"] = bg_max_i
        
    return dico_chosen_strats_t_r

#______________________________________________________________________________
#          compute bg
#           fin
#______________________________________________________________________________  


#______________________________________________________________________________
#          update p_Xs
#           debut
#______________________________________________________________________________

def update_probabilities_pXs(dico_chosen_strats_t_r, 
                             dico_insti_tperiods_players, 
                             t_period, 
                             learning_rate):
    """
    update probability for choosen strategy
    
    dico_chosen_strats_k[player_i] = {"q_plus_k_i":, "q_minus_k_i":,
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
    dico_insti_tperiods_players = {"t_j":{'player_i':dico_t_player_i, 
                                        'player_i+1':dico_t_player_i+1,
                                        ...,}, 
                             ...
                             }
        avec 
        t_period = 't_j' with j = {0,1,..., N}
        player_ = 'player_i' with i = {0,1,..., M} and 
        dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                           values=[values' list of LEARNING_PROPERTIES]
                           r_0_s = [r_0_1, r_0_2, ..., r_0_100]
                            r_1_s = [r_1_1, r_1_2, ..., r_1_100]
                            r_2_s = [r_2_1, r_2_2, ..., r_2_100]
                            f12_strategies_names = ["r_0_s", "r_1_s", "r_2_s"]
                            f12_strategies_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]        
                        }
    
    """
    is_all_strategy_probas_sup_SEUIL = False
    cpt_all_strategy_probas = 0
    for player_i in dico_chosen_strats_t_r.keys():
        u_i = dico_chosen_strats_t_r[player_i]["u_i"]
        f12_strategy_name_i = dico_chosen_strats_t_r[player_i]["f12_strategy_name_i"]
        f22_strategy_name_i = "p_0s"
        p_X = dico_insti_tperiods_players[t_period][player_i]\
                [f12_strategy_name_i][f22_strategy_name_i][-1]
        
        increasing_val = learning_rate * u_i * (1 - p_X)
        # p_X_new = round(p_X + increasing_val, csts.ARRONDI)
        p_X_new = p_X + increasing_val
        dico_chosen_strats_t_r[player_i][f22_strategy_name_i] = p_X_new
        
        is_playing = True if p_X_new < csts.MAX_STRATEGY_PROBA else False
        dico_chosen_strats_t_r[player_i]["is_playing"] = is_playing
        
        cpt_all_strategy_probas += 1 if p_X_new >= csts.MAX_STRATEGY_PROBA else 0
        
        print(f"------> {player_i} f12_strategy_name_i {f12_strategy_name_i}, f22_strategy_name_i p_Os= {p_X_new}")
        
        f22_strategy_names \
                = set(
                    filter(lambda x: x.startswith(csts.PROB_NAME), 
                           list(dico_insti_tperiods_players\
                                [t_period][player_i]\
                                [f12_strategy_name_i].keys())
                          )
                    )
                
        '''
        count the negative probability
        '''
        dico_neg_Px_name = dict()
        no_choose_strategy_names = f22_strategy_names - set({f22_strategy_name_i})
        for no_choose_strategy_name in no_choose_strategy_names:
            p_X = dico_insti_tperiods_players[t_period][player_i]\
                                             [f12_strategy_name_i]\
                                             [no_choose_strategy_name]
            if p_X - increasing_val < 0:
                dico_neg_Px_name[no_choose_strategy_name] = p_X
        
        '''
        update probas for no_choose_strategy_names
        '''
        Y_prob_sup_increaseVal = f22_strategy_names \
                                    - {f22_strategy_name_i} \
                                    - set(dico_neg_Px_name.keys())
        for no_choose_strategy_name_i in no_choose_strategy_names:
            p_y_new = None
            if no_choose_strategy_name_i in dico_neg_Px_name.keys():
                p_y_new = 0
            else:
                p_y = dico_insti_tperiods_players[t_period][player_i]\
                                            [f12_strategy_name_i]\
                                            [no_choose_strategy_name_i]
                # p_y_new = round(p_y \
                #                 - increasing_val/len(Y_prob_sup_increaseVal) \
                #                 + sum(dico_neg_Px_name.values()) / \
                #                     len(Y_prob_sup_increaseVal), csts.ARRONDI)
                p_y_new = p_y \
                            - increasing_val/len(Y_prob_sup_increaseVal) \
                            + sum(dico_neg_Px_name.values()) / \
                                len(Y_prob_sup_increaseVal)
                pass
            
            dico_chosen_strats_t_r[player_i][no_choose_strategy_name_i] = p_y_new
            
        '''
        print proba by pX_name
        '''
        sum_probas = 0
        for f22_strategy_name in f22_strategy_names:
            prob_pX_name = dico_chosen_strats_t_r[player_i][f22_strategy_name]
            print(f"{f22_strategy_name} p_0s: {prob_pX_name}") \
                if f22_strategy_name != f22_strategy_name_i \
                else print(f"{f22_strategy_name}: {prob_pX_name} <----")
            sum_probas += prob_pX_name
            
        print(f"sum_probas = {sum_probas}")
        
            
    is_all_strategy_probas_sup_SEUIL = True \
        if cpt_all_strategy_probas >= len(dico_chosen_strats_k.keys()) \
        else False
    
    return dico_chosen_strats_t_r, is_all_strategy_probas_sup_SEUIL

#______________________________________________________________________________
#          update p_Xs
#           fin
#______________________________________________________________________________ 

#______________________________________________________________________________
#          execute one learning step
#           debut
#______________________________________________________________________________
def f12_f22_execute_one_learning_step_4_one_period(dico_insti_tperiods_players, 
                                                    dico_chosen_strats_t_r,
                                                    t_period, 
                                                    k_step,
                                                    r_step,
                                                    cpt_repeated_kstep,
                                                    args):
    """
    run one learning step in one period ie
        * balancing each player
        * compute pi_EPO_plus, pi_EPO_minus
        * compute beta_sg_t_minus_1_plus, beta_sg_t_minus_1_minus
        * compute pi_0_plus, pi_0_minus
        * compute b0, c0
        * compute gamma, Siplusminus, ImM, OmM, ben_i, cst_i
        * compute bg
        * compute utility fonction u_i
        * update probabilities p_Xs, X={0,1}
        * update learning variables
    """
    
    q_plus_k, q_minus_k, \
    Out_sg_k, In_sg_k, \
    dico_chosen_strats_t_r_k \
        = choose_strategy_4_all_players(
            dico_insti_tperiods_players=dico_insti_tperiods_players, 
            dico_chosen_strats_t_r=dico_chosen_strats_t_r,
            t_period=t_period
            )
        
    # compute q_t_plus, q_t_minus
    q_plus_k = 0 if q_plus_k < 0 else q_plus_k
    q_minus_k = 0 if q_minus_k < 0 else q_minus_k
    phi_EPO_plus = fct_aux.compute_phi_EPO_plus(q_plus_k=q_plus_k, 
                                        quantity=args["quantity_a"], 
                                        a=args["a"])
    phi_EPO_minus = fct_aux.compute_phi_EPO_minus(q_minus_k=q_minus_k, 
                                          quantity=args["quantity_b"], 
                                          b=args["b"])
    pi_EPO_plus = round(phi_EPO_plus/q_plus_k, csts.ARRONDI) \
                    if q_plus_k > 0 \
                    else 0
    pi_EPO_minus = round(phi_EPO_minus/q_minus_k, csts.ARRONDI) \
                    if q_minus_k > 0 \
                    else 0
    
    beta_sg_t_minus_1_plus, beta_sg_t_minus_1_minus \
        = fct_aux.compute_beta_sg(t_period=t_period, 
                          dico_chosen_strats_t=dico_chosen_strats_t_r,
                          beta_sg_0_minus=pi_EPO_minus-1, 
                          beta_sg_0_plus=pi_EPO_plus-1,
                          quantity_a=args["quantity_a"],  
                          a=args["a"],  
                          quantity_b=args["quantity_b"],        
                          b=args["b"]
                          )

    pi_0_plus, pi_0_minus \
        = fct_aux.compute_pi_0(t_period=t_period, 
                       beta_sg_t_minus_1_plus=beta_sg_t_minus_1_plus, 
                       beta_sg_t_minus_1_minus=beta_sg_t_minus_1_minus,
                       pi_EPO_plus=pi_EPO_plus, 
                       pi_EPO_minus=pi_EPO_minus
                       )

    b0, c0 = fct_aux.compute_b0_c0(pi_0_plus=pi_0_plus, 
                       pi_0_minus=pi_0_minus, 
                       Out_sg_k=Out_sg_k, In_sg_k=In_sg_k, 
                       quantity_a=args["quantity_a"], a=args["a"], 
                       quantity_b=args["quantity_b"], b=args["b"])
    
    # compute gamma_i, Si_plus, Si_minus, pp_i for all players
    # compute also I_m, I_M, O_m, O_M 
    # compute also ben_is, cst_is
    dico_chosen_strats_t_r, I_m, I_M, O_m, O_M \
        = compute_gamma_Siplusminus_ImM_OmM_bencstis(
            dico_chosen_strats_t_r=dico_chosen_strats_t_r, 
            dico_insti_tperiods_players=dico_insti_tperiods_players, 
            t_period=t_period, 
            pp_t_i=args["pp_t_i"], 
            pi_0_minus=pi_0_minus, pi_0_plus=pi_0_plus, 
            beta_sg_t_minus_1_plus=beta_sg_t_minus_1_plus, 
            beta_sg_t_minus_1_minus=beta_sg_t_minus_1_minus,
            b0=b0, c0=c0
            )
        
    assert In_sg_k >= I_m and In_sg_k <= I_M
    print(f"Out_sg_k={Out_sg_k}, O_m={O_m}, O_M={O_M}")
    assert Out_sg_k >= O_m and Out_sg_k <= O_M

    # compute c0_M
    frac_cOM = round(((O_M - I_m)*pi_EPO_minus + I_M*pi_0_minus)/O_m, 
                      csts.ARRONDI)
    # TODO ----> IMPORTANT a voir avec Dominique c0_M = MIN(frac_cOM, pi_0_minus)
    c0_M = max(frac_cOM, pi_0_minus)
    c0_M = round(c0_M)
    print(f'c0={c0}, c0_M={c0_M}, frac_cOM={frac_cOM}')
    assert c0 <= c0_M
    
    # compute bg for all players
    # add fields bg_min_i and bg_max_i in dico_chosen_strats_k
    dico_chosen_strats_t_r \
        = compute_bg(dico_chosen_strats_t_r=dico_chosen_strats_t_r, 
                     dico_insti_tperiods_players=dico_insti_tperiods_players, 
                     t_period=t_period, 
                     c0_M=c0_M)
    
    
    # compute utility fonction
    dico_chosen_strats_t_r, bool_is_One_Ui_None \
        = fct_aux.compute_utility_fonction_ui(
            dico_chosen_strats_k=dico_chosen_strats_t_r, 
            t_period=t_period
            )
        
    # update probabilities in case of all u_i are different to None
    is_all_strategy_probas_sup_SEUIL = False
    if bool_is_One_Ui_None == False \
        or cpt_repeated_kstep >= csts.MAX_REPEATED_KSTEP-1:
        dico_chosen_strats_k, is_all_strategy_probas_sup_SEUIL \
            = update_probabilities_pXs(
                dico_chosen_strats_t_r=dico_chosen_strats_t_r, 
                dico_insti_tperiods_players=dico_insti_tperiods_players, 
                t_period=t_period, 
                learning_rate=args['learning_rate'])
        
        # update startegy keys and other variables to dico_tperiods_players
        is_repeated_kstep = False
        dico_insti_tperiods_players \
            = update_learning_variables(
                dico_chosen_strats_t_r=dico_chosen_strats_t_r, 
                dico_insti_tperiods_players=dico_insti_tperiods_players, 
                t_period=t_period, 
                k_step=k_step,
                is_repeated_kstep=is_repeated_kstep)
            
        if is_all_strategy_probas_sup_SEUIL:
            k_step = args['k_step']
        else:
            k_step += 1
        cpt_repeated_kstep = csts.MAX_REPEATED_KSTEP
        
    else:
        """ else means bool_is_One_Ui_None == True \
                        and cpt_repeated_kstep < csts.MAX_REPEATED_KSTEP
        """
        
        # update startegy keys and other variables to dico_tperiods_players
        is_repeated_kstep = True
        dico_insti_tperiods_players \
            = update_learning_variables(
                dico_chosen_strats_t_r=dico_chosen_strats_t_r, 
                dico_insti_tperiods_players=dico_insti_tperiods_players, 
                t_period=t_period, 
                k_step=k_step,
                is_repeated_kstep=is_repeated_kstep)
        
        k_step = k_step
        cpt_repeated_kstep += 1
        
    dico_chosen_strats_t[t_period] = dico_chosen_strats_t_r
    
    return dico_insti_tperiods_players, \
            dico_chosen_strats_t, \
            k_step, cpt_repeated_kstep, \
            is_all_strategy_probas_sup_SEUIL
        
    
    
    
    
#______________________________________________________________________________
#          execute one learning step
#           fin
#______________________________________________________________________________
    
#______________________________________________________________________________
#                  compute prod, cons, Ri, Si for all players: debut
#______________________________________________________________________________  
def compute_ProdConsRiSi_4_players(dico_tperiods_players,
                                   dico_chosen_strats_k, 
                                   t_period):
    """
    """
    for player_i in dico_tperiods_players[t_period].keys():
        f2_strategy_name_i = dico_chosen_strats_k[player_i]["f2_strategy_name_i"]
        Pi = dico_tperiods_players[t_period][player_i]\
                ["f2_strategies"][f2_strategy_name_i]["Pi"]
        state_i = dico_tperiods_players[t_period][player_i]\
                ["f2_strategies"][f2_strategy_name_i]["state_i"]
        Ci =  dico_tperiods_players[t_period][player_i]["Ci"]
        
        # TODO : creer une liste pour les proba de S1 et S2 
        # 3 je ne sais pas ou mettre ces 2 listes 
        
#______________________________________________________________________________
#                  compute prod, cons, Ri, Si for all players: fin
#______________________________________________________________________________  

#______________________________________________________________________________
#                       game extension 2 - f21 & f22 : debut
#______________________________________________________________________________
def game_ext2(dico_insti_tperiods_players, args):
    
    """
    * add probabilities for all f2_strategies : 
            already insert in the generate_strategies_ext2 function
    * for r in 0 to 100:
            * select f21-strategy by probabilities for all players
            * for k_step in 0 to 1000:
                    * play initial game f22
                    * at the end of k_stop due to stabilisation, update 
                        players' probabilities in f21_strategies 
    """
    for t_period in dico_insti_tperiods_players.keys():
        print(f'**** period : {t_period} ****')
        
        r_step = 0
        while r_step < args['r_steps']:
            dico_chosen_strats_t_r \
                = choose_f12_strategy_4_all_players(
                    dico_insti_tperiods_players=dico_insti_tperiods_players, 
                    t_period=t_period, 
                    r_step=r_step)
                
            k_step = 0; cpt_repeated_kstep = 0
            while k_step < args["k_steps"] \
                    or cpt_repeated_kstep < csts.MAX_REPEATED_KSTEP:
                k_step_res = None
                
                dico_tperiods_players, dico_chosen_strats_t_r, \
                k_step_res, cpt_repeated_kstep, \
                is_all_f22_strategy_probas_sup_SEUIL \
                    = f12_f22_execute_one_learning_step_4_one_period(
                        dico_insti_tperiods_players=dico_insti_tperiods_players, 
                        dico_chosen_strats_t_r=dico_chosen_strats_t_r,
                        t_period=t_period, 
                        k_step=k_step,
                        r_step=r_step,
                        cpt_repeated_kstep=cpt_repeated_kstep,
                        args=args
                        )
                k_step = k_step_res
            
            # don't forget to update P_is at ecah r_step
            
            # incrementation of r_step
            r_step += 1
            
            
                                   
        for r in range(0, args['r_steps']):
            dico_chosen_strats_k \
                = choose_f12_strategy_4_all_players(
                    dico_insti_tperiods_players=dico_insti_tperiods_players, 
                    t_period=t_period)
                
            k_step = 0; cpt_repeated_kstep = 0
            while k_step < args["k_steps"] \
                    or cpt_repeated_kstep < csts.MAX_REPEATED_KSTEP:
                   
                dico_chosen_strats_k, \
                prod_i, cons_i, r_i, S_t_plus_1_i \
                    = compute_ProdConsRiSi_4_players(
                        dico_insti_tperiods_players=dico_insti_tperiods_players,
                        dico_chosen_strats_k=dico_chosen_strats_k,
                        t_period=t_period
                        )
    
    pass

#______________________________________________________________________________
#                       game extension 2 - f21 & f22: fin
#______________________________________________________________________________

#______________________________________________________________________________
#                       test : debut
#______________________________________________________________________________
def read_json_dataset_initialise_strategies(args):
    """
    read and initialise strategies to all players
    """
    dico_insti_tperiods_players \
        = gene_strats.read_json_filename(json_filename=args['json_filename'])
        
    dico_insti_tperiods_players \
        = gene_strats.add_fields_players(dico_data=dico_insti_tperiods_players)
        
    dico_insti_tperiods_players \
        = gene_strats.generate_strategies_f21_f22_players(
            dico_insti_tperiods_players=dico_insti_tperiods_players, 
            r_steps=args["r_steps"])
        
    return dico_insti_tperiods_players

def test_chosen_strats_r(args):
    dico_insti_tperiods_players = read_json_dataset_initialise_strategies(args)
    t_period = np.random.choice(range(args["t_periods"]))
    t_period = csts.PERIOD_ROOT + str(t_period)
    
    dico_chosen_strats_k \
        = choose_f12_strategy_4_all_players(
            dico_insti_tperiods_players=dico_insti_tperiods_players, 
            t_period=t_period,
            r_step=args["r_steps"])
        
    return dico_chosen_strats_k

def test_gameExt2(args):
    
    dico_insti_tperiods_players = read_json_dataset_initialise_strategies(args) 
    

    return dico_insti_tperiods_players
    # # test generations strategies
    # dico_insti_tperiods_players \
    #     = gene_strats.generation_strategies_tperiods(
    #         dico_data=dico_insti_tperiods_players, 
    #         t_periods=args['t_periods'])
        
    # game_ext1(dico_tperiods_players=dico_insti_tperiods_players, 
    #           args=args)
    pass
#______________________________________________________________________________
#                       test : fin
#______________________________________________________________________________

if __name__ == "__main__":
    ti = time.time()
    t_periods = 5 #50 #15
    k_steps = 1000 #5000 #1000
    r_steps = 100
    n_instances = 10
    n_instance = np.random.choice(range(10))
    
    json_filename = f"data_instance_{n_instance}_tperiods_{t_periods}.json"
    
    quantity_a = 10; a = 1;
    quantity_b = 30; b = 1;
    pp_t_i = None
    learning_rate = 0.1
    path_2_save = os.path.join(".", csts.PATH_ROOT, 
                               csts.INSTANCE_ROOT+str(n_instance))
    
    args = {"json_filename": json_filename, 
            "t_periods": t_periods,
            "n_instance": n_instance, 
            "k_steps": k_steps,
            "r_steps": r_steps,
            "quantity_a": quantity_a, "a": a,
            "quantity_b": quantity_b, "b": b,
            "pp_t_i": pp_t_i, 
            "learning_rate": learning_rate, 
            "path_2_save": path_2_save
            }
    
    dico_insti_tperiods_players = test_gameExt2(args)
    # test chosen strategies
    dico_chosen_strats_k = test_chosen_strats_r(args)
    
    print(f"runtime = {time.time() - ti}")