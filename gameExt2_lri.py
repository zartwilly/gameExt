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
#                       choose strategy for all players: debut
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
    dico_chosen_strats_r = dict()     # key : player_X, values = dict()
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
            dico_prob_Pi_p_01s["f12_strategy_name_i"] = eval(f12_strategy_name_i)
            
                
            dico_chosen_strats_r[player_i] = dico_prob_Pi_p_01s
            
    return dico_chosen_strats_r
#______________________________________________________________________________
#                       choose strategy for all players: fin
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
            dico_chosen_strats_r \
                = choose_f12_strategy_4_all_players(
                    dico_insti_tperiods_players=dico_insti_tperiods_players, 
                    t_period=t_period, 
                    r_step=r_step)
                
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