
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 13:11:51 2021

@author: willy

read json file corresponding to one instance
select one period
    * there are a list of players
    * each player is respresented to a dictionary
add fields to each player (ie Sis, p1_ijks, p2_ijks, ksteps, prod_is, cons_is, etc...)

"""

import os
import json
import h5py
import time
import itertools as it
import numpy as np

import constances as csts

#______________________________________________________________________________
#                       read json filename: debut 
#______________________________________________________________________________
def read_json_filename(json_filename):
    f = open(os.path.join(csts.DATA_ROOT, json_filename))
    dico_data = eval(json.load(f))
     
    return dico_data

def add_fields_players(dico_data):
    """
    dico_data : {"t_0": {'player_0':{}, 'player_1':{}}, 
                 "t_1": {'player_0':{}, 'player_1':{}} }

    """
    for key_t, val_players in dico_data.items():
        for player_name, dico_values in val_players.items():
            for propertie in csts.LEARNING_PROPERTIES:
                if propertie == "S_t_plus_1_is":
                    dico_values[propertie] = [dico_values["Si"]]
                elif propertie == "is_playing":
                    dico_values[propertie] = [True]
                elif propertie == csts.BG_ALL_VALUES:
                    dico_values[propertie] = set()
                else:
                    dico_values[propertie] = []
                
    return dico_data
                
#______________________________________________________________________________
#                       read json filename: fin 
#______________________________________________________________________________


#______________________________________________________________________________
#                       generation of strategies: debut 
#______________________________________________________________________________
def compute_strategies_onePlayer(dico_values, m_players):
    """
    this fonction adds one field "strategies" to dico_values
    for example:
    dico_values = {'Ci':15, 'z':11, 'di':[0,5,11], 'Si':0, 'Si_max':20, ... }
    """
    possibles_strategies = []
    for Pi in dico_values["di"]:
        state_i = None
        if Pi + dico_values['Si'] <= dico_values['Ci']:
            state_i = csts.STATES[0]
            modes = csts.STATE1_STRATS
            possibles_strategies.extend( list(it.product([Pi], [state_i], modes)) )
        elif Pi + dico_values['Si'] > dico_values['Ci'] \
            and Pi <= dico_values['Ci']:
            state_i = csts.STATES[1]
            modes = csts.STATE2_STRATS
            possibles_strategies.extend( list(it.product([Pi], [state_i], modes)) )
        elif Pi >= dico_values['Ci']:
            state_i = csts.STATES[2]
            modes = csts.STATE3_STRATS
            possibles_strategies.extend( list(it.product([Pi], [state_i], modes)) )
            
    # compute probability of each player
    dico_possibles_strategies_probs = dict()
    strategies_names = set()
    for i, tuple_comb in enumerate(possibles_strategies):
        dico_values[csts.PROB_NAME+str(i)+"s"] \
            = [ 1 / len(possibles_strategies) ]
    
        strategies_names.add(csts.PROB_NAME+str(i)+"s")
        dico_possibles_strategies_probs[csts.PROB_NAME+str(i)+"s"] \
            = {"Pi": tuple_comb[0], "state_i": tuple_comb[1], 
               "mode_i": tuple_comb[2]}
        
    dico_values["strategies"] = dico_possibles_strategies_probs
    dico_values["strategies_names"] = sorted(strategies_names)
    return dico_values
            
def generation_strategies(dico_data, t_period=2):
    """
    dico_data : {"t_0": {'player_0':{}, 'player_1':{}}, 
                 "t_1": {'player_0':{}, 'player_1':{}} }

    """
    dico_data_t = dico_data[csts.PERIOD_ROOT+str(t_period)]
    m_players = len( list(dico_data_t.keys()) )
    for player_name, dico_values in dico_data_t.items():
        dico_values = compute_strategies_onePlayer(
                        dico_values=dico_values, 
                        m_players=m_players)
        
    return  dico_data_t
        
def generation_strategies_tperiods(dico_data, t_periods):
    """
    dico_data : {"t_0": {'player_0':{}, 'player_1':{}}, 
                 "t_1": {'player_0':{}, 'player_1':{}} }

    """
    for t_period in range(0, t_periods):
        dico_data[csts.PERIOD_ROOT+str(t_period)] \
            = generation_strategies(dico_data=dico_data, t_period=t_period)
        
    return dico_data

#______________________________________________________________________________
#                       generation of strategies: fin 
#______________________________________________________________________________


#______________________________________________________________________________
#                       generate strategies f12 & f22 : debut
#______________________________________________________________________________
def add_f22_strategies(dico_t_player_i, r_steps):
    """
    add f22 strategy for one player

    Parameters
    ----------
    dico_player_i : TYPE
        DESCRIPTION.
    r_step : integer
        learning step number for production_game2

    csts.F21_NAME = 'r_' 

    Returns
    -------
    dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                        values=[values' list of LEARNING_PROPERTIES],
                         r_0_s = [r_0_1, r_0_2, ..., r_0_100]
                         r_1_s = [r_1_1, r_1_2, ..., r_1_100]
                         r_2_s = [r_2_1, r_2_2, ..., r_2_100]
                         strategies_f22_names = ["r_0_s", "r_1_s", "r_2_s"]
                         strategies_f22_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]},

    """
    probs = [1/len(dico_t_player_i["di"]) for i in dico_t_player_i["di"]]
    probs = np.array(probs)
    probs /= probs.sum()
    probs = list(probs)
    
    f12_strategies = set()
    for u, di in enumerate(dico_t_player_i["di"]):
        prob = probs[u]
        for r_step in range(r_steps):
            r_u_rstep_s = {"prob":prob, "di": di, 
                           "strategies_f22_names":[csts.PROB_NAME+"0s", 
                                                   csts.PROB_NAME+"1s"],
                           csts.PROB_NAME+"0s": [0.5],
                           csts.PROB_NAME+"1s": [0.5]
                           } 
            dico_t_player_i[csts.F21_NAME+str(u)+"_s"] = r_u_rstep_s
        f12_strategies.add(csts.F21_NAME+str(u)+"_s")
      
    dico_t_player_i["strategies_f12_names"] = list(f12_strategies)
    dico_t_player_i["strategies_f12_names_is"] = []
    
    return dico_t_player_i
    
    
def generate_strategies_f21_f22_players(dico_insti_tperiods_players, r_steps):
    """
    generate strategies for production_game2 (f12) and initial_game2 (f22).
    f12 runs during r steps
    f22 runs during k steps
    
    * f12 strategy
    there are 3 strategies r_u_s = [r_u_1, ..., r_u_100] 
    with $u \in {0,1,2}$ and 
    $r_u_1 = {"prob":val1, "di":val2, "strategies_f22_names":["p_0s","p_1s"], 
              "strategies_f22_names_is": ["p_0s"^1, ..., "p_1s"^1000],
               "p_0s": [prob_1, ..., prob_1000], 
               "p_1s": [1-prob_1, ..., 1-prob_1000], 
              }
    $
    * f22 strategy
    r_0_s = [r_0_1, r_0_2, ..., r_0_100]
    r_1_s = [r_1_1, r_1_2, ..., r_1_100]
    r_2_s = [r_2_1, r_2_2, ..., r_2_100]
    strategies_f22_names = ["r_0_s", "r_1_s", "r_2_s"]
    strategies_f22_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]

    Parameters
    ----------
    dico_insti_tperiods_players : TYPE
        DESCRIPTION.

    Returns
    -------
    dico_insti_tperiods_players = {
        "t_j":{'player_i':dico_t_player_i, 
               'player_i+1':dico_t_player_i+1,
        }
    
    with 
    t_period = 't_j' with j = {0,1,..., N}
    player_ = 'player_i' with i = {0,1,..., M} and 
    dico_t_player_i = {keys=csts.LEARNING_PROPERTIES, 
                        values=[values' list of LEARNING_PROPERTIES],
                         r_0_s = [r_0_1, r_0_2, ..., r_0_100]
                         r_1_s = [r_1_1, r_1_2, ..., r_1_100]
                         r_2_s = [r_2_1, r_2_2, ..., r_2_100]
                         strategies_f22_names = ["r_0_s", "r_1_s", "r_2_s"]
                         strategies_f22_names_is = ["r_0_s"^1, "r_0_s"^2, ...,"r_2_s"^100]},
                                   
        
    """
    
    for t_period in dico_insti_tperiods_players.keys():
        dico_t_players = dico_insti_tperiods_players[t_period]
        for player_i in dico_t_players.keys():
            dico_t_players[player_i] \
                = add_f22_strategies(dico_t_player_i=dico_t_players[player_i], 
                                     r_steps=r_steps)
        
    return dico_insti_tperiods_players
#______________________________________________________________________________
#                       generate strategies f12 & f22 : fin
#______________________________________________________________________________

if __name__ == "__main__":
    t_periods = 5
    r_steps = 10
    k_steps = 100
    n_instances = 10
    n_instance = np.random.choice(range(10))
    
    json_filename = f"data_instance_{n_instance}_tperiods_{t_periods}.json"
    dico_insti_tperiods_players \
        = read_json_filename(json_filename=json_filename)
        
    dico_insti_tperiods_players \
        = add_fields_players(dico_data=dico_insti_tperiods_players)
        
    # test generations strategies
    dico_insti_tperiods_players \
        = generation_strategies_tperiods(
            dico_data=dico_insti_tperiods_players, 
            t_periods=t_periods)
        
    # test generation f12 and f22 strategies
    dico_insti_tperiods_players_f12_22 \
        = generate_strategies_f21_f22_players(dico_insti_tperiods_players, 
                                              r_steps)
    
    
    