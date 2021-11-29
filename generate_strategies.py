
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
                dico_values[propertie] = []
                
    return dico_data
                
#______________________________________________________________________________
#                       read json filename: fin 
#______________________________________________________________________________

#______________________________________________________________________________
#                       generation of strategies: debut 
#______________________________________________________________________________
def compute_strategies_onePlayer(dico_values):
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
            
    dico_values["strategies"] = possibles_strategies
    
    return dico_values
            
def generation_strategies(dico_data, t_period=2):
    """
    dico_data : {"t_0": {'player_0':{}, 'player_1':{}}, 
                 "t_1": {'player_0':{}, 'player_1':{}} }

    """
    dico_data_t = dico_data[csts.PERIOD_ROOT+str(t_period)]
    for player_name, dico_values in dico_data_t.items():
        dico_values = compute_strategies_onePlayer(dico_values)
        
    return  dico_data_t
        
def generation_strategies_tperiods(dico_data, t_periods):
    """
    dico_data : {"t_0": {'player_0':{}, 'player_1':{}}, 
                 "t_1": {'player_0':{}, 'player_1':{}} }

    """
    for t_period in range(0, t_periods):
        dico_data[csts.PERIOD_ROOT+str(t_period)] \
            = generation_strategies(dico_data, t_period=t_period)
        
    return dico_data

#______________________________________________________________________________
#                       generation of strategies: fin 
#______________________________________________________________________________





if __name__ == "__main__":
    t_periods = 5
    n_instances = 10
    n_instance = np.random.choice(range(10))
    
    json_filename = f"data_instance_{n_instance}_tperiods_{t_periods}.json"
    dico_players_insti_tperiods \
        = read_json_filename(json_filename=json_filename)
        
    dico_players_insti_tperiods \
        = add_fields_players(dico_data=dico_players_insti_tperiods)
        
    # test generation strategies t_period=2
    dico_players_insti_tperiod \
        = generation_strategies(
            dico_data=dico_players_insti_tperiods, 
            t_period=2)
        
    # test generation strategies t_periods
    dico_players_insti_tperiods \
        = generation_strategies_tperiods(
            dico_data=dico_players_insti_tperiods, 
            t_periods=t_periods)
    