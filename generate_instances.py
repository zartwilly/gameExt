#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:09:49 2021

@author: willy

Generate n instances of game. 
Each instance has to be at one situation: A, B, C
"""
import os
import json
import numpy as np
import h5py
import time

import constances as csts

def get_situation(keys_situat, num_pl, m_A, m_B, m_C):
    """
    get situation from num_pl. 
    for example, we have 
        * m_A players in situation A
        * m_B players in situation B
        * m_C players in situation C
    remember that keys_situat = [Situatiom_A, Situatiom_B, Situatiom_C]
    """
    if num_pl < m_A:
        return keys_situat[0]
    elif num_pl >= m_A and num_pl < m_A+m_B:
        return keys_situat[1]
    elif num_pl >= m_A+m_B and num_pl < m_A+m_B+m_C:
        return keys_situat[2]
    else:
        return csts.SITUATION_ERR
    
def get_values_CidizSiSimax(dico_situat):
    """
    return values of Ci, z, di, Si, Si_max

    exple: dico_situat = {'Ci': [15], 'z': [5, 16], 'di': 3, 
                          'Si': [0], 'Si_max': [20]}
    """
    Ci = dico_situat["Ci"][0] \
            if len(dico_situat["Ci"]) == 1 \
            else np.random.randint(low=dico_situat["Ci"][0], 
                                    high=dico_situat["Ci"][1]+1)
    Si = dico_situat["Si"][0] \
        if len(dico_situat["Si"]) == 1 \
        else np.random.randint(low=dico_situat["Si"][0], 
                               high=dico_situat["Si"][1]+1)
    Si_max = dico_situat["Si_max"][0] \
        if len(dico_situat["Si_max"]) == 1 \
        else np.random.randint(low=dico_situat["Si_max"][0], 
                               high=dico_situat["Si_max"][1]+1)
    z = np.random.randint(low=dico_situat["z"][0], 
                          high=dico_situat["z"][1]+1)
    z_div_2 = int(np.floor(z/2))
    di = [0, z_div_2, z]
    
    return Ci, z, di, Si, Si_max

#______________________________________________________________________________
#               generate n instanaces with t_periods periods: debut 
#______________________________________________________________________________
def generate_one_period(data, situations, n_instance, t_period, 
                        m_A=10, m_B=8, m_C=7):
    cpt_A, cpt_B, cpt_C = 0, 0, 0
    m_players = m_A + m_B + m_C
    
    # generate instances
    dico_players = dict()
    for num_pl in range(0, m_players):
        key_situat = get_situation(keys_situat=situations, 
                                   num_pl=num_pl, m_A=m_A, m_B=m_B, m_C=m_C)
        Ci, z, di, Si, Si_max = None, None, None, None, None
        lock = True
        if key_situat == situations[0]:
            Ci, z, di, Si, Si_max = get_values_CidizSiSimax(
                                        dico_situat=data[key_situat]
                                        )
            cpt_A += 1
            lock = False \
                if Ci is None or z is None or di is None or Si is None or Si_max is None \
                else True
        elif key_situat == situations[1]:
            prob = np.random.randn()
            cases_situatB = list(data[key_situat].keys())
            cases_situatB.sort(reverse=False)
            key_sitB_05 = cases_situatB[0] if prob<0.5 else cases_situatB[1]
            Ci, z, di, Si, Si_max = get_values_CidizSiSimax(
                                        dico_situat=data[key_situat][key_sitB_05]
                                        )
            cpt_B += 1
            lock = False \
                if Ci is None or z is None or di is None or Si is None or Si_max is None \
                else True
        elif key_situat == situations[2]:
            Ci, z, di, Si, Si_max = get_values_CidizSiSimax(
                                        dico_situat=data[key_situat]
                                        )
            cpt_C += 1
            lock = False \
                if Ci is None or z is None or di is None or Si is None or Si_max is None \
                else True
        else:
            print(f"---inst={num_pl},t={t_period}: key={key_situat} Unknown ----")
            
        dico_players[csts.PLAYER_ROOT+str(num_pl)] = \
                {"Ci":Ci, "z":z, "di":di, "Si":Si, "Si_max":Si_max, 
                 "Situation":key_situat.split("_")[1]} \
                    if lock \
                    else None
            
    print(f"test_instance_{n_instance}_t_{t_period}_generated --> OK") \
        if m_players == cpt_A + cpt_B + cpt_C \
        else print(f"test_instance_{n_instance}_t_{t_period}_generated --> NOK: "+
                   f"m_A={m_A} != cpt_A={cpt_A}, m_B={m_B} != cpt_B={cpt_B},"+
                   f" m_C={m_C} != cpt_C={cpt_C}")
        
    return dico_players
    
def generate_one_instance_tperiods(data, situations, n_instance, t_periods,
                                   m_A, m_B, m_C):
    dico_players_insti_tperiod = dict()
    for t_period in range(0, t_periods):
        dico_players_insti_tperiod[t_period] \
            = generate_one_period(data=data, situations=situations, 
                                  n_instance=n_instance, 
                                  t_period=t_period,
                                  m_A=m_A, m_B=m_B, m_C=m_C)
    return dico_players_insti_tperiod


def generate_n_instances_t_periods(file_name="data1_generation.json", 
                                   n_instances=10, t_periods=5,
                                   m_A=10, m_B=8, m_C=7, 
                                   path_2_save=""):
     
    # read json file
    f = open(os.path.join(csts.DATA_ROOT, file_name))
    data = json.load(f)
    situations = list(data.keys()); situations.sort(reverse=False)
    
    with h5py.File(
            os.path.join(csts.DATA_ROOT, 
                         f'data_{n_instances}_instances.hdf5'), 'w') as f:
        for n_instance in range(0, n_instances):
            dico_players_insti_tperiods = dict()
            dico_players_insti_tperiods = generate_one_instance_tperiods(
                                            data=data,
                                            situations=situations, 
                                            n_instance=n_instance,
                                            t_periods=t_periods,
                                            m_A=m_A, m_B=m_B, m_C=m_C)
            json_players_insti_tperiods = json.dumps(dico_players_insti_tperiods, 
                                                     indent=4)
            # save to json or h5
            with open(os.path.join(
                        csts.DATA_ROOT,
                        f'data_instance_{n_instance}_tperiods_{t_periods}.json'), 
                    'w') as fp:
                json.dump(json_players_insti_tperiods, fp)
                
            # insert to h5py, hdf5
            dset = f.create_dataset(
                    f'data_instance_{n_instance}_tperiods_{t_periods}.json', 
                    data=json_players_insti_tperiods
                    )
    
#______________________________________________________________________________
#           generate n instanaces with t_periods periods: fin 
#______________________________________________________________________________

    
if __name__ == "__main__":
    json_file = "data1_generation.json"
    m_A=10; m_B=8; m_C=7
    
    ti = time.time()
    # n_instances = 10
    n_instances = 10
    # generate_n_instances(file_name=json_file, 
    #                      n_instances=n_instances,
    #                      m_A=m_A, m_B=m_B, m_C=m_C)
    
    t_periods = 5
    generate_n_instances_t_periods(file_name="data1_generation.json", 
                                   n_instances=n_instances, 
                                   t_periods=t_periods,
                                   m_A=m_A, m_B=m_B, m_C=m_C, 
                                   path_2_save="")
    
    print(f"Runtime = {time.time()-ti}")

    