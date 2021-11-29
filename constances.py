#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:24:00 2021

@author: willy
"""

DATA_ROOT = "data"

SITUATION_ERR = "Situation_Error"

PLAYER_ROOT = "player_"

PERIOD_ROOT = "t_"
PROB_NAME = "p_"

STATES = ["Deficit", "Self", "Surplus"]

STATE1_STRATS = ("CONS+", "CONS-")                                             # strategies possibles pour l'etat 1 de a_i
STATE2_STRATS = ("DIS", "CONS-")                                               # strategies possibles pour l'etat 2 de a_i
STATE3_STRATS = ("DIS", "PROD")                                                # strategies possibles pour l'etat 3 de a_i

"""
bgr_is = BG_r_i
"""
LEARNING_PROPERTIES = ["Pis", "Sis", "mode_is", "ksteps", "p1_ijks", "p2_ijks",
        "state_i", "Si_minus", "Si_plus", "ppi", "gamma_i", "q_minus", "q_plus", 
        "prod_is", "cons_is", "r_is", "ben_is", "cst_is", "V_is", "u_is",
        "bg_is", "bgr_is"]