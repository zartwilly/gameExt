#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:24:00 2021

@author: willy
"""

DATA_ROOT = "data"

SITUATION_ERR = "Situation_Error"

#______________________________________________________________________________
#                       Suffixes variables: debut 
#______________________________________________________________________________

PLAYER_ROOT = "player_"

PERIOD_ROOT = "t_"
PROB_NAME = "p_"

#______________________________________________________________________________
#                       Suffixes variables: fin
#______________________________________________________________________________

#______________________________________________________________________________
#                       valeurs constantes: debut 
#______________________________________________________________________________
ARRONDI = 7

#______________________________________________________________________________
#                       valeurs constantes: fin 
#______________________________________________________________________________

#______________________________________________________________________________
#                       etats et modes possibles: debut 
#______________________________________________________________________________

STATES = ["Deficit", "Self", "Surplus"]

STATE1_STRATS = ("CONS+", "CONS-")                                             # strategies possibles pour l'etat 1 de a_i
STATE2_STRATS = ("DIS", "CONS-")                                               # strategies possibles pour l'etat 2 de a_i
STATE3_STRATS = ("DIS", "PROD")                                                # strategies possibles pour l'etat 3 de a_i

#______________________________________________________________________________
#                       etats et modes possibles: fin 
#______________________________________________________________________________

#______________________________________________________________________________
#                    variables communes a tous les joueurs: debut 
#______________________________________________________________________________

"""
bgr_is = BG_r_i
"""
LEARNING_PROPERTIES = ["mode_is","state_is", "P_is", "ksteps",
                       "S_t_plus_1_is", "S_minus_is", "S_plus_is", 
                       "pp_is", "gamma_is", "q_minus_is", "q_plus_is", 
                      "prod_is", "cons_is", "r_is", 
                      "ben_is", "cst_is", "V_is", "u_is", "bg_is", "bgr_is", 
                      "is_playing"]

#______________________________________________________________________________
#                    variables communes a tous les joueurs: fin
#______________________________________________________________________________
