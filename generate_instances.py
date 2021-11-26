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
import pandas as pd

import constances as csts

def test(name_file="data1.json"):
    f = json.open(os.path.join(csts.data_root, name_file))