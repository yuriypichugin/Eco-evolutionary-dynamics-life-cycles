#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 08:55:10 2020

@author: pichugin
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from copy import deepcopy as deepcopy
import os

for FileIndex in np.arange(0, 100):
	
    
	File2ReadName = 'Raw data/Random_screen_invasion_r17_'+str(FileIndex)+'.txt'
	if not(os.path.exists(File2ReadName)):
		print(FileIndex)
		continue

#	print(FileIndex)
	Data = pd.read_csv(File2ReadName)
	
	if FileIndex == 0:
		OutputData = deepcopy(Data)
	else:
		OutputData = pd.concat([OutputData, Data])

#OutputData = OutputData.drop(columns = 'sim_id')

File2WriteName = 'Pairwise_Invasion_calc_17_Merged.txt'
OutputData.to_csv(File2WriteName, index=True, header=True, float_format = "%.5f")