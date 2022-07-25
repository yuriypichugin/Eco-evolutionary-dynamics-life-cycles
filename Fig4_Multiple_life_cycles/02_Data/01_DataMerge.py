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


EntriesNum = 0

for FileIndex in np.arange(0, 100):
	
    
	File2ReadName = 'Raw data/Seven_Dynamics_r03_'+str(FileIndex)+'.txt'
	if not(os.path.exists(File2ReadName)):
		print(FileIndex)
		continue

#	print(FileIndex)
	Data = pd.read_csv(File2ReadName)
	
	EntriesNum += Data.shape[0]
	
	if FileIndex == 0:
		OutputData = deepcopy(Data)
	else:
		OutputData = pd.concat([OutputData, Data])

#OutputData = OutputData.drop(columns = 'sim_id')

print('Total entries = ', EntriesNum)

File2WriteName = 'Seven_dynamics_r03_Merged.txt'
OutputData.to_csv(File2WriteName, index=True, header=True, float_format = "%.5f")