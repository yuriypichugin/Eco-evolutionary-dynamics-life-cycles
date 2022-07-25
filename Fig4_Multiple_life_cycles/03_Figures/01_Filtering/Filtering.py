#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 09:28:08 2021

@author: pichugin


Filters out the entries, which has less than X repeats completed

"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from copy import deepcopy as deepcopy
from scipy import interpolate

Flag2SaveFigs = 1
Flag2ReadData = 1

""" files locations and descriptors """



DataFolder = '../../02_Data/'
SourceFilePath = DataFolder+'Seven_dynamics_r03_Merged.txt'
FileLabel = 'SD_r03'
ReplicateNum = 100

""" main script """
if Flag2ReadData == 1:
	DynamicsData = pd.read_csv(SourceFilePath)
DataSubset = DynamicsData[['Sample', ' Replicate']]
SampleSet = np.unique(DataSubset['Sample'])
df_Grouped = DataSubset.groupby(['Sample']).size()
Samples2Drop = df_Grouped[df_Grouped < ReplicateNum].index.values
Indices2Drop = [i for i in DataSubset.index if DataSubset['Sample'][i] in Samples2Drop]
print('tygydym')
# Indices2Keep = [x for x in DataSubset.index.values if x not in Indices2Drop]
# Indices2Keep = [x for x in DataSubset.index if x not in Indices2Drop]


if Flag2SaveFigs == 1:
	np.savetxt('Reports/Indices2drop_' + FileLabel + '.txt', Indices2Drop, delimiter=",")
 	# np.savetxt('Reports/Indices2keep_' + FileLabel + '.txt', Indices2Keep, delimiter=",")

