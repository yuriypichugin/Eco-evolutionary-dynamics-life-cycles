#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:32:36 2020

@author: pichugin
"""

import LifeCycleSupplementary_v2_1_1 as LCS
# import matplotlib.pyplot as plt
import pandas as pd

import numpy as np




""" Constants initialization """
Flag_SaveFigs = 0


""" Competing life cycles set up """
LCset = LCS.PartList(19)

LC_1 = LCset[0]
LC_2 = LCset[1]

MaxSize_1 = LCS.MinMaxSizes(LC_1)[1]
MaxSize_2 = LCS.MinMaxSizes(LC_2)[1]
MaxSize = np.max([MaxSize_1, MaxSize_2])

""" environment set up """
T_record = np.linspace(0, 300, 3000)
b_vals = [1.0, 0.5]
b = LCS.BirthInit('direct', MaxSize, b_vals)
d = LCS.DeathInit('', MaxSize, [])

PM_1 = LCS.ProjectionMatrix(LC_1, b, d, MaxSize)
PM_2 = LCS.ProjectionMatrix(LC_2, b, d, MaxSize)
PM_Composite = [PM_1, PM_2]


# """ Dominance 1 """
# Ksetup = [[1, 1],[1, 1]]
# FileMask = 'Dominance_1'

# """ Dominance 2 """
# Ksetup = [[3, 3],[1, 1]]
# FileMask = 'Dominance_2'

""" bi-stability """
Ksetup = [[1, 1],[0.6, 0.4]]
FileMask = 'Bi_stability'

# """ coexistence """
# Ksetup = [[1, 0.1],[0.1, 1.0]]
# FileMask = 'Coexistence'


K = LCS.InteractionInit('direct', MaxSize, Ksetup)

Init_pop_set = np.arange(0.1, 1.1, 0.1)

for count_1, start_1 in enumerate(Init_pop_set):
	print(start_1)
	for count_2, start_2 in enumerate(Init_pop_set):
		InitState_1 = np.zeros(MaxSize)
		InitState_1[0] = start_1
		InitState_2 = np.zeros(MaxSize)
		InitState_2[0] = start_2
		InitComposite = [InitState_1, InitState_2]
		Record = LCS.LifeCyclesCompetition(PM_Composite, K, InitComposite, T_record)
		Report_LC_1 = np.sum(Record[0,:,:], axis = 1)
		Report_LC_2 = np.sum(Record[1,:,:], axis = 1)
		df_Report = pd.DataFrame(np.transpose([T_record, Report_LC_1, Report_LC_2]), columns = ['Time', 'LC_1_1', 'LC_2_1'])
		df_Report.to_csv('../02 Data/Bistability/'+FileMask+'_report_'+str(count_1)+'_'+str(count_2)+'.txt')