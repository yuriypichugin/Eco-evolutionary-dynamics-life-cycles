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

LC_1 = LCset[4]
LC_2 = LCset[43]
LC_3 = LCset[45]

MaxSize_1 = LCS.MinMaxSizes(LC_1)[1]
MaxSize_2 = LCS.MinMaxSizes(LC_2)[1]
MaxSize_3 = LCS.MinMaxSizes(LC_3)[1]
MaxSize = np.max([MaxSize_1, MaxSize_2, MaxSize_3])

""" environment set up """
T_record = np.linspace(0, 200, 3000)
b_vals = MaxSize * [1.0]
b = LCS.BirthInit('direct', MaxSize, b_vals)
d = LCS.DeathInit('', MaxSize, [])

PM_1 = LCS.ProjectionMatrix(LC_1, b, d, MaxSize)
PM_2 = LCS.ProjectionMatrix(LC_2, b, d, MaxSize)
PM_3 = LCS.ProjectionMatrix(LC_3, b, d, MaxSize)
PM_Composite = [PM_1, PM_2, PM_3]



# """ bi-stability """
# Ksetup = 0.1*LCS.InteractionInit('none', MaxSize, [])
# Ksetup[1, -1] = 10
# Ksetup[-1, 1] = 10
# PopScale = 1
# FileMask = 'Bistability_4_43_45'


# """ Coexistence """
# Ksetup = 0.1*LCS.InteractionInit('none', MaxSize, [])
# dd=6
# for i in np.arange(dd+1):
#  	Ksetup[i, dd-i] = 0.0
# PopScale = 1
# FileMask = 'Coexistence_4_43_45'
# Step = 0.1

""" Hierarchical dominance """
Ksetup = 0.1*LCS.InteractionInit('none', MaxSize, [])
dd=6
for i in np.arange(dd+1):
 	Ksetup[-1, i] = 1.1
PopScale = 10.0
FileMask = 'HierDom_4_43_45'


# """ Non-hierarchical dominance """
# """ read data from file """
# FileName_anom = '../../../Data/SimulationData/DataAnalysis/Anomalous_search/AnomalousData_calc_25.txt'
# Data_anom = pd.read_csv(FileName_anom)
# Data_anom = Data_anom.drop(columns = ['Unnamed: 0', 'level_0', 'index'])
# Entry_number = 2
# Param_set = Data_anom.iloc[Entry_number,:]
# b = LCS.BirthInit('direct', MaxSize, np.asarray(Param_set[0:MaxSize]))
# d = LCS.DeathInit('direct', MaxSize, np.asarray(Param_set[MaxSize:2*MaxSize]))
# Ksetup = LCS.InteractionInit('const', MaxSize, [0])
# index = 2*MaxSize
# for i in np.arange(MaxSize):
# 	for j in np.arange(MaxSize):
# 		Ksetup[i,j] = Param_set[index]
# 		index += 1
# PopScale = 10.0
# FileMask = 'NonhierDom_4_43_45_high_step'
# T_record = np.linspace(0, 200, 3000)


K = LCS.InteractionInit('direct', MaxSize, Ksetup)


""" set up initial states """
XR_1 = LCS.StationaryState(LC_1, b, d, K, MaxSize)
XR_2 = LCS.StationaryState(LC_1, b, d, K, MaxSize)
XR_3 = LCS.StationaryState(LC_1, b, d, K, MaxSize)


Step = 0.1

Init_pop_set = np.arange(Step, 1.00001 - 2*Step, Step)

for count_1, start_1 in enumerate(Init_pop_set):
	print(start_1)
	SecondPopSet = np.arange(Step, 1.00001 - start_1 - Step, Step)
	for count_2, start_2 in enumerate(SecondPopSet):
		start_3 = 1.0 - start_1 - start_2
		
		InitState_1 = start_1 * XR_1
		InitState_2 = start_2 * XR_2
		InitState_3 = start_3 * XR_3
		
		InitComposite = [InitState_1, InitState_2, InitState_3]
		Record = LCS.LifeCyclesCompetition(PM_Composite, K, InitComposite, T_record)
		Report_LC_1 = np.sum(Record[0,:,:], axis = 1)
		Report_LC_2 = np.sum(Record[1,:,:], axis = 1)
		Report_LC_3 = np.sum(Record[2,:,:], axis = 1)
		
		Norm = Report_LC_1 + Report_LC_2 + Report_LC_3
		Frac_1 = Report_LC_1 / Norm
		Frac_2 = Report_LC_2 / Norm
		Frac_3 = Report_LC_3 / Norm
		
		
		df_Report = pd.DataFrame(np.transpose([T_record, Frac_1, Frac_2, Frac_3]), columns = ['Time', 'LC_1_1', 'LC_2_1', 'LC_1_1_1'])
		df_Report.to_csv('../02 Data/HierDom_to_plot/'+FileMask+'_report_'+str(count_1)+'_'+str(count_2)+'.txt')
		
		
		