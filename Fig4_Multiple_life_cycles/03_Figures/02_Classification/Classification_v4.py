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



# def Illustrate_class(data = [], str_parameter = 'none', str_title = 'dummy', str_par_name = 'dummy_x_axis', log_flag = 1):
# 	plt.style.use('default')
# 	fig, ax = plt.subplots(1,1)
# 	if log_flag == 1:
# 		plt.hist(np.log10(data[str_parameter]+1e-10), bins = 100)
# 	else:
# 		plt.hist(data[str_parameter], bins = 100)
# 	plt.yscale('log')
# 	ax.set_aspect(1.0/ax.get_data_ratio())
# 	ax.set_title(str_title, fontsize = 15)
# 	ax.set_xlabel(str_par_name, fontsize = 15)
# 	ax.set_ylabel('Counts', fontsize = 15, rotation = 90)
# 	plt.show()		
# 	return 0



# Flag2SaveFigs = 0

Flag2ReadData = 1
Flag2Round = 0
Flag2SaveFigs = 1
# Flag2CleanVariables = 1

# epsilon = 1e-3
# epsilon = 1e-4
epsilon = 0
Measure_par = "SumStd"
Measure_decr = 'Cumulative variance of abundances, $\log(Std)$'
""" files locations and descriptors """




DataFolder = '../../02_Data/'
SourceFilePath = DataFolder+'Seven_dynamics_r03_Merged.txt'
FilterFolder = '../01_Filtering/Reports/'
FilterFile = FilterFolder + 'Indices2drop_SD_r03.txt'
FileLabel = 'SD_r03'


""" main script """

""" read data """

OutColumns = ['Sample', ' Replicate', ' Groups_in_LC_0', ' Groups_in_LC_1', ' Groups_in_LC_2', ' Groups_in_LC_3', ' Groups_in_LC_4', ' Groups_in_LC_5', ' Groups_in_LC_6']
ValueColumns = [' Groups_in_LC_0', ' Groups_in_LC_1', ' Groups_in_LC_2', ' Groups_in_LC_3', ' Groups_in_LC_4', ' Groups_in_LC_5', ' Groups_in_LC_6']
PresenceColumns = [' presence_LC_0', ' presence_LC_1', ' presence_LC_2', ' presence_LC_3', ' presence_LC_4', ' presence_LC_5', ' presence_LC_6']
if Flag2ReadData == 1:
	DynamicsData = pd.read_csv(SourceFilePath)
	FilterInds = np.loadtxt(FilterFile)
	DynamicsData = DynamicsData.drop(FilterInds)
	ResData = DynamicsData[OutColumns]
	mindex = pd.MultiIndex.from_frame(ResData[['Sample', ' Replicate']]) 
	ResData = pd.DataFrame(ResData[ValueColumns].values, index = mindex, columns = ValueColumns)
	if Flag2Round ==1 :
		ResData = ResData.round(decimals = 3)
	del(DynamicsData)


ResData["Components"] = 0
ResData["Components_pattern"] = 0
ResData["PopSize"] = 0
CPatternWeights = [1, 2, 4, 8, 16, 32, 64]
for i, CLM in enumerate(ValueColumns):
	ResData["Components"] += (ResData[CLM]>epsilon)
	ResData["Components_pattern"] += CPatternWeights[i] * (ResData[CLM]>epsilon)
	ResData["PopSize"] += ResData[CLM]

for x,y in zip(ValueColumns, PresenceColumns):
	ResData[y] = (ResData[x]>epsilon).astype(int)


g_Data = (ResData["Components"].groupby(level = 0).min()).to_frame(name = 'Min_components')
g_Data['Max_components'] = ResData["Components"].groupby(level = 0).max()
for lc in PresenceColumns:
	g_Data[lc+'_all'] = ResData[lc].groupby(level = 0).min()
for lc in PresenceColumns:
	g_Data[lc+'_any'] = ResData[lc].groupby(level = 0).max()

g_Data['Varied_pattern'] = (ResData["Components_pattern"].groupby(level = 0).std()>0).astype(int)

g_Data['Universal_components'] = 0
g_Data['Occasional_components'] = 0
for lc in PresenceColumns:
	g_Data['Universal_components'] += g_Data[lc+'_all']
	g_Data['Occasional_components'] += g_Data[lc+'_any']


""" Classification step 1: dominance, repeating coexistence, bi-stability between single LCs, and others """
Inds_all = g_Data.index

Cond_dom = g_Data['Occasional_components'] == 1
Inds_single_dominance = g_Data[Cond_dom].index
# Inds_non_sd = g_Data[g_Data['Occasional_components'] != 1].index

Cond_repeat_coex = (g_Data['Varied_pattern']==0)&(g_Data['Universal_components'] > 1)
Inds_coex = g_Data[Cond_repeat_coex].index

Cond_simple_bistab = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] > 1)
Inds_bistab = g_Data[Cond_simple_bistab].index

Set_others = set(Inds_all) - set(Inds_single_dominance) - set(Inds_bistab) - set(Inds_coex)
Inds_others = pd.Int64Index(np.sort(list(Set_others)))

g_Data['Flag_others'] = 0
g_Data['Flag_others'][Inds_others] = 1

""" Classification step 2: finer breakdown of simple coexistence and bi-stability """

Cond_coex_2 = (g_Data['Varied_pattern']==0)&(g_Data['Occasional_components'] == 2)
Inds_coex_2 = g_Data[Cond_coex_2].index

Cond_coex_3 = (g_Data['Varied_pattern']==0)&(g_Data['Occasional_components'] == 3)
Inds_coex_3 = g_Data[Cond_coex_3].index

Cond_coex_4 = (g_Data['Varied_pattern']==0)&(g_Data['Occasional_components'] == 4)
Inds_coex_4 = g_Data[Cond_coex_4].index

Cond_coex_5 = (g_Data['Varied_pattern']==0)&(g_Data['Occasional_components'] == 5)
Inds_coex_5 = g_Data[Cond_coex_5].index

Cond_coex_6 = (g_Data['Varied_pattern']==0)&(g_Data['Occasional_components'] == 6)
Inds_coex_6 = g_Data[Cond_coex_6].index

Cond_coex_7 = (g_Data['Varied_pattern']==0)&(g_Data['Occasional_components'] == 7)
Inds_coex_7 = g_Data[Cond_coex_7].index

Cond_bistab_2 = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] == 2)
Inds_bistab_2 = g_Data[Cond_bistab_2].index

Cond_bistab_3 = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] == 3)
Inds_bistab_3 = g_Data[Cond_bistab_3].index

Cond_bistab_4 = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] == 4)
Inds_bistab_4 = g_Data[Cond_bistab_4].index

Cond_bistab_5 = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] == 5)
Inds_bistab_5 = g_Data[Cond_bistab_5].index

Cond_bistab_6 = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] == 6)
Inds_bistab_6 = g_Data[Cond_bistab_6].index

Cond_bistab_7 = (g_Data['Max_components']==1)&(g_Data['Occasional_components'] == 7)
Inds_bistab_7 = g_Data[Cond_bistab_7].index


""" Classification step 3: investigtion into others class """

Cond_others_2 = (g_Data['Flag_others']==1)&(g_Data['Occasional_components'] == 2)
Inds_others_2 = g_Data[Cond_others_2].index

Cond_others_3 = (g_Data['Flag_others']==1)&(g_Data['Occasional_components'] == 3)
Inds_others_3 = g_Data[Cond_others_3].index

Cond_others_4 = (g_Data['Flag_others']==1)&(g_Data['Occasional_components'] == 4)
Inds_others_4 = g_Data[Cond_others_4].index

Cond_others_5 = (g_Data['Flag_others']==1)&(g_Data['Occasional_components'] == 5)
Inds_others_5 = g_Data[Cond_others_5].index

Cond_others_6 = (g_Data['Flag_others']==1)&(g_Data['Occasional_components'] == 6)
Inds_others_6 = g_Data[Cond_others_6].index

Cond_others_7 = (g_Data['Flag_others']==1)&(g_Data['Occasional_components'] == 7)
Inds_others_7 = g_Data[Cond_others_7].index

Cond_bistab_21 = (g_Data['Occasional_components'] == 3)&(g_Data['Universal_components'] == 0)&(g_Data['Max_components']==2)&(g_Data['Min_components']==1)
Inds_bistab_21 = g_Data[Cond_bistab_21].index

Cond_bistab_22 = (g_Data['Occasional_components'] == 3)&(g_Data['Universal_components'] == 1)&(g_Data['Max_components']==2)&(g_Data['Min_components']==2)
Inds_bistab_22 = g_Data[Cond_bistab_22].index

Cond_coex_loss = (g_Data['Occasional_components'] == 3)&(g_Data['Universal_components'] == 2)&(g_Data['Max_components']==3)&(g_Data['Min_components']==2)
Inds_coex_loss = g_Data[Cond_coex_loss].index

""" plot results """

Num_only_LC = len(Inds_single_dominance)

Num_bistab_2 = len(Inds_bistab_2)
Num_bistab_3 = len(Inds_bistab_3)
Num_bistab_4 = len(Inds_bistab_4)
Num_bistab_5 = len(Inds_bistab_5)
Num_bistab_6 = len(Inds_bistab_6)
Num_bistab_7 = len(Inds_bistab_7)

Num_bistab_4p = Num_bistab_4 + Num_bistab_5 + Num_bistab_6 + Num_bistab_7

Num_coex_2 = len(Inds_coex_2)
Num_coex_3 = len(Inds_coex_3)
Num_coex_4 = len(Inds_coex_4)
Num_coex_5 = len(Inds_coex_5)
Num_coex_6 = len(Inds_coex_6)
Num_coex_7 = len(Inds_coex_7)

Num_coex_4p = Num_coex_4 + Num_coex_5 + Num_coex_6 + Num_coex_7

Num_others = len(Inds_others)

N_bistab = Num_bistab_2 + Num_bistab_3 + Num_bistab_4p
N_coex = Num_coex_2 + Num_coex_3 + Num_coex_4p



# CLR_dom = '#fb9a99'
# CLR_coex = '#a6cee3'
# CLR_bist = '#b2df8a'
# CLR_others = '#cccccc'
CLR_dom = '#F8E8B3'
CLR_coex = '#FA9A92'
CLR_bist = '#5C74B9'
CLR_others = '#BCC0C8'

# CLR_set = [CLR_coex, CLR_dom, CLR_bist, CLR_others]
# Counts = [N_coex, Num_only_LC, N_bistab, Num_others]

CLR_set = [CLR_others, CLR_bist, CLR_coex, CLR_dom]
Counts = [Num_others,N_bistab, N_coex, Num_only_LC]

size = 0.4

plt.style.use('default')
fig, ax = plt.subplots(1,1)
ax.pie(Counts, colors = CLR_set, startangle=270, wedgeprops=dict(width=size, edgecolor='w'))
if Flag2SaveFigs == 1:
	plt.savefig('Figs/Piechart_simple_'+FileLabel+'_v2.png', dpi=300, bbox_inches='tight')
plt.show()	



CLR_dom = '#fb9a99'
CLR_coex_2 = '#a6cee3'
CLR_coex_3 = '#63A3CC'
CLR_coex_4p = '#1F78B4'
CLR_bist_2 = '#b2df8a'
CLR_bist_3 = '#73C05B'
CLR_bist_4p = '#33A02C'
CLR_others = '#cccccc'
CLR_set = [CLR_coex_4p, CLR_coex_3, CLR_coex_2, CLR_dom, CLR_bist_2, CLR_bist_3, CLR_bist_4p, CLR_others]

Counts = [Num_coex_4p, Num_coex_3, Num_coex_2, Num_only_LC, Num_bistab_2, Num_bistab_3, Num_bistab_4p, Num_others]

size = 0.4

plt.style.use('default')
fig, ax = plt.subplots(1,1)
ax.pie(Counts, colors = CLR_set, startangle=270, wedgeprops=dict(width=size, edgecolor='w'))
if Flag2SaveFigs == 1:
	plt.savefig('Figs/Piechart_'+FileLabel+'.png', dpi=300, bbox_inches='tight')
plt.show()	



# # if Flag2SaveFigs == 1:
# # 	np.savetxt('Reports/Indices2drop_' + FileLabel + '.txt', Indices2Drop, delimiter=",")
# # 	np.savetxt('Reports/Indices2keep_' + FileLabel + '.txt', Indices2Keep, delimiter=",")

