#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:32:36 2020

@author: pichugin
"""

#import LifeCycleSupplementary_v1_3_3 as  LCS_old
import LifeCycleSupplementary_v2_1_1_m as LCS
import numpy as np
import sys


FullLCset = LCS.PartList(19)
LCset = [FullLCset[1]]
MaxSize = LCS.MinMaxSizes(LCset[-1])[1]
# JobID = int(sys.argv[1]) 
# SampleSize = 200
Replicates = 100

# JobID = 0
# SampleSize = 10
# Replicates = 2

T_record = np.arange(0, 10, 0.1)

file = open('Results/Single_Dynamics_21.txt', 'w')

headerline = ''
headerline += 'Run, Group_size, '
for i in np.arange(MaxSize):
	headerline += 'b_'+str(i+1)+', '
for i in np.arange(MaxSize):
	headerline += 'd_'+str(i+1)+', '
for i in np.arange(MaxSize):
	for j in np.arange(MaxSize):
		headerline += 'K_'+str(i+1)+'_'+str(j+1)+', '
for i in np.arange(len(T_record)):
	headerline += 'T_' + str(i) + ', '
headerline += 'dummy' + '\n'

file.write(headerline)

X_line = '-1, -1, '
for i in np.arange(MaxSize):
	X_line += '0, '
for i in np.arange(MaxSize):
	X_line += '-1, '
for i in np.arange(MaxSize):
	for j in np.arange(MaxSize):
		X_line += '0, '
for i in np.arange(len(T_record)):
	X_line += str(T_record[i]) + ', '
X_line += '-1' + '\n'
file.write(X_line)
	
""" generate control parameters (event rates) """
b = LCS.BirthInit('direct', MaxSize, [1,1,1])
d = LCS.DeathInit('', MaxSize, [])
K = LCS.InteractionInit('direct', MaxSize, [[1,0.2],[0.2,0.5]])


ProjectionMatricesSet = []
for LC in LCset:
	ProjMatrix = LCS.ProjectionMatrix(LC, b, d, MaxSize)
	ProjectionMatricesSet.append(ProjMatrix)

for r_count in np.arange(Replicates):
	
	""" initialize population """
	LCset_init_pop = []
	LCset_init_weights = np.random.uniform(low = 0.1, high = 2.0, size = len(LCset))
	for i, LC in enumerate(LCset):
		LC_pop = np.zeros(MaxSize)
		LC_pop[0] = np.random.uniform(low = 0.1, high = 2.0, size = len(LCset))
		LC_pop[1] = np.random.uniform(low = 0.1, high = 2.0, size = len(LCset))
# 		LC_max_size = LC[0][0] - 1
# 		if LC_max_size == 1:
# 			LC_pop[0] = LCset_init_weights[i]
# 		elif LC_max_size == 2:
# 			LocalWeights = np.random.uniform(low = 0.1, high = 1.0, size = 2)
# 			LocalWeights = LocalWeights/np.sum(LocalWeights)
# 			LC_pop[0] = LCset_init_weights[i] * LocalWeights[0]
# 			LC_pop[1] = LCset_init_weights[i] * LocalWeights[1]
# 		else:
# 			LocalWeights = np.random.uniform(low = 0.1, high = 1.0, size = 3)
# 			LocalWeights = LocalWeights/np.sum(LocalWeights)
# 			LC_pop[0] = LCset_init_weights[i] * LocalWeights[0]
# 			LC_pop[1] = LCset_init_weights[i] * LocalWeights[1]
# 			LC_pop[2] = LCset_init_weights[i] * LocalWeights[2]
		LCset_init_pop.append(LC_pop)
	
	""" compute dynamics """
	Output = LCS.LifeCyclesCompetition(ProjectionMatricesSet, K, LCset_init_pop, T_record)
	
	""" record output """
	for group_size in np.arange(MaxSize):
		line2write = ''
		line2write += str(r_count) + ', '
		line2write += str(group_size+1) + ', '
		for i in np.arange(MaxSize):
			line2write += str(b[i])+', '
		for i in np.arange(MaxSize):
			line2write += str(d[i])+', '
		for i in np.arange(MaxSize):
			for j in np.arange(MaxSize):
				line2write += str(K[i,j])+', '
		for i in np.arange(len(T_record)):
			line2write += str(Output[0, i, group_size]) + ', '
		line2write += str(-1)+'\n'	
		file.write(line2write)
		file.flush()

file.close()