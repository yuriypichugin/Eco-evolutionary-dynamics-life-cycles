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
LCset = [FullLCset[0], FullLCset[1], FullLCset[2], FullLCset[3], FullLCset[4], FullLCset[5], FullLCset[6]]
MaxSize = LCS.MinMaxSizes(LCset[-1])[1]

JobID = int(sys.argv[1]) 
SampleSize = 200
Replicates = 100

# JobID = 0
# SampleSize = 10
# Replicates = 2

T_end = 5000

file = open('Results/Seven_Dynamics_r03_'+str(JobID)+'.txt', 'w')

headerline = ''
headerline += 'Sample, Replicate, '
for i in np.arange(MaxSize):
	headerline += 'b_'+str(i+1)+', '
for i in np.arange(MaxSize):
	headerline += 'd_'+str(i+1)+', '
for i in np.arange(MaxSize):
	for j in np.arange(MaxSize):
		headerline += 'K_'+str(i+1)+'_'+str(j+1)+', '
for i in np.arange(len(LCset)):
	headerline += 'Groups_in_LC_' + str(i) + ', '
headerline += 'dummy' + '\n'

file.write(headerline)


for s_count in np.arange(SampleSize):
	
	""" generate control parameters (event rates) """
	b = LCS.BirthInit('rand_exp', MaxSize, [1])
	d = LCS.DeathInit('', MaxSize, [])
	K = LCS.InteractionInit('rand_exp', MaxSize, [1])


	ProjectionMatricesSet = []
	for LC in LCset:
		ProjMatrix = LCS.ProjectionMatrix(LC, b, d, MaxSize)
		ProjectionMatricesSet.append(ProjMatrix)
	
	for r_count in np.arange(Replicates):
		
		""" initialize population """
		LCset_init_pop = []
		LCset_init_weights = np.random.uniform(low = 0.1, high = 1.0, size = len(LCset))
		for i, LC in enumerate(LCset):
			LC_pop = np.zeros(MaxSize)
			LC_max_size = LC[0][0] - 1
			if LC_max_size == 1:
				LC_pop[0] = LCset_init_weights[i]
			elif LC_max_size == 2:
				LocalWeights = np.random.uniform(low = 0.1, high = 1.0, size = 2)
				LocalWeights = LocalWeights/np.sum(LocalWeights)
				LC_pop[0] = LCset_init_weights[i] * LocalWeights[0]
				LC_pop[1] = LCset_init_weights[i] * LocalWeights[1]
			else:
				LocalWeights = np.random.uniform(low = 0.1, high = 1.0, size = 3)
				LocalWeights = LocalWeights/np.sum(LocalWeights)
				LC_pop[0] = LCset_init_weights[i] * LocalWeights[0]
				LC_pop[1] = LCset_init_weights[i] * LocalWeights[1]
				LC_pop[2] = LCset_init_weights[i] * LocalWeights[2]
			LCset_init_pop.append(LC_pop)
		
		""" compute dynamics """
		Output = LCS.LifeCyclesCompetition(ProjectionMatricesSet, K, LCset_init_pop, [T_end])
		GroupCount = np.sum(np.sum(Output, axis = 1), axis = 1)
		
		""" record output """
		line2write = ''
		line2write += str(JobID * SampleSize + s_count) + ', '
		line2write += str(r_count) + ', '
		for i in np.arange(MaxSize):
			line2write += str(b[i])+', '
		for i in np.arange(MaxSize):
			line2write += str(d[i])+', '
		for i in np.arange(MaxSize):
			for j in np.arange(MaxSize):
				line2write += str(K[i,j])+', '
		for entry in GroupCount:
			record_entry = entry
			if entry < 1e-5:
				record_entry = 0
			line2write += str(record_entry)+', '
		line2write += str(-1)+'\n'
		file.write(line2write)
		file.flush()
		
# 		print(np.sum(LCset_init_pop), np.sum(GroupCount))

file.close()