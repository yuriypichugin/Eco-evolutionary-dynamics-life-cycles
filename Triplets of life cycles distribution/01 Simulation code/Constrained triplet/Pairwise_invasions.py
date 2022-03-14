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
LCset = [FullLCset[4], FullLCset[43], FullLCset[45]]
MaxSize = LCS.MinMaxSizes(LCset[-1])[1]

JobID = int(sys.argv[1]) 
# JobID = 0

repeats = 10000
#repeats = 10

file = open('../../02 Data/Raw data/Random_screen_invasion_r25_'+str(JobID)+'.txt', 'w')

headerline = ''
for i in np.arange(MaxSize):
    headerline += 'b_'+str(i+1)+', '
for i in np.arange(MaxSize):
    headerline += 'd_'+str(i+1)+', '
for i in np.arange(MaxSize):
    for j in np.arange(MaxSize):
        headerline += 'K_'+str(i+1)+'_'+str(j+1)+', '
headerline += 'pattern, ' + 'relative_pattern, ' + 'limit_inv_rate, ' + 'self_inv_rate' +'\n'

file.write(headerline)

for count in np.arange(repeats):

    b = LCS.BirthInit('rand_exp', MaxSize, [1, 0.1])
    d = LCS.DeathInit('', MaxSize, [])
    K = LCS.InteractionInit('rand_exp', MaxSize, [1, 0.1])
    
    """ normalize parameter values (b,d) are lists, K is numpy array """
    NB = np.min(b)
    for i in np.arange(len(b)):
        b[i] = b[i]/NB
        d[i] = d[i]/NB
    K = K/np.max(K)
	
#    MatrixIsWellDefined = 1

    LogString = ''
    LogString2 = ''
    MinInvRate_zscore = 1e20
    MinInvRate = 1e20
    MaxSelfInvasion = 0
    for LC_id in np.arange(len(LCset)):
        X0 = LCS.StationaryState(LCset[LC_id], b, d, K, MaxSize)
        A0 = LCS.ProjectionMatrix(LCset[LC_id], b, d, MaxSize)
        SelfInvasionRate, dummy = LCS.InvasionRate(A0, K, X0)
        for LC_id_inv in np.arange(len(LCset)):
            if LC_id!=LC_id_inv:
                A = LCS.ProjectionMatrix(LCset[LC_id_inv], b, d, MaxSize)
                InvasionRate, dummy = LCS.InvasionRate(A, K, X0)
                """---test---"""
                if np.abs(SelfInvasionRate) >= MaxSelfInvasion:
                    MaxSelfInvasion = np.abs(SelfInvasionRate)
                    MinInvRate = 1e20
#                Zscore = np.abs((InvasionRate - SelfInvasionRate)) / (np.abs(SelfInvasionRate)+1e-10) 
                    if (np.abs(InvasionRate) <  MinInvRate):
                        MinInvRate = np.abs(InvasionRate) 
                        Record_limit_invasion_rate = InvasionRate
                        Record_self_invasion_rate = SelfInvasionRate
                """---test---"""
                LogString += str(int(InvasionRate>0))+';'
                LogString2 += str(int(InvasionRate>SelfInvasionRate))+';'
    LogString = LogString[:-1]
    LogString2 = LogString2[:-1]
    print(LogString, MinInvRate_zscore)

    """ not write results from ill conditioned projection matrices """
#    if MatrixIsWellDefined:
    line2write = ''
    for i in np.arange(MaxSize):
        line2write += str(b[i])+', '
    for i in np.arange(MaxSize):
        line2write += str(d[i])+', '
    for i in np.arange(MaxSize):
        for j in np.arange(MaxSize):
            line2write += str(K[i,j])+', '
    line2write += LogString+', '+ LogString2+ ', '
#    line2write += str(MinInvRate_zscore)+', '
    line2write += str(Record_limit_invasion_rate)+', '
    line2write += str(Record_self_invasion_rate)+'\n'
    file.write(line2write)
    file.flush()

file.close()