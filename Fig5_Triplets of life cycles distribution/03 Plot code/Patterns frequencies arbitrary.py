#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:53:50 2020

@author: pichugin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

# FileName ='../../OriginalData/Pairwise invasion screen calc 25 limitative LC_4_43_45/Pairwise_Invasion_calc_25_Merged.txt'
# FileOutName = 'Figs/Patterns_calc_25.png'


FileName ='../02 Data/Pairwise_Invasion_calc_17_Merged.txt'
FileOutName = '../04 Source figs/Patterns_calc_17.png'


#flagrecordData = 0
flag_readdata = 1
Flag2SaveFigs = 1

if flag_readdata == 1:
	Data = pd.read_csv(FileName)


""" filter condition number """
Limit = 0
# zScore = Data[' z_score']
# Data = Data[zScore > Limit]
Pattern = Data[' pattern']

"""
Pattern string, in the format Resident_X Invader_Y:
	R1I2; R1I3; R2I1; R2I3; R3I1; R3I2
"""



Totalcount = 0

UniquePatterns = np.unique(Pattern)
CountArray = []
for entry in UniquePatterns:
	count = sum(Pattern == entry)
	print('Pattern ', entry, ' is found ', count, ' times')
	CountArray.append(count)
	Totalcount += count

print(Totalcount)

plt.plot(np.sort(CountArray), 'o')
plt.yscale('log')
plt.show()

ESS_num = []
for entry in UniquePatterns:
	string_0 = entry[1:]
	list_0 = string_0.split(';')
	list_1 = [ int(x) for x in list_0 ]
	
	ess_0 = (list_1[0]+list_1[1])==0
	ess_1 = (list_1[2]+list_1[3])==0
	ess_2 = (list_1[4]+list_1[5])==0
	
	ESS_num.append(ess_0+ess_1+ess_2)
	
Collection = np.transpose([CountArray, ESS_num])
PlotDataSet = pd.DataFrame(Collection, columns = ['count', 'ESS_num'])
PlotDataSet = PlotDataSet.sort_values('count')
PlotDataSet['id'] = np.arange(len(CountArray))

CLR = ['#66c2a5', '#fc8d62', '#8da0cb', '#ffd92f']


ColorSet = []
for i in np.arange(len(ESS_num)):
	ColorSet.append(CLR[PlotDataSet['ESS_num'].iloc[i]])
	
plt.style.use('default')
fig, ax = plt.subplots(1,1)	
plt.bar(PlotDataSet['id'], PlotDataSet['count'], color = ColorSet)
plt.yscale('log')
plt.ylim((5e-1, 1e5))

#plt.title('Median = '+ str(np.median(PreExcess)))
plt.xticks([],[])
plt.xlabel('Triplet pairwise invasion patterns', fontsize = 15)
plt.ylabel('Counts', fontsize = 15)


NoESS_patch = mpatches.Patch(color=CLR[0], label='0 ESS')
OneESS_patch = mpatches.Patch(color=CLR[1], label='1 ESS')
TwoESS_patch = mpatches.Patch(color=CLR[2], label='2 ESS')
ThreeESS_patch = mpatches.Patch(color=CLR[3], label='3 ESS')
# plt.legend(handles = [NoESS_patch, OneESS_patch, TwoESS_patch, ThreeESS_patch], fontsize = 15)
plt.legend(handles = [NoESS_patch, OneESS_patch, TwoESS_patch], fontsize = 15)

ax.set_aspect(1.0/ax.get_data_ratio())
if Flag2SaveFigs == 1:
	plt.savefig(FileOutName, dpi=300, bbox_inches='tight')
plt.show()