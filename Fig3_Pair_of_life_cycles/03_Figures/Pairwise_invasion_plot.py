#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:32:36 2020

@author: pichugin
"""

# import LifeCycleSupplementary_v2_1_1 as LCS
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np



""" Constants initialization """
Flag2SaveFigs = 1

FileIndices = np.arange(0, 10)


""" four blocks for four panels. Comment/uncomment as necessary: """

# FileMask = '../02_Data/Dominance 1/Dominance_1_report_'
# OutputName = 'Dominance_1.jpg'
# FlagColoration = False

# FileMask = '../02_Data/Dominance 2/Dominance_2_report_'
# OutputName = 'Dominance_2.jpg'
# FlagColoration = False

# FileMask = '../02_Data/Bistability/Bi_stability_report_'
# OutputName = 'Bistability_v2.jpg'
# FlagColoration = True

FileMask = '../02_Data/Coexistence/Coexistence_report_'
OutputName = 'Coexistence.jpg'
FlagColoration = False

plt.style.use('default')
fig, ax = plt.subplots(1,1)
for i in FileIndices:
	for j in FileIndices:
		FilePath = FileMask +str(i)+'_'+str(j)+'.txt'
		df_data = pd.read_csv(FilePath, index_col = 0)
		if FlagColoration == True:
			if df_data['LC_1_1'].iloc[-1] < 0.5:
				ColorPlot = '#F6511D'
				ColorPoint = '#ff8c00'
			else:
				ColorPlot = '#2b8cbe'
				ColorPoint = '#00008b'
		else:
			ColorPlot = '#2b8cbe'
			ColorPoint = '#00008b'
		plt.plot(df_data['LC_1_1'], df_data['LC_2_1'], color = ColorPlot, zorder=1)
		plt.scatter(df_data.iloc[-1,1], df_data.iloc[-1,2], s= 100, color = ColorPoint, zorder=2)
plt.xlim([-0.05, 1.1])
plt.ylim([-0.05, 1.1])
ax.set_aspect(1.0/ax.get_data_ratio())
ax.set_xlabel('Abundance of '+r'$\kappa = 1+1$', fontsize = 15)
ax.set_ylabel('Abundance of '+r'$\kappa = 2+1$', fontsize = 15, rotation = 90)
if Flag2SaveFigs == 1:
 	plt.savefig(OutputName, dpi=300, bbox_inches='tight')
plt.show()
