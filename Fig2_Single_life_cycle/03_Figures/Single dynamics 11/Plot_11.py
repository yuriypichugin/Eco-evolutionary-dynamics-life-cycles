#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 13:26:50 2022

@author: pichugin
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Flag2SaveFigs = 1

Path = '../../02_Data/Single_Dynamics_11.txt'
Data = pd.read_csv(Path).values

X = Data[0, 5:105]
Data = Data[1:, 5:105]

plt.style.use('default')
fig, ax = plt.subplots(1,1)
for i in np.arange(Data.shape[0]):
	plt.plot(X, Data[i, :], color = '#e41a1c')
	
ax.set_ylim(0, 2)
ax.set_xlim(0, 10)

# ax.set_ylim(0, 20)
# ax.set_xlim(0, 20)

ax.set_aspect(1.0/ax.get_data_ratio())

ax.set_xlabel('Time, '+ r'$t$', fontsize = 15)
ax.set_ylabel('Number of groups, '+ r'$x_i$', fontsize = 15)
ax.set_title('Unicellular life cycle 1+1', fontsize = 15)
# ax.set_xticks([0.5, 5.5, 10.5, 15.5])
# ax.set_xticklabels([r'$10^{-2.5}$', r'$10^{-2}$', r'$10^{-1.5}$', r'$10^{-1}$'], rotation = 'horizontal')
# ax.set_yticks([0.5, 5.5, 10.5, 15.5])
# ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$'], rotation = 'horizontal')

if Flag2SaveFigs == 1:
 	plt.savefig('Figs/Single_life_cycle_11.png', dpi=300, bbox_inches='tight')
plt.show()
