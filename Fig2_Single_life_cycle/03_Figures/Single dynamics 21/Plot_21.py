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

Path = '../../02_Data/Single_Dynamics_21.txt'

Data = pd.read_csv(Path)

df_X = Data[Data[' Group_size'] == -1].values
df_C1 = Data[Data[' Group_size'] == 1].values
df_C2 = Data[Data[' Group_size'] == 2].values

X = df_X[0, 10:110]
C1 = df_C1[:, 10:110]
C2 = df_C2[:, 10:110]

plt.style.use('default')
fig, ax = plt.subplots(1,1)
for i in np.arange(C1.shape[0]):
	plt.plot(X, C1[i, :], color = '#e41a1c')
for i in np.arange(C2.shape[0]):
	plt.plot(X, C2[i, :], color = '#4daf4a')
	
ax.set_ylim(0, 2)
ax.set_xlim(0, 10)

# ax.set_ylim(0, 20)
# ax.set_xlim(0, 20)

ax.set_aspect(1.0/ax.get_data_ratio())

ax.set_xlabel('Time, '+ r'$t$', fontsize = 15)
ax.set_ylabel('Number of groups, '+ r'$x_i$', fontsize = 15)
ax.set_title('multicellular life cycle 2+1', fontsize = 15)
# ax.set_xticks([0.5, 5.5, 10.5, 15.5])
# ax.set_xticklabels([r'$10^{-2.5}$', r'$10^{-2}$', r'$10^{-1.5}$', r'$10^{-1}$'], rotation = 'horizontal')
# ax.set_yticks([0.5, 5.5, 10.5, 15.5])
# ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$'], rotation = 'horizontal')

if Flag2SaveFigs == 1:
 	plt.savefig('Figs/Single_life_cycle_21.png', dpi=300, bbox_inches='tight')
plt.show()




plt.style.use('default')
fig, ax = plt.subplots(1,1)
for i in np.arange(C1.shape[0]):
	plt.plot(C1[i, :], C2[i, :], color = '#3182bd', zorder = 0)
for i in np.arange(C1.shape[0]):
	plt.scatter(C1[i, -1], C2[i, -1], color = '#ff7f00', s = 100, zorder = 1)

	
ax.set_ylim(0, 2)
ax.set_xlim(0, 2)

# ax.set_ylim(0, 20)
# ax.set_xlim(0, 20)

ax.set_aspect(1.0/ax.get_data_ratio())

ax.set_xlabel('Number of solitary cells, '+ r'$x_1$', fontsize = 15)
ax.set_ylabel('Number of bi-cells, '+ r'$x_2$', fontsize = 15)
ax.set_title('multicellular life cycle 2+1', fontsize = 15)
# ax.set_xticks([0.5, 5.5, 10.5, 15.5])
# ax.set_xticklabels([r'$10^{-2.5}$', r'$10^{-2}$', r'$10^{-1.5}$', r'$10^{-1}$'], rotation = 'horizontal')
# ax.set_yticks([0.5, 5.5, 10.5, 15.5])
# ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$'], rotation = 'horizontal')

if Flag2SaveFigs == 1:
 	plt.savefig('Figs/Single_life_cycle_21_traj.png', dpi=300, bbox_inches='tight')
plt.show()




