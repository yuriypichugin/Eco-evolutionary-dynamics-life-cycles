#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:32:36 2020

@author: pichugin
"""

# import LifeCycleSupplementary_v2_1_1 as LCS
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import ternary



""" Constants initialization """
Flag2SaveFigs = 1

FileIndices = np.arange(0, 10)



FolderMask = '../02 Data/Bistability_to_plot/'
OutputName = 'Bistabilty_v2.jpg'

# FolderMask = '../02 Data/Coexistence_to_plot/'
# OutputName = 'Coexistence_v2.jpg'

# FolderMask = '../02 Data/HierDom_to_plot/'
# OutputName = 'HierDom_v2.jpg'

# FolderMask = '../02 Data/NonhierDom_to_plot/'
# OutputName = 'NonhierDom_v2.jpg'

scale = 1.0
figure, tax = ternary.figure(scale=scale)
tax.boundary(linewidth=1.5)
tax.gridlines(color="black", multiple=0.1)

for FileName in os.listdir(FolderMask):
	FilePath = FolderMask + FileName
	df_data = pd.read_csv(FilePath, index_col = 0)
	tax.plot(df_data[['LC_1_1', 'LC_2_1', 'LC_1_1_1']].values, linewidth=2.0, color = '#2b8cbe', zorder=1)
	tax.scatter([[df_data.iloc[-1,1], df_data.iloc[-1,2], df_data.iloc[-1,3]]], color = 'r', zorder=2, s = 100)

tax.ticks(axis='lbr', linewidth=1, multiple=1.0, tick_formats="%.1f", offset=0.04, fontsize = 20)

tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')



tax.ax.set_aspect(0.85/tax.ax.get_data_ratio())


if Flag2SaveFigs == 1:
 	plt.savefig(OutputName, dpi=300, bbox_inches='tight')
ternary.plt.show()



