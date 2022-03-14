#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:32:36 2020

@author: pichugin
"""

#import LifeCycleSupplementary_v1_3_3 as  LCS_old
import LifeCycleSupplementary_v2_1 as LCS
import numpy as np
import sys
from copy import deepcopy as deepcopy
import matplotlib.pyplot as plt


def PlotDemDyn(T, X, ylimits = [], plot_title=''):
	""" plots the dynamics of the life cycles """
	CLR = ['#e41a1c', '#4daf4a', '#377eb8']
	
	plt.style.use('default')
	fig, ax = plt.subplots(1,1)
	for size in np.arange(np.shape(X)[1]):
		plt.plot(T, X[:,size], CLR[size], linewidth = 4)
	if len(ylimits)>0:
		plt.ylim(ylimits[0], ylimits[1])
	ax.set_aspect(1.0/ax.get_data_ratio())
# 	plt.title(plot_title, fontsize = 20)
	plt.xlabel('Time', fontsize = 20)
	plt.ylabel('Number of groups, '+r'$x$', fontsize = 20)	
	plt.xticks(fontsize = 15)	
	plt.yticks(fontsize = 15)
		
	return 0

def FinalPie(FinalState):
	plt.pie(FinalState, colors = ['#e41a1c', '#4daf4a', '#377eb8'], startangle = 90, 
							   wedgeprops=dict(width=0.5, alpha = 0.8, edgecolor = 'w', lw = 3), 
							   counterclock = False)

	return 0



Flag_SaveFigs = 1


FullLCset = LCS.PartList(19)
FocalLC = 6
LC = FullLCset[FocalLC]
MaxSize = LCS.MinMaxSizes(LC)[1]


Times2Record = np.arange(0, 6, 0.1)
b_values = [3.0, 2.0, 1.0]

InitState = np.zeros(MaxSize)
InitState[0] = 1.0


""" panel one - no competition """

b = LCS.BirthInit('direct', MaxSize, b_values)
d = LCS.DeathInit('', MaxSize, [])
K = 0*LCS.InteractionInit('none', MaxSize, [])
A = LCS.ProjectionMatrix(LC, b, d, MaxSize)

PopRecord = LCS.LifeCyclesCompetition([A], K, [InitState], Times2Record)[0]
FinalState = PopRecord[-1,:]

PlotDemDyn(Times2Record, PopRecord, plot_title = 'No competition, '+r'$K_{ij} = 0$')
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/NoComp_dyn_v2.png', dpi=300, bbox_inches='tight')
plt.show()

FinalPie(FinalState)
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/NoComp_pie.png', dpi=300, bbox_inches='tight')
plt.show()






""" panel two - neutral competition kernel """

b = LCS.BirthInit('direct', MaxSize, b_values)
d = LCS.DeathInit('', MaxSize, [])
K = 0.1*LCS.InteractionInit('none', MaxSize, [])
A = LCS.ProjectionMatrix(LC, b, d, MaxSize)

PopRecord = LCS.LifeCyclesCompetition([A], K, [InitState], Times2Record)[0]
FinalState = PopRecord[-1,:]


PlotDemDyn(Times2Record, PopRecord, plot_title = 'Constant kernel, '+r'$K_{ij} = 0.1$')
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/NeutralComp_dyn_v2.png', dpi=300, bbox_inches='tight')
plt.show()
FinalPie(FinalState)
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/NeutralComp_pie.png', dpi=300, bbox_inches='tight')
plt.show()




""" panel three - killer competition kernel """

b = LCS.BirthInit('direct', MaxSize, b_values)
d = LCS.DeathInit('', MaxSize, [])
K = LCS.InteractionInit('kernel_killer', MaxSize, [0.5, 0.1, 0.1])
A = LCS.ProjectionMatrix(LC, b, d, MaxSize)

PopRecord = LCS.LifeCyclesCompetition([A], K, [InitState], Times2Record)[0]
FinalState = PopRecord[-1,:]


PlotDemDyn(Times2Record, PopRecord, ylimits = [-.5, 10.5], plot_title = 'Killer kernel, '+r'$K_{ij} = k_j$')
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/KillerComp_dyn_v2.png', dpi=300, bbox_inches='tight')
plt.show()

FinalPie(FinalState)
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/KillerComp_pie.png', dpi=300, bbox_inches='tight')
plt.show()






""" panel four - arbitrary (victim) competition kernel """

b = LCS.BirthInit('direct', MaxSize, b_values)
d = LCS.DeathInit('', MaxSize, [])
K = LCS.InteractionInit('kernel_victim', MaxSize, [1.0, 0.0, 0.0])
A = LCS.ProjectionMatrix(LC, b, d, MaxSize)

PopRecord = LCS.LifeCyclesCompetition([A], K, [InitState], Times2Record)[0]
FinalState = PopRecord[-1,:]


PlotDemDyn(Times2Record, PopRecord, plot_title = 'Arbitrary kernel')
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/VictimComp_dyn_v2.png', dpi=300, bbox_inches='tight')
plt.show()

FinalPie(FinalState)
if Flag_SaveFigs == 1:
	plt.savefig('../SourceFigs/VictimComp_pie.png', dpi=300, bbox_inches='tight')
plt.show()



#
#
#
#LC_Binary_1 = 0
#LC_Binary_2 = 1
#LC_Mult = 2
#
##LC_Binary_1 = 1
##LC_Binary_2 = 3
##LC_Mult = 5
#
##JobID = int(sys.argv[1]) 
#JobID = 0
#
##FileName = 'Results/Anomalous_r1_LC_1_3_5_file_'+str(JobID)+'.txt'
#
#
#LCset = [FullLCset[LC_Binary_1], FullLCset[LC_Binary_2], FullLCset[LC_Mult]]
#MaxSize = LCS.MinMaxSizes(LCset[-1])[1]
#
#""" initialize by the dominance triplet """
#stop = 0
#b_r = LCS.BirthInit('bumped_rand_exp', MaxSize, [1, 0.1])
#d_r = LCS.DeathInit('', MaxSize, [])
#K_r = LCS.InteractionInit('bumped_rand_exp', MaxSize, [1, 0.1])
#
#
#
#X_1 = LCS.StationaryState(FullLCset[LC_Binary_1], b_r, d_r, K_r, MaxSize)
#X_2 = LCS.StationaryState(FullLCset[LC_Binary_2], b_r, d_r, K_r, MaxSize)
#X_M = LCS.StationaryState(FullLCset[LC_Mult], b_r, d_r, K_r, MaxSize)
#
#
#plt.plot(np.dot(K_r, X_1), '-or')
#plt.plot(np.dot(K_r, X_2), '-ob')
#plt.plot(np.dot(K_r, X_M), '-og')
#plt.show()


#""" LC_2 into LC_1 """
#X_B1 = LCS.StationaryState(FullLCset[LC_Binary_1], b_r, d_r, K_r, MaxSize)
#A_B2 = LCS.ProjectionMatrix(FullLCset[LC_Binary_2], b_r, d_r, MaxSize)
#InvasionRate_B2_in_B1 = LCS.InvasionRate(A_B2, K_r, X_B1)

#""" LC_1 into LC_2 """
#X_B2 = LCS.StationaryState(FullLCset[LC_Binary_2], b_r, d_r, K_r, MaxSize)
#A_B1 = LCS.ProjectionMatrix(FullLCset[LC_Binary_1], b_r, d_r, MaxSize)
#InvasionRate_B1_in_B2 = LCS.InvasionRate(A_B1, K_r, X_B2)

#X_Mult = LCS.StationaryState(FullLCset[LC_Mult], b_m, d_m, K_m, MaxSize)


#
#if (InvasionRate_B1_in_B2>0)and(InvasionRate_B2_in_B1<0):
#	stop = 1
#
#
#X_Mult = LCS.StationaryState(FullLCset[LC_Mult], b_r, d_r, K_r, MaxSize)
#Anomality = LCS.InvasionRate(A_B2, K_r, X_Mult)
#print("Anomality = ", Anomality)
#
#
#""" serach for an anomalous triplet """
#stop = 0
#while stop == 0:
#	b_m = np.power(10.0, np.log10(b_r) + np.random.normal(scale = 0.5, size = MaxSize))
#	d_m = deepcopy(d_r)
#	K_m = np.power(10.0, np.log10(K_r) + np.random.normal(scale = 0.5, size = (MaxSize, MaxSize) ))
#	
#	""" test LC_2 into LC_1 """
#	X_B1 = LCS.StationaryState(FullLCset[LC_Binary_1], b_m, d_m, K_m, MaxSize)
#	A_B2 = LCS.ProjectionMatrix(FullLCset[LC_Binary_2], b_m, d_m, MaxSize)
#	InvasionRate_B2_in_B1 = LCS.InvasionRate(A_B2, K_m, X_B1)
#	
#	""" test LC_1 into LC_2 """
#	X_B2 = LCS.StationaryState(FullLCset[LC_Binary_2], b_m, d_m, K_m, MaxSize)
#	A_B1 = LCS.ProjectionMatrix(FullLCset[LC_Binary_1], b_m, d_m, MaxSize)
#	InvasionRate_B1_in_B2 = LCS.InvasionRate(A_B1, K_m, X_B2)
#	
#	""" test LC_2 into LC_Mult """
#	X_Mult = LCS.StationaryState(FullLCset[LC_Mult], b_m, d_m, K_m, MaxSize)
#	InvasionRate_B2_in_Mult = LCS.InvasionRate(A_B2, K_m, X_Mult)
#	
#	if (InvasionRate_B1_in_B2>0)and(InvasionRate_B2_in_B1<0)and(InvasionRate_B2_in_Mult > Anomality):
#		b_r = deepcopy(b_m)
#		K_r = deepcopy(K_m)
#		Anomality = InvasionRate_B2_in_Mult
#		print("Anomality = ", Anomality)
#	
##	if Anomality > 0:
##		stop = 1
#	
#	
#
#file = open(FileName, 'w')
#
#headerline = ''
#for i in np.arange(MaxSize):
#    headerline += 'b_'+str(i+1)+', '
#for i in np.arange(MaxSize):
#    headerline += 'd_'+str(i+1)+', '
#for i in np.arange(MaxSize):
#    for j in np.arange(MaxSize):
#        headerline += 'K_'+str(i+1)+'_'+str(j+1)+', '
#headerline += 'dummy'+'\n'
#
#file.write(headerline)
#
#line2write = ''
#for i in np.arange(MaxSize):
#	line2write += str(b_r[i])+', '
#for i in np.arange(MaxSize):
#	line2write += str(d_r[i])+', '
#for i in np.arange(MaxSize):
#	for j in np.arange(MaxSize):
#		line2write += str(K_r[i,j])+', '
#line2write += '-'+'\n'
#file.write(line2write)
#
#file.close()







