# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Tue Aug 18

@author: Hyejin Park
"""


import numpy as np
import sys
def WriteFile(file):
	jobname = ("LC_%d" %ID)
	fp = open(file, "w")
	fp.write("#!/bin/sh\n")
	fp.write("#Set your minimum acceptable walltime, format: day-hours:minutes:seconds\n")
	fp.write("#SBATCH --time=3-00:00:00\n")
	fp.write("#Set name of job shown in squeue\n")
	fp.write("#Request CPU resources\n")
	fp.write("#SBATCH --ntasks=1\n")
	fp.write("#SBATCH --ntasks-per-node=1\n")
	fp.write("#SBATCH --cpus-per-task=1\n")
	fp.write("#SBATCH -o ./outputs/output.%s\n" %jobname)
	fp.write("#SBATCH -e ./errors/error.%s\n" %jobname)
	fp.write("python Seven_dynamics.py %d\n" %ID)
	fp.close()


#
#""" 1. allocate min max and desired number of samples in each dimension + arrays"""
#Nsample = 40
#minrho = 0.01
#maxrho = 100
#minf = 1./10.
#maxf = 10.
#Nen = 100
#idx = np.zeros(Nsample+1)
#lnrho = np.zeros(Nsample+1)
#lnf = np.zeros(Nsample+1)
#rho = np.zeros(Nsample+1)
#f = np.zeros(Nsample+1)
#for i in range(Nsample+1):
#	idx[i] = i
#
#
#""" 2. calculate equal distance in ln(f)-ln(rho) space """
#dlnrho = ( np.log10(maxrho) - np.log10(minrho) )/Nsample
#dlnf = ( np.log10(maxf) - np.log10(minf) )/Nsample
#lnrho = np.log10(minrho) + dlnrho*idx
#lnf = np.log10(minf) + dlnf*idx
#
#
#""" 3. transform lnrho-lnf to rho-f to t1-t2 """
#rho = np.exp(lnrho/np.log10(np.e))
#f = np.exp(lnf/np.log10(np.e))
#coor= []
#for i in range(len(rho)):
#	for j in range(len(f)):
#		t1 = rho[i]/(1+f[j])
#		t2 = rho[i]*f[j]/(1+f[j])
#		coor.append([t1, t2, lnrho[i], lnf[j]])

""" 4. write simulation codes """
ofp = open("run.sh", "w")
ofp.write("module load python/3.7.4 \n")
#coorid = open("id_coordi_2.txt", "w")
for ID in range(100):
	file = ( "script/EELC_seven_dyn_%d.sh" % ID)
#	coorid.write("%d\t%g\t%g\t%g\t%g\n" % (i, coor[i][0], coor[i][1], coor[i][2], coor[i][3]) )
	ofp.write("sbatch %s\n" % file )
	WriteFile(file)
ofp.close()
