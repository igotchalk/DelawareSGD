# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 14:47:26 2017

@author: ianpg
"""

#Parameters
date = time.strftime("%Y%m%d")
date = date[2:] # format 170101 for Jan 01, 2017

useHD = 1 #assign well hard data?
skipglobal = 99 #skip the global iteration number listed
pSand = [.55, .65, .75]  #global %sand bin centers (each bin is 0.1 wide)
pA = [0.03,0.64,0.33] #probability of each global sand %, calculated from 100 snesim realizations and Bayes Rule
Ntot = 1000  #total realizations
reals = [round(Ntot*pA[i],0) for i in range(len(pA))] #number of reals for each global sand %


simcount = 0 #If you want to start at realization j, set simcount to j-1. Does not reset after each global loop
SD = ['softProbJefA_C','softProbJefA']
TIs = ['TI_170112_50_stack__0__0__0__0__0__0__0__0__real0','TI_170112_67_stack__0__0__0__0__0__0__0__0__real0','TI_170112_80_stack__0__0__0__0__real0'] #must be in line with probSand
nodes = 60
#/Parameters

