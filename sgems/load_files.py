# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 14:18:04 2017

@author: ianpg
"""


import sgems
import random
import time


loadData = 1 #Should SGEMS load .sgems files into your workspace? 1/0

#Load data for SGEMS
if loadData==1:
      sgems.execute('LoadObjectFromFile  E:/SGEMS/well_data.sgems::All')
      sgems.execute('LoadObjectFromFile  E:/SGEMS/TIs_try2::All')
      sgems.execute('LoadObjectFromFile  E:/SGEMS/softProbJef_ARRgrid.sgems::All')

# Pre-loop declaration
print('Running loop_sisim. Pre-loop declaration...')

#Parameters
date = time.strftime("%Y%m%d")
date = date[2:] # format 170101 for Jan 01, 2017

useHD = 1 #assign well hard data?
skipglobal = 0 #skip the global iteration number listed
pSand = [.55, .65, .75]  #global %sand bin centers (each bin is 0.1 wide)
pA = [0.03,0.64,0.33] #probability of each global sand %, calculated from 100 snesim realizations and Bayes Rule
Ntot = 100  #total realizations
reals = [round(Ntot*pA[i],0) for i in range(len(pA))] #number of reals for each global sand %


simcount = 0 #If you want to start at realization j, set simcount to j-1. Does not reset after each global loop
SD = ['softProbJefA_C','softProbJefA']
TIs = ['TI_170112_50_stack__0__0__0__0__0__0__0__0__real0','TI_170112_67_stack__0__0__0__0__0__0__0__0__real0','TI_170112_80_stack__0__0__0__0__real0'] #must be in line with probSand
nodes = 10
#/Parameters


# assign HD or keep empty
#if useHD==1:
#      HD = '<Hard_Data_Grid value="Wellb" region=""  />     <Hard_Data_Property  value="fb - categorical"  />     <Assign_Hard_Data  value="1"  />'
#else:
#      HD = '<Hard_Data_Grid value="" region=""  />     <Hard_Data_Property  value=""  />     <Assign_Hard_Data  value="1"  />'


print('using HD? = ' + str(useHD))
# Main loop
print('Params:\nTotal realizations '+str(Ntot)+'\nRealizations per group '+str(reals)+'\nGlobal Percent in each group '+str(pSand)+'\nNumber of nodes '+str(nodes)+'\nTraining Images '+str(TIs)+'\nOutput file prefix E:\SGEMS\snesimPosteriorJ__')
print('Begin main loop...')