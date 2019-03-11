# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 20:42:46 2017

@author: ianpg
"""
import sgems
import random
import time
import timeit


# Pre-loop declaration
print('Running loop_snesim_largegrid.py; Pre-loop initialization...')

#Parameters
date = time.strftime("%Y%m%d")
date = date[2:] # format 170101 for Jan 01, 2017

useHD = 1 #assign well hard data?
skipglobal = [99] #skip the global iteration number(s) listed
pSand = [.55, .65, .75]  #global %sand bin centers (each bin is 0.1 wide)
pA = [0.03,0.64,0.33] #probability of each global sand %, calculated from 100 snesim realizations and Bayes Rule
Ntot = 100  #total realizations
reals = [round(Ntot*pA[i],0) for i in range(len(pA))] #number of reals for each global sand %


simcount = 3 #If you want to start at realization j, set simcount to j-1. Does not reset after each global loop
TIs = ['TI_170112_50_stack__0__0__0__0__0__0__0__0__real0','TI_170112_67_stack__0__0__0__0__0__0__0__0__real0','TI_170112_80_stack__0__0__0__0__real0'] #must be in line with probSand
nodes = 10

#grid details
import_dir = 'E:/SGEMS/'
import_prefix = 'export_snesimPosteriorT_170227'
sim_prefix = 'snesimPosteriorT_170227'
suffixes = [55,65,75]
export_name = 'snesimPosteriorTLarge'
largegrid_name = 'ARR_largeGrid'
# /Parameters

i = 1
j = 0

#Load data
sgems.execute('LoadObjectFromFile  E:/SGEMS/well_data.sgems::All')
sgems.execute('LoadObjectFromFile  E:/SGEMS/TIs_try2::All')



#test loop
sgems.execute('LoadObjectFromFile  E:/SGEMS/large_grid.sgems::All') #load large grid
smallgrid_name = import_prefix + '_' + str(suffixes[i])
sgems.execute('LoadObjectFromFile  ' + import_dir + smallgrid_name + '.sgems::All') #load small grid simulations
print('gl0bal loop number ' + str(i) + ': pct sand bin centered at: ' + str(pSand[i]))
#internal loop; new marginal probability for each realization

#calculate snesim parameters 
start_time = timeit.default_timer()
simcount = simcount + 1
print('simulation number ' + str(simcount) + ' out of ' + str(Ntot))
seed =  random.randrange(1000, 14071789)  #seed generator so that reals aren't identical
marg = round(float(random.randrange(int(pSand[i]*100-5),int(pSand[i]*100+4)))/100,2) #draw a random value from the histogram bin
comp_marg = round(1-marg,2) #complement of marg
print('gl0bal percent sand: ' + str(marg) + '   clay: ' + str(comp_marg))

#Part 1: copy grid over
sim_name = sim_prefix + '_' + str(suffixes[i]) + '_' + str(simcount)
sgems.execute('CopyProperty  ARR_Grid2::' + sim_name + '__real0::' + largegrid_name + '::' + sim_name + '::0::0') #copy sim from small grid to large grid
#Part 2: run snesim on large grid  
sgems.execute('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="0" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <Hard_Data  grid="Wellb" region="" property="fb - categorical"  />     <use_pre_simulated_gridded_data  value="1"  />     <Pre_Simulated_Gridded_Data  value="' + sim_name + '"  />     <Use_ProbField  value="0"  />     <ProbField_properties count="0"   value=""  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />     <GridSelector_Sim value="' + largegrid_name + '" region=""  />     <Property_Name_Sim  value="' + sim_name + '" />     <Nb_Realizations  value="1" />     <Seed  value="' + str(seed) + '" />     <PropertySelector_Training  grid="TIgrid_170112" region="" property="' + TIs[i] + '"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="' + str(comp_marg) + ' ' + str(marg) + '" />     <Max_Cond  value="' + str(nodes) + '" />     <Search_Ellipsoid  value="80 80 10  0 0 0" />  </parameters>   ')
elapsed = timeit.default_timer() - start_time
print('     ...took ' + str(round(elapsed)) + ' seconds')
 
#clean up: save grid, delete it, and reload the original version (inefficient, I know...)
#sgems.execute('SaveGeostatGrid ' + largegrid_name + '::E:\SGEMS\export_' + export_name + '_' + date + '_' + str(int(pSand[i]*100)) + '.sgems::sgems')	
#sgems.execute('DeleteObjects  ' + largegrid_name + '')
#sgems.execute('DeleteObjects  ARR_Grid2')