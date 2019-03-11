# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:53:45 2017
Runs sisim on a new HD grid with each iteration


@author: ianpg
"""

"""
Executes sisim algorithm on SGEMS in a loop, varying grid parameters.
For each iteration, the grid is saved as an .sgems file, and the grid is then deleted to save space.


Output:.sgems files in E:\SGEMS\ 
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
skipglobal = 99 #skip the global iteration number listed
pSand = [.55, .65, .75]  #global %sand bin centers (each bin is 0.1 wide)
pA = [0.03,0.64,0.33] #probability of each global sand %, calculated from 100 snesim realizations and Bayes Rule
Ntot = 50  #total realizations
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
for i in range(len(pSand)):
      if i==skipglobal:
          print('skipping iteration ' + str(i) + '...')
      else:
          print('gl0bal loop number ' + str(i) + ': pct sand bin centered at: ' + str(pSand[i]))
          #internal loop; new marginal probability for each realization
          for j in range(int(reals[i])):
                simcount = simcount + 1
                print('simulation number ' + str(simcount) + ' out of ' + str(Ntot))
                seed =  random.randrange(1000, 14071789)  #seed generator so that reals aren't identical
                marg = round(float(random.randrange(int(pSand[i]*100-5),int(pSand[i]*100+4)))/100,2) #draw a random value from the histogram bin
                comp_marg = round(1-marg,2) #complement of marg
                print('gl0bal percent sand: ' + str(marg) + '   ' + str(comp_marg))
                #normal simulation (not working) --> sgems.execute('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="0" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <Hard_Data  grid="Wellb" region="" property="fb - categorical"  />     <use_pre_simulated_gridded_data  value="0"  />     <Use_ProbField  value="1"  />     <ProbField_properties count="2"   value="' + SD[0] + ';' + SD[1] + '"  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />     <GridSelector_Sim value="ARR_Grid2" region=""  />     <Property_Name_Sim  value="snesim_' + date + '_posteriorJ_' + str(int(pSand[i]*100)) + '_' + str(simcount) + '" />     <Nb_Realizations  value="1" />     <Seed  value="' + str(seed) + '" />     <PropertySelector_Training  grid="TIgrid_170112" region="" property="' + TIs[i] + '"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="' + str(comp_marg) + ' ' + str(marg) + '" />     <Max_Cond  value="' + str(nodes) + '" />     <Search_Ellipsoid  value="80 80 10  0 0 0" />  </parameters>   ')
                print('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="0" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <Hard_Data  grid="Wellb" region="" property="fb - categorical"  />     <use_pre_simulated_gridded_data  value="0"  />     <Use_ProbField  value="1"  />     <ProbField_properties count="2"   value="' + SD[0] + ';' + SD[1] + '"  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />     <GridSelector_Sim value="ARR_Grid2" region=""  />     <Property_Name_Sim  value="snesim_' + date + '_posteriorJ_' + str(int(pSand[i]*100)) + '_' + str(simcount) + '" />     <Nb_Realizations  value="1" />     <Seed  value="' + str(seed) + '" />     <PropertySelector_Training  grid="TIgrid_170112" region="" property="' + TIs[i] + '"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="' + str(comp_marg) + ' ' + str(marg) + '" />     <Max_Cond  value="' + str(nodes) + '" />     <Search_Ellipsoid  value="80 80 10  0 0 0" />  </parameters>   ')
                #static simulation --> sgems.execute('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="0" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <GridSelector_Sim value="ARR_Grid2" region=""  />     <Property_Name_Sim  value="test" />     <Nb_Realizations  value="1" />     <Seed  value="211175" />     <PropertySelector_Training  grid="TIgrid_170112" region="" property="TI_170112_50_stack__0__0__0__0__0__0__real0"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="0.6 0.4" />     <Max_Cond  value="10" />     <Search_Ellipsoid  value="80 80 10  0 0 0" />    <Hard_Data  grid="Wellb" region="" property="fb - categorical"  />     <use_pre_simulated_gridded_data  value="0"  />     <Use_ProbField  value="1"  />     <ProbField_properties count="2"   value="softProbJefA_C;softProbJefA"  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />   </parameters>   ')
                sgems.execute('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="0" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <Hard_Data  grid="Wellb" region="" property="fb - categorical"  />     <use_pre_simulated_gridded_data  value="0"  />     <Use_ProbField  value="1"  />     <ProbField_properties count="2"   value="softProbJefA_C;softProbJefA"  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />     <GridSelector_Sim value="ARR_Grid2" region=""  />     <Property_Name_Sim  value="snesim_170212_posteriorJ_75_68" />     <Nb_Realizations  value="1" />     <Seed  value="1445868" />     <PropertySelector_Training  grid="TIgrid_170112" region="" property="TI_170112_80_stack__0__0__0__0__real0"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="0.26 0.74" />     <Max_Cond  value="10" />     <Search_Ellipsoid  value="80 80 10  0 0 0" />  </parameters>   ')
          #clean up: save grid, delete it, and reload the original version (inefficient, I know...)
          sgems.execute('SaveGeostatGrid ARR_Grid2::E:\SGEMS\snesimPosteriorJ__' + date + '_' + str(int(pSand[i]*100)) + '.sgems::sgems')	
          sgems.execute('DeleteObjects  ARR_Grid2')
          sgems.execute('LoadObjectFromFile  E:/SGEMS/softProbJef_ARRgrid.sgems::All')


          
          
          
          
          
          
          
          
          
          
          
          