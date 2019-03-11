# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 19:10:11 2017

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
import timeit

loadData = 0 #Should SGEMS load .sgems files into your workspace? 1/0

#Load data for SGEMS
if loadData==1:
      sgems.execute('LoadObjectFromFile  E:/SGEMS/well_data.sgems::All')
      sgems.execute('LoadObjectFromFile  E:/SGEMS/pSand::All')

# Pre-loop declaration
print('Running loop_sisim. Pre-loop declaration...')

date = time.strftime("%Y%m%d")
date = date[2:] # format 170101 for Jan 01, 2017


export_name = 'sisim_kmeans_Thomas' 
skipglobal = [99] #skip the global iteration number(s) listed
useHD = 1 #assign well hard data?
HD_gridname = 'sisim_forSP_kmeansT'
pSand = [.55, .65, .75]  #global %sand bin centers (each bin is 0.1 wide)
pA = [0.03,0.64,0.33] #probability of each global sand %, calculated from 100 snesim realizations and Bayes Rule
Ntot = 100  #total realizations: MUST match number of simualated HD pointsets
reals = [round(Ntot*pA[i],0) for i in range(len(pA))] #number of reals for each global sand %
simcount = 0

# assign HD or keep empty
#if useHD==1:
#      HD = '<Hard_Data_Grid value="Wellb" region=""  />     <Hard_Data_Property  value="fb - categorical"  />     <Assign_Hard_Data  value="1"  />'
#else:
#      HD = '<Hard_Data_Grid value="" region=""  />     <Hard_Data_Property  value=""  />     <Assign_Hard_Data  value="1"  />'


print('using HD? = ' + str(useHD))
# Main loop
print('Begin main loop...')
for i in range(len(pSand)):
      if i in skipglobal:
          print('skipping iteration ' + str(i) + '...')
      else:
          sgems.execute('LoadObjectFromFile  E:/SGEMS/large_grid.sgems::All')
          print('gl0bal loop number ' + str(i) + ': pct sand bin centered at: ' + str(pSand[i]))
          #internal loop; new marginal probability for each realization
          for j in range(int(reals[i])):
                start_time = timeit.default_timer()
                simcount = simcount + 1
                print('simulation number ' + str(simcount) + ' out of ' + str(Ntot))
                seed =  random.randrange(1000, 14071789)  #seed generator so that reals aren't identical
                marg = round(float(random.randrange(int(pSand[i]*100-5),int(pSand[i]*100+4)))/100,2) #draw a random value from the histogram bin
                comp_marg = round(1-marg,2) #complement of marg
                HD = '<Hard_Data_Grid value="' + HD_gridname + '" region=""  />     <Hard_Data_Property  value="simulation' + str(simcount) + '"  />'
                print('gl0bal percent sand: ' + str(marg) + '   ' + str(comp_marg))
                sgems.execute('RunGeostatAlgorithm  sisim::/GeostatParamUtils/XML::<parameters>  <algorithm name="sisim" />     <Variogram_Median_Ik  nugget="0.15" structures_count="2"  >    <structure_1  contribution="0.74"  type="Gaussian"   >      <ranges max="54"  medium="24"  min="4"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>    <structure_2  contribution="0.11"  type="Exponential"   >      <ranges max="69"  medium="20"  min="4"   />      <angles x="0"  y="0"  z="0"   />    </structure_2>  </Variogram_Median_Ik>    <Variogram_Full_Ik  nugget="0" structures_count="1"  >    <structure_1  contribution="0"  type="Spherical"   >      <ranges max="0"  medium="0"  min="0"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram_Full_Ik>    <Variogram_Full_Ik_2  nugget="0" structures_count="1"  >    <structure_1  contribution="0"  type="Spherical"   >      <ranges max="0"  medium="0"  min="0"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram_Full_Ik_2>    ' + HD + '     <coded_grid value="" region=""  />     <coded_props count="0"   value=""  />     <Max_Conditioning_Data  value="20" />     <Search_Ellipsoid  value="80 80 10  0 0 0" />    <AdvancedSearch  use_advanced_search="0"></AdvancedSearch>    <Grid_Name value="ARR_largeGrid" region=""  />     <Property_Name  value="' + export_name + '_' + date + '_' + str(int(pSand[i]*100)) + '_' + str(simcount) + '" />     <Nb_Realizations  value="1" />     <Seed  value="' + str(seed) + '" />     <Categorical_Variable_Flag  value="1"  />     <Nb_Indicators  value="2" />     <Marginal_Probabilities  value="' + str(comp_marg) + ' ' + str(marg) +'" />     <lowerTailCdf  function ="Power"  extreme ="0"  omega ="3" />    <upperTailCdf  function ="Power"  extreme ="0"  omega ="0.333" />    <Median_Ik_Flag  value="1"  />     <Full_Ik_Flag  value="0"  />   </parameters>   ')
                elapsed = timeit.default_timer() - start_time
                print('     ...took ' + str(round(elapsed)) + ' seconds')
          #clean up: save grid, delete it, and reload the original version (inefficient, I know...)
          sgems.execute('SaveGeostatGrid ARR_largeGrid::E:\SGEMS\export_' + export_name + '_' + date + '_' + str(int(pSand[i]*100)) + '.sgems::sgems')	
          sgems.execute('DeleteObjects  ARR_largeGrid')

