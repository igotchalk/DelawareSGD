import sgems
import random
import time


date = time.strftime("%Y%m%d")
date = date[2:] # format 170101 for Jan 01, 2017


#Parameters to set
useHD = 0 #  constrained to hard data? Y/N = 1/0
pSand = [0.5, 0.67, 0.72, 0.80]  #  P(sand) marginal distributions 
pClay = [0.5, 0.33, 0.28, 0.20]
skipVals = []  #Values to skip in simulation (e.g. 0.5)
reals = 25
nodes = 30
TIs = ['TI_170112_50_stack__0__0__0__0__0__0__0__0__real0','TI_170112_67_stack__0__0__0__0__0__0__0__0__real0','TI_170112_72_stack1__0__0__0__0__0__0__0__0__0__real0','TI_170112_80_stack__0__0__0__0__real0'] #must be in line with probSand
softProb1 = ['pSand_170113_50_tr1','pSand_170113_67_tr1','pSand_170113_72_tr1','pSand_170113_80_tr1']  #must be in line with probSand
softProb0 = ['complement(pSand_170113_50_tr1)','complement(pSand_170113_67_tr1)','complement(pSand_170113_72_tr1)','complement(pSand_170113_80_tr1)'] #must be in line with probSand

# assign HD or keep empty
if useHD==1:
	HD = '<Hard_Data  grid="Wellb" region="" property="fb - categorical"  />'
else:
	HD = '<Hard_Data  grid="" region="" property=""  />'



# Main loop: run snesim
for i in range(len(pSand)):
	if pSand[i] in skipVals:
		print('loop number' + str(i))
	else:
		print(str(i))
		seed =  random.randrange(1000, 14071789)
		sgems.execute('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     ' + HD + '     <use_pre_simulated_gridded_data  value="0"  />     <Use_ProbField  value="1"  />     <ProbField_properties count="2"   value="' + softProb0[i] + ';' + softProb1[i] + '"  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="0" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <GridSelector_Sim value="ARR_Grid2" region=""  />     <Property_Name_Sim  value="snesim_170113_loop' + str(int(pSand[i]*100)) + '" />     <Nb_Realizations  value="' + str(reals) + '" />     <Seed  value="' + str(seed) + '" />     <PropertySelector_Training  grid="TIgrid_170112" region="" property="' + TIs[i] + '"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="' + str(pClay[i]) + ' ' + str(pSand[i]) +'" />     <Max_Cond  value="' + str(nodes) + '" />     <Search_Ellipsoid  value="80 80 20  0 0 0" />  </parameters>')
		sgems.execute('SaveGeostatGrid ARR_Grid2::E:\SGEMS\WellProject_' + date + '_postSnesim' + str(int(pSand[i]*100)) + '.sgems::sgems')	

