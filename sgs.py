from __future__ import division
import sys
import pandas as pd
import numpy as np
import os
from subprocess import call
import distutils.spawn

def grid_string(size, cells):
    cell_size = [x/y for x, y in zip(size, cells)]
    sgems_string = 'myGrid::'+str(cells[0])+'::'+str(cells[1])+'::'+str(1)+'::'+str(cell_size[0])+'::'+str(cell_size[1])+'::1.0::0::0::0'
    return sgems_string

def make_target_file(target):
    filepath = os.getcwd() + '\\target.txt'
    filepath = filepath.replace('\\', '/')
    s = np.random.normal(target[0], target[1], 10000)
    np.savetxt('target.txt', s,  fmt='%.4f')
    return filepath

def sgs_grid(name, grid_size, grid_cells, params, target=None):
    fname = name + '_sgems'
    fname_commands = fname + '_commands.txt'
    search_ellipse = params[0]
    variogram_model = params[1]
    seed = params[2]
    if seed == 'random':
        seed = int(np.random.uniform(1, 9999999))
    search = ' '.join(str(x) for x in search_ellipse[1:])
    output = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
### No translation of origin
    sgems_grid = grid_string(grid_size, grid_cells)
    f = open(fname_commands, 'w')
    f.write('NewCartesianGrid ' + sgems_grid + '\n')
    if target:
        filepath = make_target_file(target)
        f.write('RunGeostatAlgorithm  sgsim::/GeostatParamUtils/XML::<parameters>  <algorithm name="sgsim" />     <Grid_Name value="myGrid" region=""  />     <Property_Name  value="' + fname + '" />     <Nb_Realizations  value="1" />     <Seed  value="'+str(seed)+'" />     <Kriging_Type  value="Simple Kriging (SK)"  />     <Trend  value="0 0 0 0 0 0 0 0 0 " />    <Local_Mean_Property  value=""  />     <Assign_Hard_Data  value="0"  />     <Hard_Data  grid="" region="" property=""  />     <Max_Conditioning_Data  value="'+str(search_ellipse[0])+'" />     <Search_Ellipsoid  value="'+search+'" />    <AdvancedSearch  use_advanced_search="0"></AdvancedSearch>    <Use_Target_Histogram  value="1"  />     <nonParamCdf  ref_on_file ="1"  ref_on_grid ="0"  break_ties ="1" filename ="'+filepath+'"   grid =""  property ="">  <LTI_type  function ="No extrapolation"  extreme ="0"  omega ="3" />  <UTI_type  function ="No extrapolation"  extreme ="0"  omega ="0.333" />  </nonParamCdf>    <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="' + str(variogram_model[0])+'"  type="'+str(variogram_model[1])+'"   >      <ranges max="'+str(variogram_model[2])+'"  medium="'+str(variogram_model[3])+'"  min="'+str(variogram_model[4])+'"   />      <angles x="'+str(variogram_model[5])+'"  y="'+str(variogram_model[6])+'"  z="'+str(variogram_model[7])+'"   />    </structure_1> </Variogram>  </parameters>'+'\n')
        #remove_target_file
    else:
        f.write('RunGeostatAlgorithm  sgsim::/GeostatParamUtils/XML::<parameters>  <algorithm name="sgsim" />     <Grid_Name value="myGrid" region=""  />     <Property_Name  value="' + fname + '" />     <Nb_Realizations  value="1" />     <Seed  value="'+str(seed)+'" />     <Kriging_Type  value="Simple Kriging (SK)"  />     <Trend  value="0 0 0 0 0 0 0 0 0 " />    <Local_Mean_Property  value=""  />     <Assign_Hard_Data  value="1"  />     <Hard_Data  grid="" region="" property=""  />     <Max_Conditioning_Data  value="'+str(search_ellipse[0])+'" />     <Search_Ellipsoid  value="'+search+'" />    <AdvancedSearch  use_advanced_search="0"></AdvancedSearch>    <Use_Target_Histogram  value="0"  />     <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="' + str(variogram_model[0])+'"  type="'+str(variogram_model[1])+'"   >      <ranges max="'+str(variogram_model[2])+'"  medium="'+str(variogram_model[3])+'"  min="'+str(variogram_model[4])+'"   />      <angles x="'+str(variogram_model[5])+'"  y="'+str(variogram_model[6])+'"  z="'+str(variogram_model[7])+'"   />    </structure_1> </Variogram>  </parameters>'+'\n')
    f.write('SaveGeostatGrid  myGrid::'+output+'\\'+fname+'.csv::csv::0::'+fname+'__real0'+'\n')
    f.close()
    call(['sgems', fname_commands], shell=True, env={"PATH":"C:\\SGeMS-x64-Beta\\", "GSTLAPPLIHOME":"C:\\SGeMS-x64-Beta\\"})
    sgs = pd.read_csv(fname + '.csv', header = 0, names=[name])
    os.remove(fname_commands)
    os.remove(fname + '.csv')
    return sgs



def main():
    fname = 'sgs_grid'
    grid_size = [15, 10]
    grid_cells = [100, 100]
    #parameters = [(search ellipsoid), (variogram model)]
    parameters = [(12, 10, 10, 10, 0, 0, 0), (1, 'Gaussian', 5, 1, 1, 90, 0, 0), (100)]

    perm = sgs_grid('temp', grid_size, grid_cells, parameters, [15, 5])
    print perm
    sys.exit(1)


if __name__=='__sgs_grid__':
    sgs_grid(*args)

if __name__=='__main__':
    main()
