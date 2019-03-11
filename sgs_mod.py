from __future__ import division
import sys
import pandas as pd
import numpy as np
import os
from subprocess import call
import datetime
from pathlib import Path

def grid_string(size, cells):
    cell_size = [x/y for x, y in zip(size, cells)]
    sgems_string = 'myGrid::'+str(cells[0])+'::'+str(cells[1])+'::'+str(cells[2])+'::'+str(cell_size[0])+'::'+str(cell_size[1])+'::'+str(cell_size[2])+'::0::0::0'
    return sgems_string

def make_target_file(target):
    filepath = os.getcwd() + '\\target.txt'
    filepath = filepath.replace('\\', '/')
    s = np.random.normal(target[0], target[1], 10000)
    np.savetxt('target.txt', s,  fmt='%.4f')
    return filepath

def snesim_grid(name, workdir,grid_size, grid_cells,
                TIfile=None,TIname=None,search_ellipse=(80,80,10,0,0,0),
                marg=0.5,seed=100,nodes=10,nreals=1,output=False,
                rmfiles=False,rotind=None,constrain=0):
    
    if rotind is None:
        if np.any(np.diff(search_ellipse[:3]) > 0):
            print('need to rotate grid for sgems')
            #Need flip grid around for sgems
            rotind = rot_mod2sgems(search_ellipse[:3])
        else:
            print('no rotation needed for sgems')
            rotind = np.arange(len(grid_size))
            
    if np.any(np.diff(search_ellipse[:3]) > 0):
        print('Please enter search ellipse as (max>=med>=min,az,dip,rake)... \n '
              '...Reordering to run. Please fix for next iteration')
        search_ellipse = np.r_[np.flip(np.sort(search_ellipse[:3])),
                                               np.flip(np.sort(search_ellipse[3:]))]
    grid_size = np.r_[grid_size][rotind]
    grid_cells_sgems = np.r_[grid_cells][rotind]

    fname = name
    fname_commands = workdir.joinpath(fname + '_commands.txt').as_posix()
    marg = round(marg,2)
    comp_marg = round(1-marg,2) #complement of marg
    search_str = ' '.join(str(x) for x in search_ellipse)

    #Make a grid
    sgems_grid = grid_string(grid_size, grid_cells_sgems)
    f = open(fname_commands, 'w')
    f.write('NewCartesianGrid ' + sgems_grid + '\n')
    
    #Load TI file
    if TIfile is None:
        f.write('LoadObjectFromFile  E:/SGEMS/TIs_try2::All \n')
        TIgrid = 'TIgrid_170112'
        TIname = 'TI_170112_50_stack__0__0__0__0__0__0__0__0__real0'
    else:
        f.write('LoadObjectFromFile ' + TIfile.as_posix() +'::All \n')
    
    #Write snesim XML
    f.write('RunGeostatAlgorithm  snesim_std::/GeostatParamUtils/XML::<parameters>  <algorithm name="snesim_std" />     <GridSelector_Sim value="myGrid" region=""  />     <Property_Name_Sim  value="'+ fname + '"/>     <Nb_Realizations  value="'+str(nreals)+'" />     <Seed  value="' + str(seed) + '" />    <PropertySelector_Training  grid="'+ TIgrid + '" region="" property="' + TIname + '"  />     <Nb_Facies  value="2" />     <Marginal_Cdf  value="' + str(comp_marg) + ' ' + str(marg) + '" />     <Max_Cond  value="' + str(nodes) + '" />     <Search_Ellipsoid  value="'+ search_str +'" />    <Hard_Data  grid="" region="" property=""  />     <use_pre_simulated_gridded_data  value="0"  />     <Use_ProbField  value="0"  />     <ProbField_properties count="0"   value=""  />     <TauModelObject  value="1 1" />     <use_vertical_proportion  value="0"  />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="'+str(constrain)+'" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="3" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />   </parameters>'+'\n')
    
    #Save grid
    expname = workdir.joinpath(fname+'.csv').as_posix()
    f.write('SaveGeostatGrid  myGrid::'+expname+'::csv::0::'+fname+'__real0'+'\n')
    f.close()
    
    #Call SGeMs
    call(['sgems-x64', fname_commands], shell=False)
    
    if output:
        #outmat = pd.read_csv(expname, header = 0)#, names=[fname])
        outmat = read_sgems_grid(expname,grid_cells_sgems,categorical=True,plotyn=False,rotind=rotind)
    else:
        outmat=None
    if rmfiles:
        os.remove(expname)
        os.remove(fname_commands)
        expname=None
    return Path(expname),outmat,rotind,grid_cells_sgems



def read_sgems_grid(filepath,grid_cells,grid_cells_sgems,categorical=True,plotyn=True,rotind=None):
    #takes an sgems output grid and changes it to a more plotable format
    #each column is read as a model and changed to the dimensions supplied by grid_cells
    
    def coding(col, codeDict):
        colCoded = pd.Series(col, copy=True)
        for key, value in codeDict.items():
            colCoded.replace(key, value, inplace=True)
        return colCoded
    
    df = pd.read_csv(filepath)
    if categorical:
        dtype=np.int
        for k in df.keys():
            df[k] = coding(df[k], {'category_0':0,'category_1':1})
    else:
        dtype=np.float
    if rotind is None:
        rotind = np.arange(len(grid_cells_sgems))
        
    grid_cells_exp = np.r_[grid_cells][np.argsort(grid_cells_sgems)]
    
    nmodels = len(df.columns)
    outdim = tuple(np.r_[nmodels,grid_cells])
    outmat = np.zeros(outdim,dtype=dtype)
    for i,k in enumerate(df.keys()):
        #outmat[i] = df[k].as_matrix().reshape(grid_cells2)
        rotmat =  df[k].as_matrix().reshape(grid_cells_exp)
        outmat[i] = np.rot90(rotmat,k=1,axes=(1,2))
        '''
    if plotyn:
        import matplotlib.pyplot as plt
        if len(grid_cells2)==3:
            plt.imshow(outmat[0,:,0,:])
            plt.xlabel('col number')
            plt.ylabel('lay number')
        elif len(grid_cells2)==2:
            plt.imshow(outmat[0])
        plt.title('first column of sgems csv')
        plt.show()
        '''
    return outmat

def rot_mod2sgems(inval):
    if len(np.r_[inval].shape) > 1:
        #take (nlay,nrow,ncol) format to (max,med,min), which sgems requires
        rotind = np.flip(np.argsort(inval.shape)) #sorted max to min
        sgems_grid = np.transpose(inval,axes=rotind)
        return sgems_grid,rotind
    else: #assume a list, tuple, or array of the dimensions
        rotind = np.flip(np.argsort(inval))
        return rotind

def rot_sgems2mod(inval,rotind):
    newind = np.argsort(rotind)
    if len(np.r_[inval].shape) > 1:
        return np.transpose(inval,newind)
    else:
        outval = np.r_[inval]
        return outval[newind]



def main():
    fname = 'sgs_grid'
    grid_size = [100, 100]
    grid_cells = [100, 100]
    #parameters = [(search ellipsoid), (variogram model)]
    parameters = [(12, 10, 10, 10, 0, 0, 0), (1, 'Gaussian', 5, 1, 1, 90, 0, 0), (100)]
    perm = sgs_grid('temp', grid_size, grid_cells, parameters, target=False) #[15, 5])
    print(perm)
    sys.exit(1)


if __name__=='__sgs_grid__':
    sgs_grid(*args)

if __name__=='__main__':
    main()
