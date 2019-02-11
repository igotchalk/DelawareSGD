
# coding: utf-8

import os
from pathlib import Path
import sys
import numpy as np
import numpy.matlib as matlib
import matplotlib.pyplot as plt
import matplotlib.colors
import warnings
import scipy.stats as sts

#%% Useful functions

def load_obj(dirname,name):
    import pickle
    with open(Path(dirname).joinpath(name + '.pkl').as_posix(), 'rb') as f:
        return pickle.load(f)

def save_obj(dirname,obj,name):
    import pickle
    with open(Path(dirname).joinpath(name + '.pkl').as_posix(), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

#Create new MC_file
def create_MC_file():
    import datetime
    ts = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M')
    MC_dir = Path(os.path.join(m.model_ws, 'MC_expt_' + ts))
    if not MC_dir.exists():
        MC_dir.mkdir()
    m.MC_file = MC_dir.joinpath('expt.txt')
    with m.MC_file.open('w') as wf:
        wf.close
    print(m.MC_file)
    return

#nearest value in array
def find_nearest(array,value):
    import numpy as np
    idx = (np.abs(array-value)).argmin()
    idx.astype('int')
    return array[idx]

#take distance in meters, return column in model
def loc_to_col(locs):
    cols = [int(find_nearest((np.arange(ncol)*delc),loc)) for loc in locs]
    return cols

#make a line across the grid
def get_line(start, end,allrows=1,nrow=None):
    """Bresenham's Line Algorithm
    Produces a list of tuples from start and end

    >>> points1 = get_line((0, 0), (3, 4))
    >>> points2 = get_line((3, 4), (0, 0))
    >>> assert(set(points1) == set(points2))
    >>> print points1
    [(0, 0), (1, 1), (1, 2), (2, 3), (3, 4)]
    >>> print points2
    [(3, 4), (2, 3), (1, 2), (1, 1), (0, 0)]
    """
    # Setup initial conditions
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1

    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1

    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1

    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        if allrows==1:
            if nrow is None:
                nrow = m.nrow
            for row in range(nrow):
                coord = (y, row, x) if is_steep else (x, row, y)
                points.append(coord)
        else:
            coord = (y, x) if is_steep else (x, y)
            points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()
    return points


#make all cells=0 above the line from get_line()
#Calculate freshwater head based on column of saltwater above each node (rho*g*z)
def shade_above(nlay,nrow,ncol,point_list,third_dim=1):
    import numpy as np
    grd = np.ones((nlay,nrow,ncol),dtype='int')
    ocean_hf = []
    if len(point_list)==0:
        return grd,ocean_hf
    for (lay,row,col) in point_list:
        grd[lay,:,col] = -1 #assign ocean ibound to -1
        grd[:lay,:,col] = 0 #assign cells above ocean to 0
        hf = densefresh/densesalt*ocean_elev - (densesalt - densefresh)/densefresh*(henry_botm[lay] +.5*delv)
        for irow in range(nrow):
            ocean_hf.append((int(lay),int(irow),int(col),hf))
    ocean_hf = tuple(np.array(ocean_hf).T)
    ocean_hf = (ocean_hf[0].astype('int'),
                ocean_hf[1].astype('int'),
                ocean_hf[2].astype('int'),
                ocean_hf[3])
    return grd,ocean_hf

def get_ocean_right_edge(m,ocean_line_tuple,startlay=None,col=None):
    import numpy as np
    point_list = []
    if col is None:
        col = m.ncol-1
    #If there is no vertical side boundary, return bottom-right corner node
    if len(ocean_line_tuple)==0:
        if startlay is None:
            startlay = 0
    elif max(ocean_line_tuple[0])==m.nlay:
        startlay = m.nlay
    elif max(ocean_line_tuple[0])<m.nlay:
        startlay = max(ocean_line_tuple[0])
    for lay in range(startlay,m.nlay):
        for row in range(m.nrow):
            point_list.append((lay,row,col))
    point_list = tuple(np.array(point_list).T)
    return point_list

def add_pumping_wells(wel_data,ssm_data,n_wells,flx,rowcol,kper):
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    new_weldata = wel_data
    new_ssmdata = ssm_data
    wel_cells = []
    for k in range(n_wells):
        row,col = rowcol[k]
        for i in range(nper):
            if i in kper:
                for j in range(nlay):
                    #WEL {stress_period: [lay,row,col,flux]}
                    new_weldata[i].append([j,row,col,-flx[k]*delv_weight[j]])
                    wel_cells.append((j,row,col))
                    #SSM: {stress_period: [lay,row,col,concentration,itype]}
                    new_ssmdata[i].append([j,row,col,Cfresh,itype['WEL']]) #since it's a sink, conc. doesn't matter
            else:
                for j in range(nlay):
                    #WEL {stress_period: [lay,row,col,flux]}
                    new_weldata[i].append([j,row,col,0])
                    #SSM: {stress_period: [lay,row,col,concentration,itype]}
                    new_ssmdata[i].append([j,row,col,Cfresh,itype['WEL']]) #since it's a sink, conc. doesn't matter
                    wel_cells.append((j,row,col))
                continue
    wel_cells = tuple(np.array(list(set(wel_cells))).T)
    return new_weldata, new_ssmdata,wel_cells

#Add recharge if desired
def make_rech_array(low=1e-2,high=1e0):
    import scipy.stats as sts
    llow,lhigh = np.log10((low,high))
    rech = np.exp(sts.uniform.rvs(size=1,loc=llow,scale=lhigh-llow)[0])
    return rech/(nrow*ncol)

def add_recharge_cells(recharge_generator,const=1,*args):
    if const==1:
        rech_data = recharge_generator(*args)
    else:
        rech_data = {}
        for i in range(nper):
            rech_array = recharge_generator(*args)
        rech_data[i] = rech_array
    return rech_data

def sample_dist(distclass,size,*args):
    smp = distclass.rvs(*args,size=size)
    if size==1:
        smp=smp[-1]
    return smp

def write_sample(fname,varname,distclass,sample):
    fout= open(fname,"a")
    fout.write(varname + ',' + str(type(distclass)) + ',' + str(sample) + '\n')
    fout.close()
    return

def truncate_grf(grid,lith_props,hk_vals,log10trans=True,plotyn=False):
    grid_cutoffs = []
    for q in np.cumsum(lith_props):
        grid_cutoffs.append(np.quantile(grid,q))

    if plotyn:
        h = plt.hist(grid.flatten())
        for cutoff in grid_cutoffs:
            plt.vlines(cutoff,0,14000)
        plt.show()

    outgrid = np.ones(grid.shape,dtype=np.float32)
    for i,cutoff in reversed(list(enumerate(grid_cutoffs))):
        outgrid[np.where(grid<cutoff)] = hk_vals[i]

    if plotyn:
        f,axs = plt.subplots(2,1,sharex=True)
        axs[0].imshow(grid[:,0,:])
        axs[1].imshow(outgrid[:,0,:])
    if log10trans:
        return np.power(10,outgrid)
    else:
        return outgrid
#%%
#Name model
modelname = 'passive_active'
tot_it = 2
# run installed version of flopy or add local path
try:
    import flopy
except:
    fpth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(fpth)
    import flopy

print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('flopy version: {}'.format(flopy.__version__))


if sys.platform == "darwin":
    repo = Path('/Users/ianpg/Documents/ProjectsLocal/DelawareSGD')
    model_ws = os.path.join('/Users/ianpg/Documents/ProjectsLocal/DelawareSGD','work',modelname)
elif sys.platform == "win32":
    repo = Path('E:\Projects\DelawareSGD')
    model_ws = os.path.join('E:\Projects\DelawareSGD','work',modelname)

sys.path.append(repo)
import SGD
import config

if not os.path.exists(model_ws):
    os.makedirs(model_ws)
sys.path.append(os.path.join(model_ws,'..','..'))

sw_exe = config.swexe #set the exe path for seawat
print('Model workspace:', os.path.abspath(model_ws))


#%%
#Model discretization
Lx = 3000.
Ly = 600.
Lz = 80.
nlay = int(Lz/3)
nrow = int(Ly/30)
ncol = int(Lx/30)

henry_top = 5
ocean_elev = 0

delv_first = Lz/nlay
botm_first = henry_top-delv_first

delv = (Lz-delv_first) / (nlay-1)
delr = Lx / ncol
delc = Ly / nrow

henry_botm = np.hstack(([botm_first],np.linspace(botm_first-delv,henry_top-Lz,nlay-1)))
delv_vec = np.hstack((delv_first,np.repeat(delv,nlay-1)))
delv_weight = [x/np.sum(delv_vec) for x in delv_vec]

beachslope = .05
ocean_col = [np.floor(ncol-1).astype('int'),ncol-1] #Manually done to make sure it's in the right place rn
inland_elev = beachslope*ocean_col[0]*delr
offshore_elev = -beachslope*(ocean_col[1]-ocean_col[0])*delr


#Period data
Lt = 360*20 #Length of time in days
perlen = list(np.repeat(180,int(Lt/180)))
nstp = list(np.ones(np.shape(perlen),dtype=int))

nper = len(perlen)
steady = [False for x in range(len(perlen))] #Never steady
itmuni = 4 #time unit 4= days
lenuni = 2 #length unit 2 = meter
tsmult = 1.8
ssm_data = None
verbose = True

print('Model setup: \n'
      'nlay: {}\n'
      'nrow: {}\n'
      'ncol: {}\n'
      'Total cells: {}\n'
      'Total time: {} days\n'
      'nper: {}\n'.format(nlay,nrow,ncol,nlay*nrow*ncol,Lt,nper))
# In[4]:

#Create basic model instance and dis pacakge
m = flopy.seawat.Seawat(modelname, exe_name=sw_exe, model_ws=model_ws,verbose=verbose)
SGD.ModelSGD.Seawat2SGD(m)  #convert to subclass ModelSGD
print(m.namefile)

# Add DIS package to the MODFLOW model
dis = flopy.modflow.ModflowDis(m, nlay, nrow, ncol, nper=nper, delr=delr,
                               delc=delc,
                               laycbd=0, top=henry_top,
                               botm=henry_botm, perlen=perlen, nstp=nstp,
                               steady=steady,itmuni=itmuni,lenuni=lenuni,
                               tsmult=tsmult)

# In[5]:

#Hydraulic conductivity field
hkSand = 10.  #horizontal hydraulic conductivity m/day
hkClay = hkSand*.01

heterogenous = True
mu = np.log(hkSand)
sill = .1
modeltype = 'Exponential'
llay = int(20/np.mean(delv))
lrow = int(2000/delc)
lcol = int(2000/delr)


if heterogenous:
    import simulationFFT
    fft_grid = np.exp(simulationFFT.simulFFT(nrow, nlay, ncol, mu, sill, modeltype, lrow , llay, lcol))
    #hk[0:int(np.where(henry_botm==find_nearest(henry_botm,ocean_elev))[0])+1,:,:] = hkSand
else:
    hk = hkSand*np.ones((nlay,nrow,ncol), dtype=np.int32)

grid = np.log10(fft_grid)
#lith_props = [0.2,0.5,0.3]
#hk_vals = [-1,0,2]
lith_props = [0.2,0.8]
hk_vals = [0,2]

log10trans = True
plotyn= True
hk = truncate_grf(grid,lith_props,hk_vals,log10trans=True,plotyn=plotyn)


#plt.figure(),plt.imshow((hk[:,0,:])),plt.colorbar(),plt.title('Sill:{}'.format(sill)),plt.show()

#Set Hydraulic properties
sy = 0.24
ss = 1e-5
por = 0.3
vka = 1 # = vk/hk
al = 1 #longitudinal dispersivity (m) from Walther et al. 2017
dmcoef = 2e-9 #m2/day

#Variable density parameters
Csalt = 35.0001
Cfresh = 0.
densesalt = 1025.
densefresh = 1000.
denseslp = (densesalt - densefresh) / (Csalt - Cfresh)
#denseslp = 0 #trick for testing constant density

# In[8]:

#BCs
bc_ocean = 'GHB'
bc_right_edge = 'GHB'
bc_inland = 'GHB'
add_wells = 0
n_wells = 0
rech_on = 0

#Inland
calc_inland_head = 0 #calculate from hgrad
manual_inland_head = 0.3184
start_fresh_yn = 1
ocean_shead = [ocean_elev for x in range(len(perlen))]
ocean_ehead = ocean_shead

# save cell fluxes to unit 53
ipakcb = 53

#Create ocean boundary at top of model
ocean_col_vec = (0,0,np.arange(ocean_col[0],ocean_col[1]+1))
ocean_coords = (0,slice(0,nrow),slice(ocean_col[0],ocean_col[1]+1)) #top of the model
ocean_bool = np.zeros((nlay,nrow,ncol), dtype=bool)
ocean_bool[0,:,np.arange(ocean_col[0],ocean_col[1]+1)] = 1
m.ocean_arr = ocean_bool


if calc_inland_head == 1:
    head_inland = ocean_col[0]*delc*hgrad + ocean_elev
else:
    head_inland = manual_inland_head

####IN TESTING#####
#Create a line of where the ocean is, and any nodes on right edge below ocean
offshore_lay = (np.abs(henry_botm-offshore_elev)).argmin().astype('int')
if ocean_col[0] == ncol-1:
    ocean_line = []
    bc_ocean = 'XXX'
else:
    ocean_line = get_line((0,ocean_col[0]),(offshore_lay,ocean_col[1]),allrows=1,nrow=1)

ocean_line_tuple = tuple(np.array(ocean_line).T) #use this for indexing numpy arrays
right_edge = get_ocean_right_edge(m,ocean_line_tuple,
                                  int(np.where(henry_botm==find_nearest(henry_botm,ocean_elev))[0]))
left_edge = get_ocean_right_edge(m,ocean_line_tuple,
                                  int(np.where(henry_botm==find_nearest(henry_botm,head_inland))[0]),
                                col=0)

#Create ibound
ibound,ocean_hf = shade_above(nlay,nrow,ncol,ocean_line) #don't set ibound of ocean
#ibound[:right_edge[0][0],right_edge[1][0],right_edge[2][0]] = 0
#ibound[:right_edge[0][0],right_edge[1][0],0] = 0

if bc_ocean == 'GHB':
    ibound[ocean_line_tuple]=1


#Set starting heads
strt = np.zeros((nlay,nrow,ncol),dtype=np.int)


#Transport BCs
if start_fresh_yn == 1:
    sconc = Cfresh*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin fresh
elif start_fresh_yn == 0:
    sconc = Csalt*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin SW-saturated
else:
    sconc = Cfresh*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin SW-saturated
    sconc[:,:,int(np.floor(ncol/2)):-1] = Csalt

if ocean_hf:
    sconc[ocean_hf[0:3]] = Csalt
sconc[right_edge] = Csalt
sconc[:,:,0] = Cfresh

icbund = np.ones((nlay, nrow, ncol), dtype=np.int)
icbund[np.where(ibound==-1)] = -1

def make_bc_dicts():
    #Ocean and inland boundary types
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    chd_data = {}
    ssm_data = {}
    ghb_data = {}
    wel_data = {}
    for i in range(nper):
        dat_chd = []
        dat_ssm = []
        dat_ghb = []
        dat_wel = []
        #Ocean boundary
        if ocean_hf:
            for j in range(np.size(ocean_hf[0])):
                if bc_ocean=='CHD':
                    #CHD: {stress_period: [lay,row,col,starthead,endhead]}
                    dat_chd.append([ocean_line_tuple[0][j],
                                ocean_line_tuple[1][j],
                                ocean_line_tuple[2][j],
                                ocean_shead[i],
                                ocean_ehead[i]])
                    #SSM: {stress_period: [lay,row,col,concentration,itype]}
                    dat_ssm.append([ocean_line_tuple[0][j],
                                ocean_line_tuple[1][j],
                                ocean_line_tuple[2][j],
                                Csalt,
                                itype['CHD']])
                elif bc_ocean=='GHB':
                    #GHB: {stress period: [lay,row,col,head level,conductance]}
                    #conductance c = K*A/dL; assume horizontal flow at outlet,
                    #and calculate length to be at edge of ocean cell, as opposed to mipoint
                    # c = (K*dy*dz)/(dx/2) = 2*K*delr*delv/delc
                    dat_ghb.append([ocean_hf[0][j],
                                   ocean_hf[1][j],
                                   ocean_hf[2][j],
                                   #ocean_hf[3][j],
                                    ocean_elev,
                                   2*hkSand*delc*delv_vec[ocean_hf[0][j]]/delr])
                    #SSM: {stress_period: [lay,row,col,concentration,itype]}
                    dat_ssm.append([ocean_hf[0][j],
                                   ocean_hf[1][j],
                                   ocean_hf[2][j],
                                   Csalt,
                                   itype['GHB']])
        else:
            pass
        #Right edge boundary
        if bc_right_edge=='GHB':
            for j in range(np.size(right_edge[0])):
                #GHB: {stress period: [lay,row,col,head level,conductance]}
                #conductance c = K*A/dL; assume horizontal flow at outlet,
                #and calculate length to be at edge of ocean cell, as opposed to mipoint
                # c = (K*dy*dz)/(dx/2) = 2*K*delr*delv/delc
                dat_ghb.append([right_edge[0][j],
                               right_edge[1][j],
                               right_edge[2][j],
                               #ocean_hf[3][j],
                                ocean_elev,
                               2*hkSand*delc*delv_vec[right_edge[0][j]]/delr])
                #SSM: {stress_period: [lay,row,col,concentration,itype]}
                dat_ssm.append([right_edge[0][j],
                               right_edge[1][j],
                               right_edge[2][j],
                               Csalt,
                               itype['GHB']])
        else:
            pass
        #Inland boundary
        if bc_inland=='GHB':
            for j in range(np.size(left_edge[0])):
                dat_ghb.append([left_edge[0][j],
                               left_edge[1][j],
                               left_edge[2][j],
                                head_inland,
                               2*hkSand*delc*delv_vec[left_edge[0][j]]/delr])
                #SSM: {stress_period: [lay,row,col,concentration,itype]}
                dat_ssm.append([left_edge[0][j],
                               left_edge[1][j],
                               left_edge[2][j],
                               Cfresh,
                               itype['GHB']])
        elif bc_inland=='WEL':
            for j in range(nlay):
                for k in range(nrow):
                    #WEL: {stress_period: [lay,row,col,flux]}
                    dat_wel.append([j,k,0,influx*delv_weight[j]/nrow])
                    #SSM: {stress_period: [lay,row,col,concentration,itype]}
                    dat_ssm.append([j,k,0,Cfresh,itype['WEL']])
        chd_data[i] = dat_chd
        ssm_data[i] = dat_ssm
        ghb_data[i] = dat_ghb
        wel_data[i] = dat_wel

    #saving concentrations at specified times
    #timprs = [k for k in range(1,np.sum(perlen),50)]
    return chd_data, ssm_data, ghb_data, wel_data

chd_data, ssm_data, ghb_data, wel_data = make_bc_dicts()
wel_data_base,ssm_data_base = wel_data,ssm_data
timprs = np.round(np.linspace(1,np.sum(perlen),20),decimals=0)

create_MC_file()
save_obj(m.MC_file.parent,wel_data_base,'wel_data_base')
save_obj(m.MC_file.parent,ssm_data_base,'ssm_data_base')

# In[9]:

#### ADD WELL AND RECHARRGE DATA####
#Winter is even stress periods, summer is odd SP.
#Winter= wells OFF, natural precip (rech) ON, irrigation rech OFF,
#Summer = wells ON, irrigation rech (farm_rech) ON,  precip (rech) OFF
kper_odd = list(np.arange(1,nper,2))
kper_even = list(np.arange(0,nper,2))

##Add recharge data
rech = 1e-6

#Assign the location of the farms
farm_leftmargin = 10
farm_uppermargin = 1
nfarms = 4
farm_size = (200,200) #m in row,col direction
farm_size_rowcol = (int(farm_size[0]/delc),int(farm_size[1]/delr)) #size of farm in number of row,col

farm_loc_r = []
farm_loc_c = []
farm_orig = []
for x in range(int(nfarms/2)):
    for y in range(2):
        for z1 in range(farm_size_rowcol[0]):
            for z2 in range(farm_size_rowcol[1]):
                farm_loc_r.append(farm_uppermargin + y*(farm_size_rowcol[0]+2) + z1)
                farm_loc_c.append(farm_leftmargin + x*(farm_size_rowcol[1]+2) + z2)
                if (z1==0) and (z2==0):
                    farm_orig.append((farm_loc_r[-1],farm_loc_c[-1])) #upper left of ea. farm=loc of well
farm_loc = (np.array(farm_loc_r),np.array(farm_loc_c))

## Add well data
n_wells = nfarms
wel_flux = list(np.ones(n_wells)*10)
top_lay = int(np.where(henry_botm==find_nearest(henry_botm,ocean_elev))[0])+1
wel_data,ssm_data,wel_cells = add_pumping_wells(wel_data_base.copy(),ssm_data_base.copy(),n_wells,wel_flux,farm_orig,kper_odd)


hk[wel_cells] = hkSand


## Add farm recharge data
farm_rech_flux = wel_flux[0]*0.2
farm_rech = np.zeros((nrow,ncol),dtype=np.float)
farm_rech[farm_loc] = farm_rech_flux/np.prod(farm_size)
#Set rech_data for winter and summer
rech_data = {}
for i in range(len(perlen)):
    if i in kper_even:
        rech_data[i] = rech
    elif i in kper_odd:
        rech_data[i] = farm_rech

# In[10]:

riv_loc = get_line((0,0),(0,ncol-1),allrows=1,nrow=nrow)
riv_loc = [x for x in riv_loc if x[1]==int(nrow/2)]
riv_loc = tuple(np.array(riv_loc).T)

riv_grad = .001
rbot_vec = np.linspace(riv_grad*Lx,ocean_elev,ncol)

#Stage and conductance:
stage = 1
cond = 10
riv_grad = .001

def write_river_data(riv_loc,stage,cond,riv_grad,kper,ssm_data):

    ####ADD RIVER DATA####
    rbot_vec = np.linspace(riv_grad*Lx,ocean_elev,ncol)

    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    riv_data = {}
    new_ssm_data = ssm_data
    for i in range(nper):
        dat_riv = []
        if i in kper:
            for j in range(np.size(riv_loc[0])):
                #RIV: {stress_period:[lay, row, col, stage, cond, rbot],...}
                dat_riv.append([riv_loc[0][j],
                                    riv_loc[1][j],
                                    riv_loc[2][j],
                                    stage+rbot_vec[riv_loc[2][j]],
                                    cond,
                                    rbot_vec[riv_loc[2][j]]])
                #SSM: {stress_period: [lay,row,col,concentration,itype]}
                new_ssm_data[i].append([riv_loc[0][j],
                                        riv_loc[1][j],
                                        riv_loc[2][j],
                                        Cfresh,
                                        itype['RIV']])
        else:
            for j in range(np.size(riv_loc[0])):
                #RIV: {stress_period:[lay, row, col, stage, cond, rbot],...}
                dat_riv.append([riv_loc[0][j],
                                    riv_loc[1][j],
                                    riv_loc[2][j],
                                    rbot_vec[riv_loc[2][j]], #set stage as bottom of river
                                    cond,
                                    rbot_vec[riv_loc[2][j]]])
                #SSM: {stress_period: [lay,row,col,concentration,itype]}
                new_ssm_data[i].append([riv_loc[0][j],
                                        riv_loc[1][j],
                                        riv_loc[2][j],
                                        Cfresh,
                                        itype['RIV']])
        riv_data[i] = dat_riv
    return riv_data,new_ssm_data

riv_data,ssm_data = write_river_data(riv_loc,stage,cond,riv_grad,kper_even,ssm_data)

# In[9]:

#Output control
oc_data = {}
for kper in range(nper):
    oc_data[(kper,0)] = ['save head','save budget']


# In[10]:

#Create instances in flopy
bas = flopy.modflow.ModflowBas(m, ibound, strt=strt)
if bc_ocean=='CHD' or bc_inland=='CHD' :
    chd = flopy.modflow.ModflowChd(m, stress_period_data=chd_data)
if bc_ocean=='GHB' or bc_inland=='GHB'or bc_right_edge=='GHB':
    ghb = flopy.modflow.ModflowGhb(m, stress_period_data=ghb_data)

rch = flopy.modflow.ModflowRch(m, rech=rech_data)
wel = flopy.modflow.ModflowWel(m, stress_period_data=wel_data, ipakcb=ipakcb)
riv = flopy.modflow.ModflowRiv(m, stress_period_data=riv_data)
# Add LPF package to the MODFLOW model
lpf = flopy.modflow.ModflowLpf(m, hk=hk, vka=vka, ipakcb=ipakcb,laytyp=1,laywet=1,
                              ss=ss,sy=sy)

# Add PCG Package to the MODFLOW model
pcg = flopy.modflow.ModflowPcg(m, hclose=1.e-8)

# Add OC package to the MODFLOW model
oc = flopy.modflow.ModflowOc(m,
                             stress_period_data=oc_data,
                             compact=True)

#Create the basic MT3DMS model structure
btn = flopy.mt3d.Mt3dBtn(m,
                         laycon=lpf.laytyp, htop=henry_top,
                         dz=dis.thickness.get_value(), prsity=por, icbund=icbund,
                         sconc=sconc, nprs=1,timprs=timprs)
adv = flopy.mt3d.Mt3dAdv(m, mixelm=-1)
dsp = flopy.mt3d.Mt3dDsp(m, al=al, dmcoef=dmcoef)
gcg = flopy.mt3d.Mt3dGcg(m, iter1=50, mxiter=1, isolve=1, cclose=1e-5)
ssm = flopy.mt3d.Mt3dSsm(m, stress_period_data=ssm_data)

#vdf = flopy.seawat.SeawatVdf(m, iwtable=0, densemin=0, densemax=0,denseref=1000., denseslp=0.7143, firstdt=1e-3)
vdf = flopy.seawat.SeawatVdf(m, mtdnconc=1, mfnadvfd=1, nswtcpl=0, iwtable=1,
                             densemin=0., densemax=0., denseslp=denseslp, denseref=densefresh)

# In[11]:

printyn = 0
gridon=0
rowslice=0
rowslice=farm_orig[0][0]
m.plot_hk_ibound(rowslice=rowslice,printyn=printyn,gridon=gridon);


# In[12]:

#Write input
m.write_input()

# Try to delete the output files, to prevent accidental use of older files
try:
    os.remove(os.path.join(model_ws,'MT3D.CNF'))
    os.remove(os.path.join(model_ws,'MT3D001.MAS'))
    os.remove(os.path.join(model_ws, 'MT3D001.UCN'))
    os.remove(os.path.join(model_ws, modelname + '.hds'))
    os.remove(os.path.join(model_ws, modelname + '.cbc'))
except:
    pass


# In[13]:

#Run model
import datetime
ts = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M')
v = m.run_model(silent=False, report=True)
for idx in range(-3, 0):
    print(v[1][idx])

# In[14]:

#Post-processing functions
def plotdischarge(modelname,model_ws,color='w',per=-1,scale=50,rowslice=0):
    fname = os.path.join(model_ws, '' + modelname + '.cbc')
    budobj = flopy.utils.CellBudgetFile(fname)
    qx = budobj.get_data(text='FLOW RIGHT FACE')[per]
    qz = budobj.get_data(text='FLOW LOWER FACE')[per]

    # Average flows to cell centers
    qx_avg = np.empty(qx.shape, dtype=qx.dtype)
    qx_avg[:, :, 1:] = 0.5 * (qx[:, :, 0:ncol-1] + qx[:, :, 1:ncol])
    qx_avg[:, :, 0] = 0.5 * qx[:, :, 0]
    qz_avg = np.empty(qz.shape, dtype=qz.dtype)
    qz_avg[1:, :, :] = 0.5 * (qz[0:nlay-1, :, :] + qz[1:nlay, :, :])
    qz_avg[0, :, :] = 0.5 * qz[0, :, :]

    y, x, z = dis.get_node_coordinates()
    X, Z = np.meshgrid(x, z[:, 0, 0])
    iskip = 1 #how many cells to skip, 1 means plot every cell

    ax = plt.gca()
    cpatchcollection = ax.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip],
              qx_avg[::iskip, rowslice, ::iskip], -qz_avg[::iskip, rowslice, ::iskip],
              color=color, scale=scale, headwidth=4, headlength=2,
              headaxislength=1, width=0.0025)
    return cpatchcollection

def permute_kstpkper(ucnobj):
    kstpkper = ucnobj.get_kstpkper()
    kstpkper_unique = []
    index_unique = []
    niter = 0
    for entry in kstpkper:
        if not entry in kstpkper_unique:
            kstpkper_unique.append(entry)
            index_unique.append(niter)
        niter += 1
    return kstpkper_unique, index_unique

def kstpkper_from_time(ucnobj,tottim):
    kstpkpers = ucnobj.get_kstpkper()
    times = ucnobj.get_times()
    timeind = times.index(tottim)
    kstpkper = kstpkpers[timeind]
    return kstpkper

def kstpkper_ind_from_kstpkper(ucnobj,kstpkper=(0,0)):
    kstpkpers = ucnobj.get_kstpkper()
    kstpkper_unique = permute_kstpkper(ucnobj)[0]
    kstpkper_ind = kstpkper_unique.index(kstpkper)
    return kstpkper_ind

def get_salt_outflow(m,kstpkper=None,totim=None):
    fname = os.path.join(m.model_ws, 'MT3D001.UCN')
    ucnobj = flopy.utils.binaryfile.UcnFile(fname)
    totim = ucnobj.get_times()[-1]
    if kstpkper==None:
        kstpkper = ucnobj.get_kstpkper()[-1]
    ocean_conc = ucnobj.get_data(kstpkper=kstpkper)
    return ocean_conc

def plot_background(mm,array,label=None):
    if label==None:
        label = [k for k,v in globals().items() if v is array][-1]
    if label=='hk':
        norm=matplotlib.colors.LogNorm()
        vmin=hkClay
        vmax=hkSand
        cmap='jet'
    else:
        norm = None
        vmin=None
        vmax=None
        cmap='jet'
    cpatchcollection = mm.plot_array(array,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    cpatchcollection.set_label(label)
    return cpatchcollection,label

def plot_mas(m):
    # Load the mas file and make a plot of total mass in aquifer versus time
    fname = os.path.join(m.model_ws, 'MT3D001.MAS')
    mas = flopy.mt3d.Mt3dms.load_mas(fname)
    f = plt.figure()
    ax = f.add_subplot(1, 1, 1)
    plt.xlabel('Time (d)')
    plt.ylabel('Mass (kg)')
    plt.title('Mass of salt within model through time')
    lines = ax.plot(mas.time, mas.total_mass)
    plt.show()
    return mas

def extract_hds_conc(per):
    fname = os.path.join(model_ws, '' + modelname + '.hds')
    hdobj = flopy.utils.binaryfile.HeadFile(fname)
    times = hdobj.get_times()
    hds = hdobj.get_data(totim=times[per])
    hds[np.where(ibound != 1)] = np.nan
    hds[np.where((hds>1e10) | (hds<-1e10))] = np.nan

    # Extract final timestep salinity
    fname = os.path.join(model_ws, 'MT3D001.UCN')
    ucnobj = flopy.utils.binaryfile.UcnFile(fname)
    times = ucnobj.get_times()
    kstpkper = ucnobj.get_kstpkper()
    conc = ucnobj.get_data(totim=times[per])
    conc[np.where(ibound != 1)] = np.nan
    conc[np.where((conc>1e10) | (conc<-10))] = np.nan
    return conc,hds

# Make head and quiver plot
import utils
def basic_plot(per,backgroundplot,rowslice=0,printyn=0,contoursyn=1,**kwargs):
    printyn = 1

    f, axs = plt.subplots(1, figsize=(6, 2))

    plt.tight_layout()

    #Plot discharge and ibound
    mm = flopy.plot.ModelCrossSection(ax=axs, model=m, line={'row':rowslice})

    #Plot background
    backgroundpatch,lbl = cpatchcollection,label = plot_background(mm,backgroundplot,'conc(g/L)')
    lvls = Cfresh + (Csalt-Cfresh)*np.array([.05,.5,.95])
    if contoursyn==1:
        CS = mm.contour_array(conc,head=hds,levels=lvls,colors='white')
        plt.clabel(CS, CS.levels, inline=True, fontsize=10)

    #mm.contour_array(hds,head=hds)
    mm.plot_ibound()
    mm.plot_bc(ftype='GHB',color='blue')
    if m.Wel:
        mm.plot_bc(ftype='WEL',color='black')
    #Plot discharge
    utils.plotdischarge(m,color='white',per=per,rowslice=rowslice,**kwargs);
    plt.xlabel('Distance (m)')
    plt.ylabel('Elevation (m)')
    plt.title('Flow during period {} of {}'.format(np.arange(nper)[per],nper-1))
    plt.subplots_adjust(bottom=.1)

    #align plots and set colorbar
    f.subplots_adjust(left=.1,right=0.88)
    cbar_ax = f.add_axes([0.90, 0.1, 0.02, 0.7])
    cb = f.colorbar(cpatchcollection,cax=cbar_ax)
    cb.set_label(label)
    if printyn == 1:
        plt.savefig(os.path.join(m.model_ws, m.name + '_' + ts + '_flowvec_row' + str(rowslice) +
                                 '_per' + str(per) + '_' + lbl[:3] + '.png'),dpi=150)
    plt.show()
    return


# In[19]:
per = [-1,-2]
mas = plot_mas(m)
rowslice = 1
for p in per:
    conc,hds = extract_hds_conc(p)
    basic_plot(p,hk,rowslice=rowslice,scale=70,iskip=3,printyn=1,contoursyn=1)
m.plot_hk_ibound(rowslice=rowslice,gridon=0)
# In[21]:

def add_to_paramdict(paramdict,paramname,val):
    if paramdict is None:
        paramdict = {}
    if  paramname in list(paramdict.keys()):
        paramdict[paramname].append(val)
    else:
        #paramdict.update(paramname=[val])
        paramdict[paramname] = [val]
    return


def record_salinity(m,totim=None,fname_write=None,ts_hms=None):
    from pathlib import Path
    if ts_hms is None:
        ts_hms = datetime.datetime.now().strftime('%H-%M-%S')
    # Extract final timestep salinity
    fname = os.path.join(m.model_ws, 'MT3D001.UCN')
    ucnobj = flopy.utils.binaryfile.UcnFile(fname)
    if totim is None:
        totim = ucnobj.get_times()[-1]
    conc = ucnobj.get_data(totim=totim)
    if fname_write is None:
        fname_write = m.MC_file.parent.joinpath('conc_' + str(int(totim)) + '_' + ts_hms + '.npy')
        print(fname_write)
        np.save(fname_write,conc)
    return conc

def copy_rename(src_file, dst_file):
    import shutil
    from pathlib import Path
    shutil.copy(str(Path(src_file)),str(Path(dst_file)))
    return


# In[22]:

def get_yn_response(prompt):
    while True:
        try:
            resp = str(input(prompt))
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        if resp[0] is 'y':
            value = True
            break
        elif resp[0] is 'n':
            value = False
            break
        else:
            print('This didnt work right. Try again')
            continue
    return value

def get_value(prompt):
    while True:
        try:
            resp = str(input(prompt))
            break
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
    return resp

def check_MC_inputParams():
    if m.MC_file is not None:
        use_existing_MCfile = get_yn_response("m.MC_file already exists, continue using this experiment?")
    else:
        use_existing_MCfile = False
    if use_existing_MCfile:
        if m.inputParams is not None:
            if len(m.inputParams)>0:
                add_to_inputParams = get_yn_response("m.inputParams already has entries, do you want to add to it?")
            else:
                add_to_inputParams =False
            if add_to_inputParams:
                pass
            else:
                m.inputParams = {}
        else:
            m.inputParams = {}
    else:
        load_existing_MCfile = get_yn_response("load MC file?")
        if load_existing_MCfile:
            f = get_value("path to MC_file (path/to/test.expt): ")
            m.inputParams = load_obj(Path(f),'inputParams')
            print('loaded .pkl file!')
        else:
            create_MC_file()
            m.inputParams = {}
    return

def idx2centroid(node_coord_tuple,idx_tuple):
    z_pt = node_coord_tuple[2][idx_tuple]
    x_pt = node_coord_tuple[1][idx_tuple[2]]
    y_pt = node_coord_tuple[0][idx_tuple[1]]
    return (z_pt,y_pt,x_pt)


# ### Run the MC experiment:

# In[23]:

from scipy.io import savemat,loadmat

def run_MC(tot_it):
    #### MAKE NEW/ADD TO OLD EXPT ####
    check_MC_inputParams()

    #### VARY PARAMS ####
    it = 0
    while it < tot_it:
        ssm_data = {}
        it += 1
        
        
        ''' HOMOGENOUS ONLY
        ##hk: hk
        #      Uniform (1,100)
        low= 1e-1
        high = 100
        parname='hk'
        hk = sample_dist(sts.uniform,1,1,m,'hk',0,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,hk)
        '''
#########HETEROGENOUS ONLY ##############
        #hk1
        low= -3
        high = 0
        parname='hk1'
        val = sample_dist(sts.uniform,1,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,val)
        hk1 = 10**val
        
        #hk2
        low= 1
        high = 3
        parname='hk2'
        val = sample_dist(sts.uniform,1,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,val)
        hk2 = 10**val
        
        #lith_prop
        low= 0
        high = 0.5
        parname='lith_prop'
        val = sample_dist(sts.uniform,1,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,val)
        lith_prop = val
        
        #vario_type
        parname='vario_type'
        val = int(round(np.random.rand()))
        add_to_paramdict(m.inputParams,parname,val)
        if val==1:
            vario_type = 'Gaussian'
        elif val==0:
            vario_type = 'Exponential'
        else:
            pass

        #corr_len
        low= 250
        high = 1000
        parname='corr_len'
        val = sample_dist(sts.uniform,1,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,val)
        corr_len = val
        
        #corr_len_zx
        # equal to lz/lx
        low= .01
        high = .1
        parname='corr_len_zx'
        val = sample_dist(sts.uniform,1,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,val)
        corr_len_zx = val
        
        #corr_len_yx
        # equal to ly/lx
        low= 0.1
        high = 1
        parname='corr_len_yx'
        val = sample_dist(sts.uniform,1,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,val)
        corr_len_yx = val

        
        #Create hk grid
        mu = np.log(hkSand)
        sill = 1
        lcol = int(corr_len/delr)
        llay = int(corr_len*corr_len_zx/np.mean(delv))
        lrow = int(corr_len*corr_len_yx/delc)
        fft_grid = np.exp(simulationFFT.simulFFT(nrow, nlay, ncol, mu, sill, vario_type, lrow , llay, lcol))
        grid = np.log10(fft_grid)
        lith_props = [lith_prop,1-lith_prop]
        hk_vals = [hk1,hk2]
        plotyn= True
        hk = truncate_grf(grid,lith_props,hk_vals,log10trans=False,plotyn=plotyn)


######## END OF HETEROGENOUS BLOCK ############

''' NOTE TO SELF: NEED TO UPDATE sample_dist() arguments below since changing the function'''

        ##vka: ratio of vk/hk
        #      Uniform (1/20,1)
        low= 1/20
        high = 1
        parname='vka'
        vka = sample_dist(sts.uniform,1,1,m,'vka',0,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,vka)

        ##al: #longitudinal dispersivity (m)
        #      Uniform [0.1,20] #10 from Walther et al
        low= 0.1
        high = 20
        parname='al'
        al = sample_dist(sts.uniform,1,1,m,'al',0,*(low,high-low))
        add_to_paramdict(m.inputParams,parname,al)

        ##dmcoef: #dispersion coefficient (m2/day)
        #      log-uniform [1e-10,1e-5] #2e-9 from Walther et al
        lowhigh = np.log10([1e-10,1e-5])
        parname='dmcoef'
        dmcoef = sample_dist(sts.uniform,1,1,m,'dmcoef',1,*(lowhigh[0],lowhigh[1]-lowhigh[0]))
        add_to_paramdict(m.inputParams,parname,dmcoef)

        ##rech
        rechargs = tuple(np.log10((1e-6/(nrow*ncol),1e-1/(nrow*ncol))))
        rechargs = (rechargs[0],rechargs[1]-rechargs[0])
        rech = sample_dist(sts.uniform,1,0,m,'rech',1,*rechargs)

        farm_rechargs = (rechargs[0],rechargs[1]) #note: this is in log space
        farm_rech_flux = sample_dist(sts.uniform,1,0,m,'rech',1,*farm_rechargs)
        farm_rech[farm_loc] = farm_rech_flux
        rech_data = {}
        for i in range(len(perlen)):
            if i%2==0:
                rech_data[i] = np.ones((nrow,ncol),dtype=np.float)*rech
            else:
                rech_data[i] = farm_rech
        parname = 'rech'
        add_to_paramdict(m.inputParams,parname,rech)
        parname = 'farm_rech'
        add_to_paramdict(m.inputParams,parname,farm_rech_flux)

        ##wel
        lowhigh = np.log10((1e1,1e3))
        wel_flux = sample_dist(sts.uniform,n_wells,0,m,'wel',1,*(lowhigh[0],lowhigh[1]-lowhigh[0]))
        parname = 'wel'
        for i in range(n_wells):
            parname_temp = parname+str(i)
            add_to_paramdict(m.inputParams,parname_temp,wel_flux[i])

        #write wel data
        ssm_data_base = load_obj(m.MC_file.parent,'ssm_data_base')
        wel_data_base = load_obj(m.MC_file.parent,'wel_data_base')
        wel_data,ssm_data,_ = add_pumping_wells(wel_data_base,ssm_data_base,n_wells,wel_flux,farm_orig,kper_odd)

        ##riv_stg
        parname = 'riv_stg'
        lowhigh = (.5,1.5)
        stage = sample_dist(sts.uniform,1,0,m,'riv_stg',0,*(lowhigh[0],lowhigh[1]-lowhigh[0]))
        add_to_paramdict(m.inputParams,parname,stage)

        ##riv_cond
        parname = 'riv_cond'
        lowhigh = np.log10((.1,100))
        cond = sample_dist(sts.uniform,1,0,m,'riv_cond',1,*(lowhigh[0],lowhigh[1]-lowhigh[0]))
        add_to_paramdict(m.inputParams,parname,cond)

        #Write river data--take SSM data from WEL!!
        riv_grad = .0005
        riv_data,ssm_data = write_river_data(riv_loc,stage,cond,riv_grad,kper_even,ssm_data)

        #Create instances in flopy
        bas = flopy.modflow.ModflowBas(m, ibound, strt=strt)
        if bc_ocean=='CHD' or bc_inland=='CHD' :
            chd = flopy.modflow.ModflowChd(m, stress_period_data=chd_data)
        if bc_ocean=='GHB' or bc_inland=='GHB'or bc_right_edge=='GHB':
            ghb = flopy.modflow.ModflowGhb(m, stress_period_data=ghb_data)

        rch = flopy.modflow.ModflowRch(m, rech=rech_data)
        wel = flopy.modflow.ModflowWel(m, stress_period_data=wel_data, ipakcb=ipakcb)
        riv = flopy.modflow.ModflowRiv(m, stress_period_data=riv_data)
        # Add LPF package to the MODFLOW model
        lpf = flopy.modflow.ModflowLpf(m, hk=hk, vka=vka, ipakcb=ipakcb,laytyp=1,laywet=1,
                                      ss=ss,sy=sy)

        # Add PCG Package to the MODFLOW model
        pcg = flopy.modflow.ModflowPcg(m, hclose=1.e-8)

        # Add OC package to the MODFLOW model
        oc = flopy.modflow.ModflowOc(m,
                                     stress_period_data=oc_data,
                                     compact=True)

        #Create the basic MT3DMS model structure
        btn = flopy.mt3d.Mt3dBtn(m,
                                 laycon=lpf.laytyp, htop=henry_top,
                                 dz=dis.thickness.get_value(), prsity=por, icbund=icbund,
                                 sconc=sconc, nprs=1,timprs=timprs)
        adv = flopy.mt3d.Mt3dAdv(m, mixelm=-1)
        dsp = flopy.mt3d.Mt3dDsp(m, al=al, dmcoef=dmcoef)
        gcg = flopy.mt3d.Mt3dGcg(m, iter1=50, mxiter=1, isolve=1, cclose=1e-5)
        ssm = flopy.mt3d.Mt3dSsm(m, stress_period_data=ssm_data)

        #vdf = flopy.seawat.SeawatVdf(m, iwtable=0, densemin=0, densemax=0,denseref=1000., denseslp=0.7143, firstdt=1e-3)
        vdf = flopy.seawat.SeawatVdf(m, mtdnconc=1, mfnadvfd=1, nswtcpl=0, iwtable=1,
                                     densemin=0., densemax=0., denseslp=denseslp, denseref=densefresh)

        #Write input
        m.write_input()

        # Try to delete the output files, to prevent accidental use of older files
        flist = [os.path.join(model_ws,'MT3D.CNF'),
                  os.path.join(model_ws,'MT3D001.MAS'),
                  os.path.join(model_ws, modelname + '.hds'),
                  os.path.join(model_ws, 'MT3D001.UCN'),
                  os.path.join(model_ws, 'MT3D001.UCN'),
                  os.path.join(model_ws, modelname + '.cbc')]
        for f in flist:
            try:
                os.remove(f)
            except:
                pass

        #Make timestamp
        import datetime
        sep = '-'
        ts = datetime.datetime.now().strftime('%m'+sep+'%d'+sep+'%H'+sep+'%M'+sep+'%S')
        ts_hms = ts.split(sep)[2:]
        ts_hms = sep.join(ts_hms)

        #Run model
        print('Running iteration {} of {}...'.format(it,tot_it))
        v = m.run_model(silent=True, report=True)
        for idx in range(-3, 0):
            print(v[1][idx])

        #Record final salinity as .npy, also move full CBC and UCN files to expt folder
        _ = record_salinity(m,ts_hms=ts_hms);
        copy_rename(os.path.join(m.model_ws,'MT3D001.UCN'),m.MC_file.parent.joinpath('conc_'+ts_hms+'.UCN').as_posix())
        #copy_rename(os.path.join(m.model_ws,m.name+'.cbc'),m.MC_file.parent.joinpath('cbc_'+ts_hms+'.cbc').as_posix())

        print('Finished iteration {} of {}'.format(it,tot_it))
    #Save inputParams immediately to prevent accidental destruction of them
    savemat(m.MC_file.parent.joinpath('inputParams.mat').as_posix(),m.inputParams)
    np.save(m.MC_file.parent.joinpath('inputParams.npy'),m.inputParams)
    save_obj(m.MC_file.parent,m.inputParams,'inputParams')
    save_obj(m.MC_file.parent,m.dis.get_node_coordinates(),'yxz')
    return m.inputParams


# In[24]:

####Run the MC experiment ####
inputParams = run_MC(tot_it)



#%% Calculate hausdorff matrix and export
import hausdorff_from_dir
#importlib.reload(hausdorff_from_dir)
fnames = hausdorff_from_dir.create_concmat_from_ucndir(m.MC_file.parent)
yxz = load_obj(m.MC_file.parent,'yxz')
export_ind = -1 #which of the saved fnames to use for conc_mat
hausdorff_from_dir.compute_export_hausdorff(m.MC_file.parent,
                                            conc_mat=np.load(m.MC_file.parent.joinpath(fnames[export_ind])),
                                            yxz=yxz)

#%% Load from Matlab
medoids = loadmat(dirname.joinpath('medoids.mat').as_posix())['medoids'][0]
mdict = {}
for i in range(len(medoids)):
    mdict['medoid_mat'+str(i)] = conc_mat[i]
savemat(dirname.joinpath('medoid_mats.mat').as_posix(),mdict)
