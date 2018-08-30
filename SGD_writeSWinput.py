import os
from pathlib import Path
import sys
import numpy as np
from __main__ import *
import flopy

base_dir = os.path.abspath(os.getcwd())
model_ws = Path(os.path.join(base_dir,'work',modelname)).as_posix()
sys.path.append(base_dir)
import SGD

def write_swt_input(modelname):
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)
    sys.path.append(os.path.join(model_ws,'..','..'))
    import config
    sw_exe = config.swexe

    os.listdir(model_ws)
    print('Model workspace:',os.path.abspath(model_ws))

    #nearest value in array
    def find_nearest(array,value):
        import numpy as np
        idx = (np.abs(array-value)).argmin()
        idx.astype('int')
        return array[idx]

    def loc_to_col(locs):
        cols = [int(find_nearest((np.arange(ncol)*delc),loc)) for loc in locs]
        return cols

    #Model grid
    #Grid: 0.5 *1*1m • Size:70m*20m
    Lx = 70.
    Lz = 20.
    Ly = 10
    nlay = 40
    nrow = 3
    ncol = 70
    delr = Lx / ncol
    delc = Ly / nrow
    delv = Lz / nlay
    henry_top = Lz
    henry_botm = np.linspace(henry_top - delv, 0., nlay)

    #Period data
    perlen = [1,100]
    nstp = [100,100]
    nper = len(perlen)
    steady = [True,False]
    itmuni = 4 #time unit 4= days
    lenuni = 2 #length unit 2 = meter
    tsmult = 1
    ssm_data = None
    verbose = True


    # In[4]:

    #Create basic model instance and dis pacakge
    m = flopy.seawat.Seawat(modelname, exe_name=config.swexe, model_ws=model_ws,verbose=verbose)
    SGD.ModelSGD.Seawat2SGD(m)
    print(type(m))
    print('Name file: ' + m.namefile)

    # Add DIS package to the MODFLOW model
    dis = flopy.modflow.ModflowDis(m, nlay, nrow, ncol, nper=nper, delr=delr,
                                   delc=delc, laycbd=0, top=henry_top,
                                   botm=henry_botm, perlen=perlen, nstp=nstp,
                                   steady=steady,itmuni=itmuni,lenuni=lenuni,
                                   tsmult=tsmult)


    # In[5]:

    #Hydraulic conductivity field 

    hkSand = 80.  #horizontal hydraulic conductivity m/day
    hkClay = 1. 
    lithmat = hkSand*np.ones((nlay,nrow,ncol), dtype=np.int32) #sandy background

    def rand_clay_blocks(lithmat,hkClay,numblocks,sizeblocks):
        nlay,nrow,ncol = lithmat.shape
        lay_block = np.random.randint(1,nlay-sizeblocks[1],numblocks)
        col_block = np.random.randint(1,ncol-sizeblocks[0],numblocks)
        lithmat_blocks = lithmat.copy()
        for i in range(numblocks):
            block_coords = [slice(lay_block[i],lay_block[i]+sizeblocks[1]),
                            0,
                            slice(col_block[i],col_block[i]+sizeblocks[0])]
            lithmat_blocks[block_coords] = hkClay
        return lithmat_blocks

    #add low conductivity regions
    lithmat = rand_clay_blocks(lithmat,hkClay,100,(5,2))
    low_k_loc = (20,30)
    low_k_col = loc_to_col(low_k_loc)


    #Set Hydraulic properties

    hk = lithmat
    sy = 0.15
    ss = 0.00005
    por = 0.2
    vka = 1/10 # = vk/hk
    al = 3 #longitudinal dispersivity, in m
    #dmcoef = 0.57024 #m2/day  Could also try 1.62925 as another case of the Henry problem
    dmcoef = 0. #test for numerical dispersion

    #Variable density parameters
    Csalt = 35.0001
    Cfresh = 0.
    densesalt = 1025.
    densefresh = 1000.
    denseslp = (densesalt - densefresh) / (Csalt - Cfresh)


    # In[7]:

    #BCs
    calc_inland_head = 1 #don't calculate from hgrad
    manual_inland_head = Lz + .1
    ocean_shead = [Lz,Lz-.5]
    ocean_ehead = ocean_shead
    start_fresh_yn = 1

    # save cell fluxes to unit 53
    ipakcb = 53

    #MODFLOW BCs
    hgrad = 0.0033 #hydraulic gradient, m/m
    qinflow = 0  #influent FW m3/day

    #Create ocean boundary at top of model
    ocean_loc = [30,70] # location of ocean [start,end] in m
    ocean_col = [int(find_nearest((np.arange(ncol)*delc),loc)) for loc in ocean_loc]
    ocean_col_vec = (0,0,np.arange(ocean_col[0],ocean_col[1]+1))
    ocean_coords = (0,slice(0,nrow),slice(ocean_col[0],ocean_col[1]+1)) #top of the model
    ocean_bool = np.zeros((nlay,nrow,ncol), dtype=bool)
    ocean_bool[0,:,np.arange(ocean_col[0],ocean_col[1])] = True
    m.ocean_bool = ocean_bool
    m.ocean_arr = ocean_bool

    if calc_inland_head ==1:
        head_inland = ocean_loc[0]*hgrad + Lz
    else:
        head_inland = manual_inland_head
    #Set ibound
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, 0] = -1 #first column (FW boundary) has fixed head
    ibound[ocean_bool] = -1
    ibound[0:2,:,ocean_col[0]-3:ocean_col[0]] = 0

    #Set starting heads
    strt = head_inland*np.ones((nlay, nrow, ncol)) #starting heads (for fixed head BC, this is will assign the head)
    strt[ocean_bool] = ocean_shead[0] #head of ocean boundary

    #Transport BCs
    if start_fresh_yn == 1:
        sconc = Cfresh*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin fresh
    elif start_fresh_yn == 0:
        sconc = Csalt*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin SW-saturated

    sconc[ocean_bool] = Csalt
    sconc[:,:,0] = Cfresh

    icbund = np.ones((nlay, nrow, ncol), dtype=np.int) 
    icbund[(ibound < 0)] = -1 #constant concentration cells where also constant head

    #Constant head boundary
    stress_period_data = {}
    ocean_sub = np.where(ocean_bool)
    for i in range(nper):
        dat = []
        for j in range(np.sum(ocean_bool)):
            dat.append([ocean_sub[0][j],ocean_sub[1][j],ocean_sub[2][j],ocean_shead[i],ocean_ehead[i]])
        stress_period_data[i] = dat


    #Create instances in flopy
    bas = flopy.modflow.ModflowBas(m, ibound, strt=strt)

    # Add LPF package to the MODFLOW model
    lpf = flopy.modflow.ModflowLpf(m, hk=hk, vka=vka, ipakcb=ipakcb, ss=ss, sy=sy)

    # Add PCG Package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(m, hclose=1.e-8)

    # Add OC package to the MODFLOW model
    oc = flopy.modflow.ModflowOc(m,
                                 stress_period_data={(0, 0): ['save head', 'save budget']},
                                 compact=True)
    #wel = flopy.modflow.ModflowWel(m, stress_period_data=wel_data, ipakcb=ipakcb)

    #Create the basic MT3DMS model structure
    btn = flopy.mt3d.Mt3dBtn(m, 
                             #nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, 
                             laycon=lpf.laytyp, htop=henry_top, 
                             dz=dis.thickness.get_value(), prsity=0.2, icbund=icbund,
                             sconc=sconc, nprs=-10)
    adv = flopy.mt3d.Mt3dAdv(m, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(m, al=al, trpt=1., trpv=1., dmcoef=dmcoef)
    gcg = flopy.mt3d.Mt3dGcg(m, iter1=500, mxiter=1, isolve=1, cclose=1e-7)
    ssm = flopy.mt3d.Mt3dSsm(m, stress_period_data=ssm_data)
    chd = flopy.modflow.ModflowChd(m, stress_period_data=stress_period_data)

    #vdf = flopy.seawat.SeawatVdf(m, iwtable=0, densemin=0, densemax=0,denseref=1000., denseslp=0.7143, firstdt=1e-3)
    vdf = flopy.seawat.SeawatVdf(m, mtdnconc=1, mfnadvfd=1, nswtcpl=0, iwtable=0, 
                                 densemin=0., densemax=0., denseslp=denseslp, denseref=densefresh)
    #Write input
    m.write_input()
    #Create a storage dictionary to pass to SGD model
    storage_dict = {'ocean_col': ocean_col,
    'model_ws':m.model_ws,
    'modelname':m.name,
    'ocean_bool': ocean_bool,
    'head_inland': head_inland,
    'ocean_shead': ocean_shead,
    'start_fresh_yn': start_fresh_yn
    }
    m.ocean_bool = ocean_bool
    m.set_ocean_arr(ocean_col_vec)
    m.set_storage_dict(storage_dict)
    return m, ocean_col

def get_model():
    fname = os.path.join(model_ws,modelname + '.nam')
    m = flopy.seawat.Seawat.load(fname,exe_name = sw_exe, model_ws = model_ws)
    return m
