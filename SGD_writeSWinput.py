import os
import sys
import numpy as np
from __main__ import *
import flopy

base_dir = os.path.abspath(os.getcwd())
model_ws = os.path.join(base_dir,'work',modelname)
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

            
    Lx = 70.
    Lz = 20.
    nlay = 40
    nrow = 1
    ncol = 70
    delr = Lx / ncol
    delc = 1.0
    delv = Lz / nlay
    henry_top = Lz
    henry_botm = np.linspace(henry_top - delv, 0., nlay)

    #Period data
    nper = 1
    perlen = [100]
    nstp = [100]
    steady = [True]
    itmuni = 4 #time unit 4= days
    lenuni = 2 #length unit 2= meter
    tsmult = 1
    ssm_data = None
    verbose=True


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
    calc_inland_head = 0 #don't calculate from hgrad
    manual_inland_head = Lz + .1
    ocean_head = Lz
    start_fresh_yn = 0

    def calc_head_inland(calc_head_yn,manual=1):
        if calc_head_yn==1:
            head_inland = ocean_loc[0]*hgrad #calculate inland head from measured head gradient
        else:
            head_inland = manual
        return head_inland


    # save cell fluxes to unit 53
    ipakcb = 53

    #MODFLOW BCs
    hgrad = 0.0033 #hydraulic gradient, m/m
    qinflow = 0  #influent FW m3/day
    head_inland = calc_head_inland(calc_inland_head,manual_inland_head)

    #Create ocean boundary at top of model
    ocean_loc = [30,70] # location of ocean [start,end] in m
    ocean_col = [int(find_nearest((np.arange(ncol)*delc),loc)) for loc in ocean_loc]
    ocean_col_vec = (0,0,np.arange(ocean_col[0],ocean_col[1]+1))
    ocean_bool = np.zeros((nlay,nrow,ncol), dtype=bool)
    ocean_bool[ocean_col_vec] = 1
    ocean_coords = (0,0,slice(ocean_col[0],ocean_col[1]+1)) #top of the model

    #Set ibound
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, 0] = -1 #first column (FW boundary) has fixed head
    ibound[ocean_col_vec] = -1
    ibound[0:2,0,ocean_col[0]-3:ocean_col[0]] = 0

    #Set starting heads
    strt = head_inland*np.ones((nlay, nrow, ncol)) #starting heads (for fixed head BC, this is will assign the head)
    strt[ocean_col_vec] = ocean_head #head of ocean boundary

    #Transport BCs
    if start_fresh_yn == 1:
        sconc = Cfresh*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin fresh
    elif start_fresh_yn == 0:
        sconc = Csalt*np.ones((nlay, nrow, ncol), dtype=np.float32) #Begin SW-saturated

    sconc[ocean_col_vec] = Csalt
    sconc[:,0,0] = Cfresh

    icbund = np.ones((nlay, nrow, ncol), dtype=np.int) 
    icbund[(ibound < 0)] = -1 #constant concentration cells where also constant head


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

    #vdf = flopy.seawat.SeawatVdf(m, iwtable=0, densemin=0, densemax=0,denseref=1000., denseslp=0.7143, firstdt=1e-3)
    vdf = flopy.seawat.SeawatVdf(m, mtdnconc=1, mfnadvfd=1, nswtcpl=0, iwtable=0, 
                                 densemin=0., densemax=0., denseslp=denseslp, denseref=densefresh)
    #Write input
    m.write_input()
    #Create a storage dictionary to pass to SGD model
    storage_dict = {'ocean_col': ocean_col,
    'ocean_bool': ocean_bool,
    'head_inland': head_inland,
    'ocean_head': ocean_head,
    'start_fresh_yn': start_fresh_yn
    }

    m.set_ocean_arr(ocean_col_vec)
    m.set_storage_dict(storage_dict)
    return m, ocean_col

def get_model():
    fname = os.path.join(model_ws,modelname + '.nam')
    m = flopy.seawat.Seawat.load(fname,exe_name = sw_exe, model_ws = model_ws)
    return m
def plot(m):
    import matplotlib.pyplot as plt
    import matplotlib.colors
    import numpy as np
    # Make plot of the grid
    f = plt.figure(figsize=(15, 5))
    plt.clf()
    ax = f.add_subplot(1, 1, 1)
    mm = flopy.plot.ModelCrossSection(ax=ax, model=m, line={'row':0});

    lpf = m.get_package('lpf')
    hk = lpf.hk.array
    hkpatchcollection = mm.plot_array(hk, norm=matplotlib.colors.LogNorm(),vmin=np.min(hk), vmax=np.max(hk));
    linecollectdsion = mm.plot_grid();
    patchcollection = mm.plot_ibound();
    cb = plt.colorbar(patchcollection);
    cb.set_label('Boundary condition',rotation=90)

    cb.set_ticks((1.5,2.5))
    #cb.set_ticklabels(('No flow','Const head'))
    cb.ax.set_yticklabels(('No flow','Const head'),rotation=90)


    cb2 = plt.colorbar(hkpatchcollection,ax=ax);
    cb2.set_label('Kh (m/d)', rotation=90)
    plt.title('K-field & Boundary conditions');
    plt.show()
    return