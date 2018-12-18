
# coding: utf-8

# # Delaware Case Study
# ## Simple Model of SGD
#
#

# In[1]:

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import warnings

#Name model
modelname = 'test'

# run installed version of flopy or add local path
try:
    import flopy
except:
    fpth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(fpth)
    import flopy
import SGD

print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('flopy version: {}'.format(flopy.__version__))

if sys.platform == "darwin":
    model_ws = os.path.join('/Users/ianpg/Documents/ProjectsLocal/DelawareSGD','work',modelname)
elif sys.platform == "win32":
    model_ws = os.path.join('E:\Projects\DelawareSGD','work',modelname)

if not os.path.exists(model_ws):
    os.makedirs(model_ws)
sys.path.append(os.path.join(model_ws,'..','..'))
import config
sw_exe = config.swexe
print('Model workspace:',os.path.abspath(model_ws))



#Other functions and imports

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



#Model grid
#Grid: 0.5 *1*1m • Size:70m*20m

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
perlen = [1,100]
nstp = [100,100]
nper = len(perlen)
steady = [True,False]
itmuni = 4 #time unit 4= days
lenuni = 2 #length unit 2 = meter
tsmult = 1
ssm_data = None
verbose = True



#Create basic model instance and dis pacakge
m = flopy.seawat.Seawat(modelname, exe_name=sw_exe, model_ws=model_ws,verbose=verbose)
SGD.ModelSGD.Seawat2SGD(m)
print(m.namefile)

# Add DIS package to the MODFLOW model
dis = flopy.modflow.ModflowDis(m, nlay, nrow, ncol, nper=nper, delr=delr,
                               delc=delc, laycbd=0, top=henry_top,
                               botm=henry_botm, perlen=perlen, nstp=nstp,
                               steady=steady,itmuni=itmuni,lenuni=lenuni,
                               tsmult=tsmult)


#Hydraulic conductivity field

hkSand = 80.  #horizontal hydraulic conductivity m/day
hkClay = 1.
lithmat = hkSand*np.ones((nlay,nrow,ncol), dtype=np.int32) #sandy background
addclay_yn = 1


def rand_clay_blocks(lithmat,hkClay,numblocks,sizeblocks):
    nlay,nrow,ncol = lithmat.shape
    lay_block = np.random.randint(1,nlay-sizeblocks[0],numblocks)
    row_block = np.random.randint(0,nrow-sizeblocks[1]+1,numblocks)
    col_block = np.random.randint(1,ncol-sizeblocks[2],numblocks)
    lithmat_blocks = lithmat.copy()
    for i in range(numblocks):
        block_coords = [slice(lay_block[i],lay_block[i]+sizeblocks[0]),
                        slice(row_block[i],row_block[i]+sizeblocks[1]),
                        slice(col_block[i],col_block[i]+sizeblocks[2])]
        lithmat_blocks[block_coords] = hkClay
    return lithmat_blocks

#add low conductivity regions
if addclay_yn == 1:
    lithmat = rand_clay_blocks(lithmat,hkClay,100,(2,1,5))
low_k_loc = (20,30)
low_k_col = loc_to_col(low_k_loc)
#lithmat[1:3,0,0:65] = hkClay





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
ocean_bool[0,:,np.arange(ocean_col[0],ocean_col[1]+1)] = 1
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

#saving concentrations at specified times
timprs = [k for k in range(1,np.sum(perlen),5)]


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
                         laycon=lpf.laytyp, htop=henry_top,
                         dz=dis.thickness.get_value(), prsity=0.2, icbund=icbund,
                         sconc=sconc, nprs=1,timprs=timprs)
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

# Try to delete the output files, to prevent accidental use of older files
try:
    os.remove(os.path.join(model_ws,'MT3D.CNF'))
    os.remove(os.path.join(model_ws,'MT3D001.MAS'))
    os.remove(os.path.join(model_ws, 'MT3D001.UCN'))
    os.remove(os.path.join(model_ws, modelname + '.hds'))
    os.remove(os.path.join(model_ws, modelname + '.cbc'))
except:
    pass



#Run model
v = m.run_model(silent=False, report=True)
for idx in range(-3, 0):
    print(v[1][idx])


# ## Post-processing results

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
        label = [ k for k,v in globals().items() if v is array][-1]
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


#

# Extract heads
fname = os.path.join(model_ws, '' + modelname + '.hds')
hdobj = flopy.utils.binaryfile.HeadFile(fname)
times = hdobj.get_times()
hdobj.get_kstpkper()
hds = hdobj.get_data(totim=times[-1])
hds[np.where(ibound != 1)] = np.nan

# Extract salinity
fname = os.path.join(model_ws, 'MT3D001.UCN')
ucnobj = flopy.utils.binaryfile.UcnFile(fname)
times = ucnobj.get_times()
conc = ucnobj.get_data(totim=times[-1])
conc[np.where(ibound != 1)] = np.nan







# In[14]:
# Make head and quiver plot
import utils
printyn = 0
rowslice = 0
f, axs = plt.subplots(2,1, sharex=True, figsize=(15, 8),
                     gridspec_kw = {'height_ratios':[1, 7]})
plt.tight_layout()

#Plot discharge and ibound
mm = flopy.plot.ModelCrossSection(ax=axs[1], model=m, line={'row':rowslice})
cpatchcollection,label = plot_background(mm,hds)
mm.plot_ibound()
plotdischarge(m.name,m.model_ws,color='white',per=-1,scale=1,rowslice=rowslice);
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')

#plot SGD
ocean_flow = 100 * np.asarray(utils.get_ocean_outflow_chd(m,ocean_bool))

plt.sca(axs[0]) #set current axes
plt.bar(np.arange(ocean_col[0],ocean_col[1]+1)+delc/2,
        ocean_flow[:,rowslice,:][ocean_bool[:,rowslice,:]])
axs[0].axhline()
axs[0].set_ylim(np.min(ocean_flow)*0.9,np.max(ocean_flow)*1.1)
plt.ylabel('SGD (cm/d)')
axs[0].annotate('total flow = ' + str(np.sum(ocean_flow)) + ' cm/d', xy=(.1,.5),
            xytext=(0.1, 0.5), xycoords='axes fraction', textcoords='axes fraction')

#Plot salt content of cells beneath outflow
ocean_below_ind = [np.where(ocean_bool)[0]+1, #layer below ocean
                   np.where(ocean_bool)[1],
                   np.where(ocean_bool)[2]]
ind = np.where(ocean_below_ind[1] == rowslice)
ocean_below_ind_slice = (ocean_below_ind[0][ind],ocean_below_ind[1][ind],ocean_below_ind[2][ind])
ocean_salt_outflow = get_salt_outflow(m,kstpkper = (nstp[-1]-1,nper-1))[ocean_below_ind_slice]


axs02 = axs[0].twinx()
plt.plot(ocean_col_vec[2]+delc/2 , ocean_salt_outflow,'r.')
plt.ylabel('conc of underlying (g/L)')

#align plots and set colorbar
f.subplots_adjust(right=0.93)
cbar_ax = f.add_axes([0.95, 0.1, 0.02, 0.7])
cb = f.colorbar(cpatchcollection,cax=cbar_ax)
cb.set_label(label)
if printyn == 1:
    plt.savefig(os.path.join(m.model_ws, m.name + '_flow_vectors.png'),dpi=300)
plt.show()


# In[16]:
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

writeyn = 1
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
rowslice = 0
times_ucn = ucnobj.get_times()
times_hds = hdobj.get_times()
kstpkper_ucn = ucnobj.get_kstpkper()
kstpkper_hds = hdobj.get_kstpkper()


#Plot discharge and ibound
f, axs = plt.subplots(3,1, sharex=True, figsize=(10, 6),
                     gridspec_kw = {'height_ratios':[1, 4,4]})
plt.tight_layout()

plt.sca(axs[1]) #set current axes
pos = axs[1].get_position()

mm = flopy.plot.ModelCrossSection(ax=axs[1], model=m, line={'row':rowslice})

#plot background
conc = sconc
conc[np.where(ibound != 1)] = np.nan

background_patch,label = plot_background(mm,conc)
mm.plot_ibound()

scale = 1
discharge_arrows = utils.plotdischarge(m,color='white',per=0,scale=scale,rowslice=rowslice);
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')

#second background
plt.sca(axs[2]) #set current axes
mm2 = flopy.plot.ModelCrossSection(ax=axs[2], model=m, line={'row':rowslice})
background_patch2,label = plot_background(mm2,hk,'K')
mm2.plot_ibound()
discharge_arrows = utils.plotdischarge(m,color='white',per=0,scale=scale,rowslice=rowslice);


#plot SGD
ocean_flow = 100 * np.asarray(utils.get_ocean_outflow_chd(m,ocean_bool))

plt.sca(axs[0]) #set current axes
flux_bar = plt.bar(np.arange(ocean_col[0],ocean_col[1]+1)+delc/2,
        ocean_flow[:,rowslice,:][ocean_bool[:,rowslice,:]])
axs[0].axhline()
axs[0].set_ylim(np.min(ocean_flow)*0.9,np.max(ocean_flow)*1.1)
plt.ylabel('SGD (cm/d)')
axs[0].annotate('total flow = ' + str(np.sum(ocean_flow)) + ' cm/d', xy=(.1,.5),
            xytext=(0.1, 0.5), xycoords='axes fraction', textcoords='axes fraction')

#Plot salt content of cells beneath outflow
ocean_below_ind = [np.where(ocean_bool)[0]+1, #layer below ocean
                   np.where(ocean_bool)[1],
                   np.where(ocean_bool)[2]]
ind = np.where(ocean_below_ind[1] == rowslice)
ocean_below_ind_slice = (ocean_below_ind[0][ind],ocean_below_ind[1][ind],ocean_below_ind[2][ind])
ocean_salt_outflow = get_salt_outflow(m,kstpkper = (nstp[-1]-1,nper-1))[ocean_below_ind_slice]


axs02 = axs[0].twinx()
salt_line, = plt.plot(ocean_col_vec[2]+delc/2 , ocean_salt_outflow,'r.')
plt.ylabel('conc \nunderlying (g/L)')

#align plots and set colorbar
f.subplots_adjust(right=0.88)
cbar_ax = f.add_axes([0.90, 0.1, 0.02, 0.7])
cb = f.colorbar(cpatchcollection,cax=cbar_ax)
cb.set_label(label)
plt.show()

#animation function
def animate(i):

    #background patch
    plt.sca(axs[1]) #set current axes
    conc = ucnobj.get_data(kstpkper=i)
    conc[np.where(ibound != 1)] = np.nan
    background_patch,label = utils.plot_background(mm,conc,'Concentration (g/L)')

    #working
    #discharge arrows
    tot_stp = kstpkper_hds.index(i)
    discharge_arrows = utils.plotdischarge(m,color='white',per=tot_stp,scale=scale,rowslice=rowslice);

    #flux bar
    plt.sca(axs[0]) #set current axes
    ocean_flow = 100 * np.asarray(utils.get_ocean_outflow_chd(m,ocean_bool,tot_stp=tot_stp))
    ocean_flow_slice = ocean_flow[:,rowslice,:][ocean_bool[:,rowslice,:]]
    for rect, ind in zip(flux_bar, range(len(flux_bar))):
        rect.set_height(ocean_flow_slice[ind])

    #salt line
    ocean_salt_outflow = get_salt_outflow(m,kstpkper = i)[ocean_below_ind_slice]
    salt_line.set_ydata(ocean_salt_outflow)
    plt.title('kstkper = ' + str(i))
    print('kstpkper =  ' + str(i))

    #second background
    if len(axs==3):
        plt.sca(axs[2]) #set current axes
        discharge_arrows = utils.plotdischarge(m,color='white',per=tot_stp,scale=scale,rowslice=rowslice);
    return (background_patch,salt_line,flux_bar,discharge_arrows,)

if writeyn==1:
    # Set up formatting for the movie files
    writer = animation.FFMpegWriter(fps=2, bitrate=-1)
    anim = animation.FuncAnimation(f,animate,frames=ucnobj.get_kstpkper(),
                                   blit=False)
                                   #interval=50,blit=False,fargs=fargs)
    anim.save("test.mp4", writer=writer)

#%%
## Reading in real data


import pandas
from pathlib import Path
ws = Path(m.model_ws)
basedir = ws.joinpath('..','..')
datadir = basedir.joinpath('data')
fname = datadir.joinpath('seepage_flow_2015.xlsx')
flowsheet = pandas.read_excel(fname,skiprows=1)
print('Columns: ',flowsheet.columns)


# In[ ]:


# #Issue: which "flow" is the correct flow?
#     #Can get flow from the cell-by-cell file in either FLOW RIGHT FACE and FLOW LOWER FACE
#     # Or from the constant head flow
#
# ocean_flow = get_ocean_outflow(m,ocean_col);
# print('total flow into ocean cells from const. head flow:',-np.sum(ocean_flow),'m^3/d')
#
# fname = os.path.join(model_ws, '' + modelname + '.cbc')
# budobj = flopy.utils.CellBudgetFile(fname)
# qx = budobj.get_data(text='FLOW RIGHT FACE')[-1]
# qz = budobj.get_data(text='FLOW LOWER FACE')[-1]
# tot_flow = np.sum( np.sqrt(np.square(-qz[ocean_coords]) + np.square(qx[ocean_coords])))
# print('Total flow from lower-right face in cbc file', tot_flow ,'m^3/d' )
