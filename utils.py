#Post-processing functions
import os
import numpy as np
import flopy
from pathlib import Path

def plotdischarge(modelname,model_ws,color='w',kstpkper_ind=-1,scale=50):
    fname = os.path.join(model_ws, '' + modelname + '.cbc')
    budobj = flopy.utils.CellBudgetFile(fname)
    qx = budobj.get_data(text='FLOW RIGHT FACE')[kstpkper_ind]
    qz = budobj.get_data(text='FLOW LOWER FACE')[kstpkper_ind]
    
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
              qx_avg[::iskip, 0, ::iskip], -qz_avg[::iskip, 0, ::iskip],
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

#calculate ocean flow from the cell-by-cell budget file
def get_ocean_outflow(m,ocean_coords):
    fname = os.path.join(m.model_ws, m.name + '.cbc') 
    budobj = flopy.utils.CellBudgetFile(fname)
    qx = budobj.get_data(text='FLOW RIGHT FACE')[-1]
    qz = budobj.get_data(text='FLOW LOWER FACE')[-1]
    ocean_flow = np.asarray(qz[ocean_coords])
    return ocean_flow

#calculate the ocean outflows from the constant head file
def get_ocean_outflow_chd(m,ocean_bool):
    fname = os.path.join(m.model_ws,m.name + '.cbc')
    budobj = flopy.utils.CellBudgetFile(fname)
    ch = budobj.get_data(text='CONSTANT HEAD')
    ocean_ind_MF = np.ravel_multi_index(np.where(ocean_bool==1),(m.nlay,m.nrow,m.ncol))+1 #ones-based
    #Get flux from .cbc file
    flx = []
    for node,val in ch[-1]:
        if node in ocean_ind_MF:
            flx.append(-val)
    #Assign to grid
    ocean_outflow = np.zeros((m.nlay,m.nrow,m.ncol))
    ocean_outflow[np.where(ocean_bool==1)] = flx
    return ocean_outflow

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
    print('label: ',label)
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
    
def read_ref(fname='ref_file.txt'):
    import re
    #Load ref file to get info about model
    reffile = os.path.join('.',fname)
    reftext = open(reffile, 'r').read()    
    beg = [m.start() for m in re.finditer('<<<', reftext)]
    betw = [m.start() for m in re.finditer('>>>', reftext)]
    end = [m.start() for m in re.finditer('\n', reftext)]
    d = {}
    for i in range(len(beg)):
        d[str(reftext[beg[i]+3:betw[i]])] =  reftext[betw[i]+3:end[i]] 
    return d
    #[exec(reftext[beg[i]+3:betw[i]]) for i in range(len(beg))]