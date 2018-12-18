#Post-processing functions
import os
import numpy as np
import flopy
from pathlib import Path

def plotdischarge(m,color='w',per=-1,scale=50,rowslice=0,iskip=1):
    import matplotlib.pyplot as plt
    fname = os.path.join(m.model_ws, '' + m.name + '.cbc')
    budobj = flopy.utils.CellBudgetFile(fname)
    qx = budobj.get_data(text='FLOW RIGHT FACE')[per]
    qz = budobj.get_data(text='FLOW LOWER FACE')[per]

    # Average flows to cell centers
    qx_avg = np.empty(qx.shape, dtype=qx.dtype)
    qx_avg[:, :, 1:] = 0.5 * (qx[:, :, 0:m.ncol-1] + qx[:, :, 1:m.ncol])
    qx_avg[:, :, 0] = 0.5 * qx[:, :, 0]
    qz_avg = np.empty(qz.shape, dtype=qz.dtype)
    qz_avg[1:, :, :] = 0.5 * (qz[0:m.nlay-1, :, :] + qz[1:m.nlay, :, :])
    qz_avg[0, :, :] = 0.5 * qz[0, :, :]

    y, x, z = m.dis.get_node_coordinates()
    X, Z = np.meshgrid(x, z[:, 0, 0])

    ax = plt.gca()
    cpatchcollection = ax.quiver(X[::iskip, ::iskip], Z[::iskip, ::iskip],
              qx_avg[::iskip, rowslice, ::iskip], -qz_avg[::iskip, rowslice, ::iskip],
              color=color, scale=scale, headwidth=8, headlength=3,
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
def get_ocean_outflow_chd(m,ocean_bool=None,tot_stp=None):
    if ocean_bool is None:
        ocean_bool = m.ocean_arr
    fname = os.path.join(m.model_ws,m.name + '.cbc')
    budobj = flopy.utils.CellBudgetFile(fname)
    ch = budobj.get_data(text='CONSTANT HEAD')

    #ocean_ind_MF = np.ravel_multi_index(np.where(ocean_bool==1),(m.nlay,m.nrow,m.ncol))+1 #ones-based
    ocean_ind_MF = np.ravel_multi_index(ocean_bool,(m.nlay,m.nrow,m.ncol))+1 #ones-based
    if tot_stp is None:
        tot_stp = len(ch)-1
    #Get flux from .cbc file
    flx = []
    for node,val in ch[tot_stp]:
        if node in ocean_ind_MF:
            flx.append(-val)
    #Assign to grid
    ocean_outflow = np.zeros((m.nlay,m.nrow,m.ncol))
    #ocean_outflow[np.where(ocean_bool==1)] = flx
    ocean_outflow[ocean_bool] = flx
    return ocean_outflow

def get_salt_outflow(m,kstpkper=None):
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