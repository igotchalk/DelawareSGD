#%%

def load_obj(dirname,name):
    import pickle
    with open(dirname.joinpath(name + '.pkl').as_posix(), 'rb') as f:
        return pickle.load(f)

def mod_hausdorff(u,v):
    from scipy.spatial.distance import directed_hausdorff
    return max(directed_hausdorff(u,v)[0], directed_hausdorff(v,u)[0])

def compute_export_hausdorff(dirname,saveyn=1):
    import glob
    import flopy
    import numpy as np
    from pathlib import Path
    from scipy.io import savemat
    Csalt = 35.0001
    Cfresh = 0
    #dirname = Path('/Users/ianpg/Documents/ProjectsLocal/DelawareSGD/work/homogenous/MC_expt_2018-12-21-01-47')
    dirname=Path(dirname)
    mname = dirname.parts[-2]
    m = flopy.modflow.Modflow.load(dirname.parent.joinpath(mname + '.nam'))
    m.inputParams = load_obj(dirname,'inputParams')


    #make a point cloud for each conc array
    conc_fnames = glob.glob(dirname.joinpath('conc*.npy').as_posix())

    pct50 = (Csalt+Cfresh)/2
    tol = .20 #percentage
    conc_mat = np.zeros((len(conc_fnames),m.nlay,m.nrow,m.ncol),dtype=float)

    idx_dict = {}
    for i,fname in enumerate(conc_fnames):
        conc_mat[i] = np.load(fname)
        idx_dict[i] = np.where((conc_mat[i]<pct50*(1+tol)) & (conc_mat[i]>pct50*(1-tol)))

    #Filter out keys that are empty and make new dict with filled sequential keys
    filt = [k for k, v in idx_dict.items() if len(v[0])==0]
    i=0
    idx_dict_filt = {}
    for k,v in idx_dict.items():
        if k in filt:
            print('iteration ',k,'is empty and will be filtered out')
        else:
            idx_dict_filt[i] = v
            i+=1

    #Do the same with inputParams dictionary
    i=0
    inputParams_filt = {}
    for k,v in m.inputParams.items():
        vnew = [x for i, x in enumerate(v) if i not in filt]
        inputParams_filt[k] = vnew
        if i==0:
            N = len(vnew)
        i+=1

    #Calc modified hausdorff matrix
    #idx_dict = {k: v for k, v in idx_dict.items() if len(v[0])!=0} #Filter out empty keys. Might want to tie this to the pertinent files
    hdorf_mat = np.zeros((len(idx_dict_filt),len(idx_dict_filt)),dtype=float)
    hdorf_list= []
    for i in range(len(idx_dict_filt)):
        for j in range(i+1,len(idx_dict_filt)):
            hdorf_list.append(mod_hausdorff(idx_dict_filt[j],idx_dict_filt[i]))

    hdorf_mat = hdorf_mat + hdorf_mat.T
    #Save hausdorff matrix, Input Parameter Values
    ParametersValues = np.zeros((N,len(inputParams_filt)))
    for i,key in enumerate(inputParams_filt):
        ParametersValues[:,i] = np.asarray(inputParams_filt[key])

    if saveyn==1:
        hdorf_matdict = {}
        hdorf_matdict['hausdorff_mat'] = hdorf_list #list of distance calculations (not squareformed)
        hdorf_matdict['InputParams'] = inputParams_filt #Dict of the names and values, mostly for error checking
        hdorf_matdict['ParametersValues'] = ParametersValues
        np.save(dirname.joinpath('hausdorff.npy'),hdorf_list)
        savemat(dirname.joinpath('hausdorff.mat').as_posix(),hdorf_matdict)
        return
    else:
        return hdorf_matdict
