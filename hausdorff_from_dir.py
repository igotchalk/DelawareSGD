#%%

def load_obj(dirname,name):
    import pickle
    with open(dirname.joinpath(name + '.pkl').as_posix(), 'rb') as f:
        return pickle.load(f)

def save_obj(dirname,obj, name ):
    import pickle
    with open(dirname.joinpath(name + '.pkl').as_posix(), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def mod_hausdorff(u,v):
    from scipy.spatial.distance import directed_hausdorff
    return max(directed_hausdorff(u,v)[0], directed_hausdorff(v,u)[0])


def idx2centroid(node_coord_tuple,idx_tuple):
    z_pt = node_coord_tuple[2][idx_tuple]
    x_pt = node_coord_tuple[1][idx_tuple[2]]
    y_pt = node_coord_tuple[0][idx_tuple[1]]
    return (z_pt,y_pt,x_pt)


def show_concplots(conc_mat,startind=0,rowslice=0,numplots=10,saveyn=0,dirname=None,contoursyn=True):
    import numpy.ma as ma
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    Csalt = 35.0001
    Cfresh = 0
    #pctage = np.array([.05,.5,.95])
    pctage = np.array([.5])
    lvls = Cfresh + (Csalt-Cfresh)*np.array([.05,.5,.95])
    for i in range(startind+numplots-1,startind-1,-1):
        if i>conc_mat.shape[0]:
            break
        for k in range(len(pctage)):
            f, axs = plt.subplots()
            mskd = ma.masked_where((conc_mat[i]>1e10),conc_mat[i])
            arr_plot = mskd[:,rowslice,:]
            cpatch = plt.imshow(arr_plot,cmap='jet')
            plt.title('iter. {}, slice {} \nMasked points shown in white'.format(i,rowslice))
            plt.xlabel('Col. number')
            plt.ylabel('Lay. number')
            if contoursyn:
                CS = plt.contour(arr_plot,levels=lvls,colors='white')
                plt.clabel(CS, CS.levels, inline=True, fontsize=10)
            cbar_ax = f.add_axes([0.95, 0.3, 0.02, 0.4])
            cb = f.colorbar(cpatch,cax=cbar_ax)
            #cb = plt.colorbar(cpatch)
            cb.set_label('Salt concentration (g/L)')
            plt.show()
            if saveyn==1:
                plt.savefig(Path(dirname).joinpath('conc_it{}_slice{}.png'.format(i,rowslice)),dpi=150)
    return


def compute_export_hausdorff(dirname,conc_mat=None,yxz=None,saveyn=1,pct=(.05,.95)):
    import glob
    import os
    import flopy
    import numpy as np
    from pathlib import Path
    from scipy.io import savemat

    hdorf_matdict = {}
    Csalt = 35.0001
    Cfresh = 0
    pct50 = (Csalt+Cfresh)/2
    tol = .20 #percentage
    dirname=Path(dirname)
    mname = dirname.parts[-2]

    try:
        m = flopy.modflow.Modflow.load(str(dirname.parent.joinpath(mname + '.nam')))
    except:
        print('cant load modflow model... assuming model is size 26x20x100')
        m = flopy.modflow.Modflow('temp')
        flopy.modflow.ModflowDis(m, 26,20,100)

    #Load a conc_mat from .npy files in directory if not specified
    if conc_mat is None:
        #make a point cloud for each conc array
        conc_fnames = sorted(glob.glob(dirname.joinpath('conc*.npy').as_posix()),
                             key=os.path.getctime)

        conc_mat = np.zeros((len(conc_fnames),m.nlay,m.nrow,m.ncol),dtype=float)
        for i,fname in enumerate(conc_fnames):
            conc_mat[i] = np.load(fname)

    from skimage import measure
    idx_dict = {}
    face_dict = {}
    if pct is None:
        pct = (.05,.95)
    for i in range(len(conc_mat)):
        print('it',i)
        idx_dict[i] = np.where((conc_mat[i]<pct50*(1+tol)) & (conc_mat[i]>pct50*(1-tol)))
        faces = []
        for p in pct:
            verts, f, normals, values = measure.marching_cubes_lewiner(conc_mat[i],(Csalt+Cfresh)*p)
            faces.append(f[np.arange(1,len(f),3),:])
        face_dict[i] = list(map(tuple, np.vstack(faces)))


    #Filter out keys that have errors and make new dict with filled sequential keys
    filt_mat = np.where(conc_mat>1e10,np.nan,conc_mat) #get rid of region above water table (value = 1e30)
    filt = [k for k, v in idx_dict.items() if len(v[0])==0 #which indicies to get rid of
        or (filt_mat[k].max() > Csalt*1.05)
        or (filt_mat[k].min() < -.1)]
    i=0
    idx_dict_filt = {}
    for k,v in idx_dict.items():
        if k in filt:
            print('iteration ',k,'doesnt meet reqs and will be filtered out')
        else:
            idx_dict_filt[i] = v
            i+=1

    i=0
    face_dict_filt = {}
    for k,v in face_dict.items():
        if k in filt:
            print('iteration ',k,'doesnt meet reqs and will be filtered out')
        else:
            face_dict_filt[i] = v
            i+=1


    try:
        #Do the same with inputParams dictionary
        m.inputParams = load_obj(dirname,'inputParams')
        i=0
        inputParams_filt = {}
        for k,v in m.inputParams.items():
            vnew = [x for i, x in enumerate(v) if i not in filt]
            inputParams_filt[k] = vnew
            if i==0:
                N = len(vnew)
            i+=1
        ParametersValues = np.zeros((N,len(inputParams_filt)))
        for i,key in enumerate(inputParams_filt):
            ParametersValues[:,i] = np.asarray(inputParams_filt[key])
        hdorf_matdict['InputParams'] = inputParams_filt #Dict of the names and values, mostly for error checking
        hdorf_matdict['ParametersValues'] = ParametersValues
    except:
        print('Not able to import inputParams.pkl, cannot export inputParams file')
        pass

    #create pointset dictionary: (x1,y1,z1),...
    #   from the numpy tuple format: (x1,x2,...),(y1,y2,...),(z1,z2,...)
    if yxz is None:
        yxz = m.dis.get_node_coordinates()
    ptset_dict = {}
    for i in range(len(idx_dict_filt)):
        print('iteration',i)
        ptset_tuple = idx2centroid(yxz,idx_dict_filt[i])
        ptset_dict[i] = [(ptset_tuple[0][i],ptset_tuple[1][i],ptset_tuple[2][i]) for i in range(len(ptset_tuple[0]))]

    hdorf_mat = np.zeros((len(idx_dict_filt),len(idx_dict_filt)),dtype=float)
    hdorf_list= []
    for i in range(len(idx_dict_filt)):
        for j in range(i+1,len(idx_dict_filt)):
            #hdorf_list.append(mod_hausdorff(ptset_dict[j],ptset_dict[i]))
            hdorf_list.append(mod_hausdorff(face_dict_filt[j],face_dict_filt[i]))
            hdorf_mat[i,j] = hdorf_list[-1]
            print('Calculating hausdorff row ',i,' col ',j)
    hdorf_mat = hdorf_mat + hdorf_mat.T
    hdorf_matdict['hausdorff_mat'] = hdorf_list #list of distance calculations (not squareformed)
    hdorf_matdict['culled_conc_mat'] = np.delete(conc_mat, filt, 0)
    #Save hausdorff matrix, Input Parameter Values
    if saveyn==1:
        np.save(dirname.joinpath('hausdorff.npy'),hdorf_list)
        savemat(dirname.joinpath('hausdorff.mat').as_posix(),hdorf_matdict)
        return
    else:
        return hdorf_matdict,conc_mat

def create_concmat_from_ucndir(dirname,pattern='*.UCN',totims=(2340.0,4860.0,7200.0),modsize=(26,20,100),saveyn=1):
    import glob
    import os
    import flopy
    import numpy as np
    from pathlib import Path

    dirname=Path(dirname)
    ucn_fnames = sorted(glob.glob(dirname.joinpath(pattern).as_posix()),
                         key=os.path.getctime)
    if modsize is None:
        #check out concentration mat size
        ucnobj = flopy.utils.binaryfile.UcnFile(ucn_fnames[0])
        modsize =  np.shape(ucnobj.get_data(totim=ucnobj.get_times()[-1]))

    conc_mat = np.zeros((len(totims),len(ucn_fnames),modsize[0],modsize[1],modsize[2]),dtype=float)
    filt = []
    for i,fname in enumerate(ucn_fnames):
        print('file {} of {}'.format(i,len(ucn_fnames)))
        ucnobj = flopy.utils.binaryfile.UcnFile(fname)
        if ucnobj.get_times()[-1] != 360*20:
            filt.append(i)
            print('Simulation does not appear to have completed in file:\n    ' + Path(fname).parts[-1] +
                      '\n...removing file from conc_mat')
            continue

        for j,tim in enumerate(totims):
            try:
                conc_mat[j,i,:,:,:] = ucnobj.get_data(totim=tim)
            except:
                filt.append(i)
                print('requested totim not found in file:\n    ' + Path(fname).parts[-1] +
                      '\n...removing file from conc_mat')

    filt = list(np.unique(filt))
    print('Removing {} files out of {}...'.format(len(filt),len(ucn_fnames)))
    conc_mat_filt = np.zeros((len(totims),len(ucn_fnames)-len(filt),modsize[0],modsize[1],modsize[2]),dtype=float)
    it = 0
    for i in range(len(ucn_fnames)):
        if i not in filt:
            print('it',it)
            print('i',i)
            conc_mat_filt[:,it,:,:,:] = conc_mat[:,i,:,:,:]
            it+=1
    ucn_names_filt = [fname for i,fname in enumerate(ucn_fnames) if i not in filt]

    if saveyn==1:
        print('saving...')
        fnames = []
        for i,tim in enumerate(totims):
            fnames.append('conc_mat_totim' + str(int(tim)) + '.npy')
            np.save(dirname.joinpath(fnames[i]),conc_mat_filt[i,:,:,:,:])
        save_obj(dirname,ucn_names_filt,'ucn_fnames')
        save_obj(dirname,filt,'filt')
        print('...done!')
        return fnames
    else:
        return conc_mat




















