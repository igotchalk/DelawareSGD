#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 16:20:56 2019

@author: ianpg
"""



from scipy.spatial.distance import pdist, squareform
from scipy import exp
from scipy.linalg import eigh
import numpy as np

def stepwise_kpca(X, gamma, n_components):
    """
    Implementation of a RBF kernel PCA.

    Arguments:
        X: A MxN dataset as NumPy array where the samples are stored as rows (M),
           and the attributes defined as columns (N).
        gamma: A free parameter (coefficient) for the RBF kernel.
        n_components: The number of components to be returned.

    Returns the k eigenvectors (alphas) that correspond to the k largest
        eigenvalues (lambdas).

    """
    # Calculating the squared Euclidean distances for every pair of points
    # in the MxN dimensional dataset.
    sq_dists = pdist(X, 'sqeuclidean')

    # Converting the pairwise distances into a symmetric MxM matrix.
    mat_sq_dists = squareform(sq_dists)

    # Computing the MxM kernel matrix.
    K = exp(-gamma * mat_sq_dists)

    # Centering the symmetric NxN kernel matrix.
    N = K.shape[0]
    one_n = np.ones((N,N)) / N
    K_norm = K - one_n.dot(K) - K.dot(one_n) + one_n.dot(K).dot(one_n)

    # Obtaining eigenvalues in descending order with corresponding
    # eigenvectors from the symmetric matrix.
    eigvals, eigvecs = eigh(K_norm)

    # Obtaining the i eigenvectors (alphas) that corresponds to the i highest eigenvalues (lambdas).
    alphas = np.column_stack((eigvecs[:,-i] for i in range(1,n_components+1)))
    lambdas = [eigvals[-i] for i in range(1,n_components+1)]

    return alphas, lambdas
#%%
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



nlay=m.dis.nlay
nrow=m.dis.nrow
ncol=m.dis.ncol

if heterogenous:
    import simulationFFT
    fft_grid = np.exp(simulationFFT.simulFFT(nrow, nlay, ncol, mu, sill, modeltype, lrow , llay, lcol))
    #hk[0:int(np.where(henry_botm==find_nearest(henry_botm,ocean_elev))[0])+1,:,:] = hkSand
else:
    hk = hkSand*np.ones((nlay,nrow,ncol), dtype=np.int32)

grid = np.log10(fft_grid)
#lith_props = [0.2,0.5,0.3]
#hk_vals = [-1,0,2]
lith_props = [0.5,0.5]
hk_vals = [0,2]

log10trans = True
plotyn= False
hk = truncate_grf(grid,lith_props,hk_vals,log10trans=True,plotyn=plotyn)
#%%

def make_hk():
    #hk1
    low= -2
    high = 0
    parname='hk1'
    val = sample_dist(sts.uniform,1,*(low,high-low))
    hk1 = 10**val

    #hk2
    low= 0
    high = 2
    parname='hk2'
    val = sample_dist(sts.uniform,1,*(low,high-low))
    hk2 = 10**val

    #lith_prop
    low= 0
    high = 1
    parname='lith_prop'
    val = sample_dist(sts.uniform,1,*(low,high-low))
    lith_prop = val

    #vario_type
    parname='vario_type'
    val = int(round(np.random.rand()))
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
    corr_len = val

    #corr_len_zx
    # equal to lz/lx
    low= .01
    high = .1
    parname='corr_len_zx'
    val = sample_dist(sts.uniform,1,*(low,high-low))
    corr_len_zx = val

    #corr_len_yx
    # equal to ly/lx
    low= 0.1
    high = 1
    parname='corr_len_yx'
    val = sample_dist(sts.uniform,1,*(low,high-low))
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
    hk_vals = [0,1]
    hk = truncate_grf(grid,lith_props,hk_vals,log10trans=False,plotyn=False)
    hk[wel_cells] = hk.max()
    return hk


def change_lith_prop_only(grf,ngrids,lith_prop_range):
    hk_mat = np.zeros(np.hstack((ngrids,np.shape(grf))))
    for i in range(ngrids):
        lith_prop = np.random.uniform(low=lith_prop_range[0],high=lith_prop_range[1])
        lith_props = [lith_prop,1-lith_prop]
        hk_mat[i] = truncate_grf(grf,lith_props,[0,1],log10trans=False,plotyn=False)
    return hk_mat

#%%



#%%

from sklearn.decomposition import KernelPCA

#Make hk_mat
ngrids=500
#hk_mat = change_lith_prop_only(grf,ngrids,[0,1])


hk_mat = np.zeros((ngrids,nlay,nrow,ncol),dtype=np.float)
for i in range(ngrids):
    hk_mat[i] = make_hk()


ndims=20
X = hk_mat.reshape((ngrids,-1),order='C')


print('Shape of X:',X.shape)

scikit_kpca = KernelPCA(n_components=2, kernel='rbf', gamma=5)
X_skernpca = scikit_kpca.fit_transform(X)




#alphas,lambdas = stepwise_kpca(X,1,ndims)

#%%
ndims=10

from sklearn.decomposition import PCA
scikit_pca = PCA(n_components=ndims)
X_spca = scikit_pca.fit_transform(X)
sortind = np.argsort(X_spca[:,0])


#%%
def scatter_minmax(alphas):
    minx,miny = np.argmin(alphas,axis=0)[:2]
    maxx,maxy = np.argmax(alphas,axis=0)[:2]
    minsmaxs = np.hstack((minx,miny,maxx,maxy))
    ttls = ['min x','min y','max x','max y']
    marker = ['x','o','v','s']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(alphas[:,0],alphas[:,1])
    for i in range(len(minsmaxs)):
        ax.scatter(alphas[minsmaxs[i],0],alphas[minsmaxs[i],1],color='red',marker=marker[i])
    plt.show()

    plt.figure()
    for i in range(4):
        plt.subplot(2, 2, i+1, aspect='equal')
        plt.imshow(hk_mat[minsmaxs[i],:,0,:],vmin=0,vmax=1)
        plt.title(ttls[i])
    return minsmaxs

def plot_across_axis(alphas,component=0,numplots=10):
    sortind = np.argsort(alphas[:,component])
    plotinds = np.linspace(0,len(alphas)-1,numplots,dtype=np.int)
    for i in range(numplots):
        print(i)
        f,ax= plt.subplots(1)
        plt.imshow(hk_mat[sortind[plotinds[i]],:,0,:],vmin=hk_mat.min(),vmax=hk_mat.max())
        plt.title('Component no. {} \nPlot no. {} of {}'.format(component,i+1,numplots))
        plt.show()
    return





minsmaxs = scatter_minmax(X_spca)
plot_across_axis(X_spca,component=0,numplots=2)



#%%
X = X_spca.copy()
size = ngrids
som = MiniSom(size, 1, len(X[0]),
              neighborhood_function='gaussian', sigma=.5,
              random_seed=1)

#som.pca_weights_init(X)
som.random_weights_init(X)
wghts_init = som.get_weights().copy().squeeze()
som.train_random(X, 1000, verbose=True)
#%%
rnk = np.zeros((ngrids),dtype=np.int)
grd = []
for cnt, xx in enumerate(X):
    #rnk[cnt] = som.winner(xx)  # getting the winner
    grd.append(som.winner(xx))

alpha_gridcell = np.r_[grd].T


#%% More plotting
ranks = np.arange(0,size,10,dtype=np.int)
inds=[]
pltranks=[]
for ind in ranks:
    val = np.where(alpha_gridcell[0]==ind)[0]
    if len(val>0):
            inds.append(int(val[0]))
            pltranks.append(ind)

plt.figure(figsize=(7, 7))
# Plotting the response for each pattern in the iris dataset
plt.scatter(X[:,0],X[:,2],c=alpha_gridcell[0])
plt.colorbar()
plt.show()

#%%
select_rank = np.where(rnks==rnks.max())[0]

for i in range(len(select_rank)):
    print(i)
    f,ax= plt.subplots(1)
    plt.imshow(hk_mat[select_rank[i],:,0,:],vmin=hk_mat.min(),vmax=hk_mat.max())
    plt.title('Component no. {} \nPlot no. {} of {}'.format(component,i+1,numplots))
    plt.show()

#%% If map is 2-D
numplots = 20

sortind = np.argsort(rnks)
plotinds = np.linspace(0,len(rnks)-1,numplots,dtype=np.int)
for i in range(numplots):
    print(i)
    f,ax= plt.subplots(1)
    plt.imshow(hk_mat[sortind[plotinds[i]],:,0,:],vmin=hk_mat.min(),vmax=hk_mat.max())
    plt.title('Component no. {} \nPlot no. {} of {}'.format(component,i+1,numplots))
    plt.show()




