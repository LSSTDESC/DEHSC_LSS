#Imports
%matplotlib inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import random
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, hstack
from scipy.optimize import curve_fit
from sympy import *
import math
from mpmath import *
from pandas.compat import StringIO
from numpy.linalg import inv
from scipy import linalg
from scipy import interpolate
from scipy import integrate
from scipy.integrate import quad
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KDTree
import scipy.spatial as spatial



def make_samples(cosmos_matched_data, hsc_cosmos_data):
    cosmos_matched_df= pd.DataFrame(cosmos_matched_data)
    hsc_cosmos_df= pd.DataFrame(hsc_cosmos_data)

    #Training set, or HSC photometry that matches to COSMOS 30-band
    Train_sample= np.zeros([len(cosmos_matched_df),5])
    for i in range(len(cosmos_matched_df)):
        Train_sample[i]= np.array([cosmos_matched_df.iloc[i,536:541]]) #These are the magnitude values
        
    #COSMOS 30-band photometric redshifts
    Train_z = np.array(cosmos_matched_df['PHOTOZ'])
    
    #All HSC photometry in COSMOS field
    Photoz_sample= np.zeros([len(hsc_cosmos_df),5])
    for i in range(len(hsc_cosmos_df)):
        Photoz_sample[i]= np.array([hsc_cosmos_df.iloc[i,88:93]]) #These are the magnitude values
        
    return Train_sample, Train_z, Photoz_sample


def find_nn(N_nei, Train_sample):
    #Find the nearest neigbours for each of the training galaxies 
    N_nbrs = NearestNeighbors(n_neighbors=N_nei, algorithm='kd_tree', metric='euclidean').fit(Train_sample)

    #Get the distances to each of the neighbours
    distances, indices = N_nbrs.kneighbors(Train_sample)
    return distances, indices

def furthest_dist(Train_sample, distances):
    d_dist= np.zeros(len(Train_sample))
    #Make an array with the distances to the furthest neighbour for each training galaxy
    for i in range(len(Train_sample)):
        d_dist[i] = distances[i][len(distances[0])-1]
    return d_dist


def photoz_within_dist(Photoz_sample, Train_sample, d_dist):
    #Make the spatial.cKDTree object
    tree_NN_lookup = spatial.cKDTree(Photoz_sample, leafsize=40)

    #Make a corresponding array of the number of photo-z galaxies in the same area as N_nei training gals, for each train gal 
    num_photoz= np.zeros(len(Train_sample))

    #Query the indicies of the points that lie within a distance of the point given 
    for i in range(len(Train_sample)):
        num_photoz[i]= len(tree_NN_lookup.query_ball_point(Train_sample[i], d_dist[i]+0.00000100))
        
    return num_photoz


def make_weights(Train_sample, num_photoz, N_nei, N_photo_tot):
    
    #Make a corresponding array of the weights for each of the training galaxies
    weights= np.zeros(len(Train_sample))

    #Make the weights
    weights = np.true_divide(num_photoz,N_nei*N_photo_tot)
    return weights


def weights_per_bin(num_bin, Train_sample, Train_z, bins):
    
    #Weights per redshift bin  
    weights_bin= np.zeros(num_bin)

    for j in range(num_bin):
        for i in range(len(Train_sample)):
            if Train_z[i]>bins[j] and Train_z[i]<=bins[j+1]:
                weights_bin[j]+=weights[i]
                
    return weights_bin


def Nz(weights_bin, weights_tot, N_photo_tot):
    P_T_z_delta_z =  np.true_divide(weights_bin, weights_tot)
    N_p_est = P_T_z_delta_z*N_photo_tot
    return N_p_est



def bootstrap_resample(Train_sample):
    n= len(Train_sample)
    resample = np.random.choice(n,n)
    X_resample = Train_sample[resample]
    return X_resample

def Nz_bootstrap(n_bootstrap):
    N_p_est_boot= np.zeros([n_bootstrap, num_bin])
    for i in range(n_bootstrap):
        Train_resample= bootstrap_resample(Train_sample)
        distances, indices = find_nn(N_nei, Train_resample)
        d_dist = furthest_dist(Train_resample, distances)
        num_photoz= photoz_within_dist(Photoz_sample, Train_resample, d_dist)
        weights = make_weights(Train_resample, num_photoz, N_nei, N_photo_tot)
        weights_tot= np.sum(weights)  #Sum weights over all training galaxies
        weights_bin= weights_per_bin(num_bin, Train_resample, Train_z, bins)
        N_p_est_boot[i]= Nz(weights_bin, weights_tot, N_photo_tot)
    return N_p_est_boot


def make_uncertainty_limits(num_bin, N_p_est_boot):
    bin_lower= np.zeros(num_bin)
    bin_upper= np.zeros(num_bin)
    for i in range(num_bin):
        bin_lower[i], bin_upper[i]= np.percentile(N_p_est_boot[:,i], [16,84])
    return bin_lower, bin_upper


#The COSMOS 30-band data matched to some of the HSC data in the COSMOS field
cosmos_matched_path= '../../../cosmos_matched_hsc.fits'
hdulist = fits.open(cosmos_matched_path)
cosmos_matched_data = hdulist[1].data
cosmos_matched_columns= hdulist[1].columns

#All the HSC data in COSMOS field
hsc_cosmos_path= '/global/cscratch1/sd/damonge/HSC/HSC_processed/DEEP_COSMOS/DEEP_COSMOS_Catalog_i24.50.fits'
hdulist = fits.open(hsc_cosmos_path)
hsc_cosmos_data = hdulist[1].data
hsc_cosmos_columns= hdulist[1].columns

#Input the total number of training galaxies
N_train_tot= len(Train_sample)

#Input the total number of photometric galaxies
N_photo_tot= len(Photoz_sample)

#Define the fixed N_nei that I want to use 
N_nei = 20

#Input redshift bin info - weights will be summed in each redshift bin
num_bin= 150
bins= np.linspace(0, 6.0, num=150, endpoint=False)
bins=np.append(bins,6)

Train_sample, Train_z, Photoz_sample = make_samples(cosmos_matched_data, hsc_cosmos_data)
n_bootstrap= 1000
N_p_est_boot= Nz_bootstrap(n_bootstrap)



bin_lower, bin_upper= make_uncertainty_limits(num_bin, N_p_est_boot)
errors= np.vstack([abs(N_p_est-bin_lower), abs(bin_upper-N_p_est)])


plt.errorbar(bins[:-1], N_p_est, yerr=errors, xerr=None ,ecolor= 'k', elinewidth= 1)
plt.xlabel('Redshift, z')
plt.ylabel('N(z)')
plt.savefig('N(z)_errorbars.fig')

plt.plot(bins[:-1], N_p_est, 'k', linewidth= 0.5)
plt.fill_between(bins[:-1], bin_lower, bin_upper, color='b')
plt.xlabel('Redshift, z')
plt.ylabel('N(z)')
plt.savefig('N(z)_errorfilled.fig')