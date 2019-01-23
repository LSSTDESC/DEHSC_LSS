#Imports
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

### Input weights file
#This file also has COSMOS photo-zs, HSC photo-z's, positions and  hsc magnitudes). 
#Note: these weights were found based on the density of the whole HSC sample in the COSMOS field
#The COSMOS 30-band data matched to some of the HSC data in the COSMOS field

weights_path= '../../../cosmos_hsc_weights.fits'
hdulist = fits.open(weights_path)
weights_data = hdulist[1].data
weights_columns= hdulist[1].columns

weights= weights_data['weight']
hsc_z= weights_data['hsc_pz_best_eab']
COSMOS_z= weights_data['cosmos_photoz']

#Input the number of galaxies in the whole HSC field
N_photo_tot = 537547

#Input redshift bin info - weights will be summed in each redshift bin
num_bin= 100
bins= np.linspace(0, 4.0, num=100, endpoint=False)
bins=np.append(bins,4)

weights_tot= np.sum(weights)  #Sum weights over all matched COSMOS galaxies



### N(z) with whole galaxy sample

def weights_per_bin(num_bin,COSMOS_z, bins, weights):
    
    #Weights per redshift bin  
    weights_bin= np.zeros(num_bin)

    for j in range(num_bin):
        for i in range(len(COSMOS_z)):
            if COSMOS_z[i]>bins[j] and COSMOS_z[i]<=bins[j+1]:
                weights_bin[j]+=weights[i]
                
    return weights_bin
    
    
def Nz(weights_bin, weights_tot, N_photo_tot):
    P_T_z_delta_z =  np.true_divide(weights_bin, weights_tot)
    N_p_est = P_T_z_delta_z*N_photo_tot
    return N_p_est
    
    
    
    
weights_bin= weights_per_bin(num_bin, COSMOS_z, bins, weights)

N_p_est= Nz(weights_bin, weights_tot, N_photo_tot)


