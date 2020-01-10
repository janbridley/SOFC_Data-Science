# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 14:58:55 2020

@author: jjb247
"""

# %% Load Packages ============================================================
# Package Import
#from threshold import threshold           # File I created to threshold images
from scipy.io import loadmat              # loads .mat file
#import matplotlib.pyplot as plt           # plots data: includes plot & imshow
#from mpl_toolkits import mplot3d          # supports 3d plotting
#import pandas as pd                       # easy to use 2D data structure
import numpy as np                        # contains a variety of numerical
                                          # methods and data structures
from collections import Counter           # structure that counts occurences
from smt.sampling_methods import LHS      # latin hypercube sampling function
#import pyDOE as pydoe                     # contains sampling methods
import os                                 # interacts with PATH to find files
import time                               # used for timer function
import tqdm                               # progress bar for iterables
import math                               # mathematical functions
import re                                 # regex bayBEEE
#import operator                           # adds named operators (i.e. add) in
                                          # order to use with mapping
alltime = time.time() # start timer for total run time
BIGN = 1.15e6

#import cubefinder                         # note that this is self-programmed
                                          # cubefinder finds a spherical and 
                                          # cubic point set around the origin
from cubefinder import cubicset
from cubefinder import sphereset
# cubefinder adds the cubicset and sphereset lists
# %% define target folder and load files ======================================
path = 'C:\\Users\\User\\Documents\\SOFC_Data_Science\\Data_OG'


def load_files(path):
    """
    param: path: full path to the folder containing all data files
    return: loaded: dict of all files loaded, with keys equal to the important 
            part of the file names: format "PCT195-A2-Anode"
    return: filenames: list of full file paths for each file
    """
    tic = time.time()
    filenames = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if '.mat' in file:
                filenames.append(os.path.join(r, file))

    print('files found: ')
    for f in filenames:
        print('.' + f[49:])

    # load all .mat files from target folder
    loaded={}
    for filename in filenames:
        fn = 'Amira_' + re.sub(r"[^\w\s\\:]", '_', filename)[50:]
        loaded.update({filename[50:-4] : loadmat(filename)[fn]})
    print(str(time.time()-tic) + ' seconds elapsed')
    return loaded, filenames
loaded,filenames = load_files(path)

# %% Load only as reduced button cell data ====================================
# this loads only the as-reduced anode file data [4]
# onefile = loaded[filenames[4][50:-4]]
filenum = 0
onefile = loaded[filenames[filenum][50:-4]]

# filenames[1] = 6.012 # could be Button cell, tested for 493 h?
# filenames[2] = 6.931
# filenames[3] = 8.257
print(filenames[filenum][50:-4])
# %% LHS Sampling =============================================================
tic = time.time()
matrixlimits = np.array([ [0, onefile.shape[0]-1], [0, onefile.shape[1]-1], [0, onefile.shape[2]-1]])
sampling = LHS(xlimits=matrixlimits)

samps = sampling(int(BIGN)) # control number of points generated
samps = np.round(samps,decimals=0).astype(int)
print('\n')
print('Cube sampled: ', str(time.time()-tic)[0:5], ' seconds elapsed.')
print('Sample size (points, dimensions):', samps.shape)
print('\n')



# %% Create search algorithm for TPBs'
def quicksearch(boundarycubes, lhspoint, filematrix, searchdirections):  
  
    """
    searchdirections = [(-1, 1, 1),
                        (0 , 1, 1),
                        (1 , 1, 1),
                        (-1, 0, 1),
                        (0 , 0, 1),
                        (1 , 0, 1),
                        (-1,-1, 1),
                        (0 ,-1, 1),
                        (1 ,-1, 1),
                        (-1, 1, 0),
                        (0 , 1, 0),
                        (1 , 1, 0),
                        (-1, 0, 0),
                        (0 , 0, 0),
                        (1 , 0, 0),
                        (-1,-1, 0),
                        (0 ,-1, 0),
                        (1 ,-1, 0),
                        (-1, 1,-1),
                        (0 , 1,-1),
                        (1 , 1,-1),
                        (-1, 0,-1),
                        (0 , 0,-1),
                        (1 , 0,-1),
                        (-1,-1,-1),
                        (0 ,-1,-1),
                        (1 ,-1,-1)]
    """
    
    
    # add a test point for every direction inside searchdirections
    #pointcube = [np.array(lhslist[i]) + np.array(searchdirections[i]) for i in range(len(searchdirections))]
    pointcube = [lhspoint + np.array(searchdirections[i]) for i in range(len(searchdirections))]
    pointcube = list(map(tuple, pointcube))
    pointcube = [poss for poss in pointcube if not(any(item < 0 for item in poss))] # remove all items with negative values
    pointcube = [poss for poss in pointcube if  not(any([poss[i] >= matrixlimits[i][1]-1 for i in [0,1,2]]))]
    
    # convert pointcubes to the corresponding phases
    phases = [filematrix[poss] for poss in pointcube] 
    phases = list(dict(Counter(phases)).keys())
    
    # check to ensure three phases occur inside the shape
    if len(phases) == 3:
        boundarycubes += 1
        
    return boundarycubes
# %% Run through random samples ===============================================
boundaryvols = 0
for point in tqdm.tqdm(samps): 
    boundaryvols = quicksearch(boundaryvols, tuple(point), onefile, sphereset)
    #boundaryvols = quicksearch(boundaryvols, tuple(point), onefile, cubicset)
    
    
# %% count up total number of cubes ===========================================
# We are aiming for 17.1 um^-2
# print(boundaryvols*0.05629/4350,'um^-2') # answer, if r = 1 and cubic
# print(boundaryvols*0.394943/4350,'um^-2') # answer, if r = 3 (s = 7) and cubic
print(boundaryvols*0.195446/4350,'um^-2') # answer, if r = 3 and spheric

print('Total time elapsed: ', time.time()-alltime, 'seconds')
print(str(BIGN) + ' points scanned.')

# to make this faster, we can increase the size of the cubes

# %% Results (assuming cubic, r=1) ============================================
# BIGN  # Time (s) # TBP V # dT/dN (WRT 1e5)
# 2e5   # 39       # 0.16  # .....
# 1e6   # 209      # 0.80  # 
# 2e6   # 428      # 1.60  # 
# 1e7   # 1909     # 7.98  # 

# %% Results (assuming cubic, r=3) ============================================
# BIGN  # Time (s) # TBP V # dT/dN (WRT 1e5)
# 2.5e5 # 603      # 9.912 # .....
# 5e5   # 1229     # 19.83 # 

# %% Results (assuming spheric, r=3) ==========================================
# BIGN  # Time (s) # TBP V # % of total volume covered
# 1e6   # 889      # 15.39 # 2.5e-3 %
# 1.1e6 # 1063     # 16.86 # 
# 1.15e6# 979      # 17.69 # 2.9e-3 % 
# 1.25e6# 1130     # 19.21 # 

# %% Results (assuming spheric, r=2) ==========================================
# BIGN  # Time (s) # TBP V # % of total volume covered
# 2e6   # 529      # 15.04 # 
# 2.5e6 # 1055*    # 18.80 # 
# 5e6   # 1991     # 37.69 # 

# %% Notes ====================================================================

# https://drive.google.com/drive/u/1/folders/0B7jcpLXoa_I9fmplOHBYMzVXV09fOTJJWWYxZjRQR2pmczlUeVViTjVfcFRrREhOQl90ck0
# sample volume (from H's data):
# 4350 um^3
# note that the data shows a volume of 125856000 units^3, so we know each 
# integer step in the matrix is 1/30.699 um [0.0325 um]
# total TPB (from H's data)
# 17.1 um^-2: WE WANT TO GET THIS VALUE

# diagonal = sqrt(3)*a = 1.7320508075688772 * a
# diag = math.sqrt(3) * (0.0325743509560572*(2r+1))
# diag(soherical) = diameter = 0.0325743509560572 * 2 * r


# LEAP OF FAITH: the A2 file SHOULD correspond to either the 72hr or 493hr tests














