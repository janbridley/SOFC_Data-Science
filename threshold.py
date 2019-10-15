# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:58:28 2019

@author: User
"""
import numpy as np


def threshold(img,highThreshold=2,lowThreshold=1.5):
    # this function seperates out all of the pixels 
    # high/low of 2/1.5 seperates out all phases
    # high/low of 3/2.5 finds only the phase labelled '3'
    
    #highThreshold = 2
    #lowThreshold = 1.5#highThreshold * lowThresholdRatio;
    
    M, N = img.shape
    res = np.zeros((M,N), dtype=np.uint8)
    
    weak = np.uint8(1)
    strong = np.uint8(2)
    
    strong_i, strong_j = np.where(img >= highThreshold)
    zeros_i, zeros_j = np.where(img < lowThreshold)
    
    weak_i, weak_j = np.where((img <= highThreshold) & (img >= lowThreshold))
    
    res[strong_i, strong_j] = strong
    res[weak_i, weak_j] = weak
    
    return res