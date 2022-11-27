#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:59:52 2022

@author: mvecsei
"""

import os.path
import h5py
import numpy as np
from functools import partial
from scipy.optimize import minimize
from scipy.ndimage import median_filter
from functools import partial
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.image import imread
from collections import OrderedDict
import flap
import flap_w7x_abes


def merge1821apdcam():
    try:
        currdir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        currdir = os.path.abspath(os.getcwd())
        
    currdir = os.path.dirname(currdir)
    
    dir2018 = os.path.join(currdir, 'spatcal', '2018')
    fileID = np.arange(56)+4
    padded_fileID = ["{:08d}".format(number) for number in fileID]
    all_files = os.listdir(dir2018)
    rel_files = np.sort(np.asarray([filename for filename in all_files if (filename[:8] in padded_fileID)]))
    
    fileindex = 0
    for file in rel_files:
        image = np.asarray(imread(os.path.join(dir2018, file)))
        image[np.where(image<np.max(image)*0.2)] = 0
        # image[np.where(image<np.median(image)*2)] = 0
        if fileindex == 0:
            sumimage2018 = image/np.max(image)
        else:
            sumimage2018 = np.abs(image/np.max(image)-sumimage2018)
        fileindex += 1
            
    dir2021 = os.path.join(currdir, 'spatcal', '2021')
    fileID = ["16","17", "18", "19", "20", "21", "22"]
    all_files = os.listdir(dir2021)
   
    rel_files = np.sort(np.asarray([filename for filename in all_files if (filename.split("_")[0] in fileID)]))
    
    fileindex = 0
    for file in rel_files:
        image = np.asarray(imread(os.path.join(dir2021, file)))
        image[np.where(image<np.max(image)*0.5)] = 0
        # image[np.where(image<np.median(image)*2)] = 0
        if fileindex == 0:
            sumimage2021 = image/np.max(image)
        else:
            sumimage2021 += image/np.max(image)
        fileindex += 1
        
    image = np.load("202211_gasshot.npy")
    
    merged = np.zeros([np.shape(sumimage2021)[0], np.shape(sumimage2021)[1], 3])
    merged[:,:,0] = np.fliplr(image/np.max(image)*2)
    merged[:,:,1] = np.fliplr(sumimage2021/np.max(sumimage2021))
    merged[:,:,2] = np.fliplr(sumimage2018/np.max(sumimage2018)*255)
    plt.imshow(merged)

    plt.figure()

    merged[:,:,2]*=0
    plt.imshow(merged)
    filename = os.path.join(currdir, "spatcal", "2021",
                            'Geometry', 'apdcam_to_cmos.hdf5')
    apdcam_to_cmos = flap.load(filename)
    for channel in apdcam_to_cmos['APDCAM channel centers'].keys():
        try:
            plt.scatter(apdcam_to_cmos['APDCAM channel centers'][channel][1][2],
                        apdcam_to_cmos['APDCAM channel centers'][channel][1][3],
                        color='blue')
            plt.text(apdcam_to_cmos['APDCAM channel centers'][channel][1][2]+5,
                        apdcam_to_cmos['APDCAM channel centers'][channel][1][3]+5, 
                        channel, color='blue')
        except:
            plt.scatter(apdcam_to_cmos['APDCAM channel centers'][channel][1][0][2],
                        apdcam_to_cmos['APDCAM channel centers'][channel][1][0][3],
                        color='tab:blue')
            plt.text(apdcam_to_cmos['APDCAM channel centers'][channel][1][0][2]+7,
                        apdcam_to_cmos['APDCAM channel centers'][channel][1][0][3]+7, 
                        channel, color='blue')
    
merge1821apdcam()

