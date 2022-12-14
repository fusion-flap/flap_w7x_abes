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
            
class CMOS2021:
    
    def __init__(self):
    
        try:
            currdir = os.path.dirname(os.path.abspath(__file__))
        except NameError:
            currdir = os.path.abspath(os.getcwd())
            
        currdir = os.path.dirname(currdir)
        self.geomfolder = os.path.join(currdir, 'spatcal', '2021', 'Geometry')
        
        self.points_cmos = np.array([])
    
    def start(self):
        plt.ion()
        f = plt.figure(figsize=[8,5])
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 6, 1])
        # Adding the 'buttons'
        ax = plt.subplot(gs[0])
        buttons = OrderedDict()
        buttons['Save list'] = [True, self.save_list]
        buttons['Delete Last \nPoint'] = [True, self.delete_last_selected_point]
        buttons['Add Point'] = [False, None]
        self.clickbuttons = ClickButtonList(ax, buttons)
        
        calib_image = os.path.join(self.geomfolder, 'calib_image.png')
        self.calib_image = np.asarray(imread(calib_image))
        plt.subplot(gs[1])
        self.image_plot = plt.imshow(self.calib_image)
        self.image_plot.figure.canvas.mpl_connect('button_press_event', self.add_to_selected_points)
    
    def plot_points(self):
        calib_image = os.path.join(self.geomfolder, 'calib_image.png')
        self.calib_image = np.asarray(imread(calib_image))
        self.image_plot = plt.imshow(self.calib_image)
        
        with open(os.path.join(self.geomfolder, 'calib_image_points.dat'), "r") as f:
            all_data = f.readlines()
        self.points_cmos = list()
        for point in all_data:
            self.points_cmos += [[float(point.split("\n")[0].split("\t")[0]), float(point.split("\n")[0].split("\t")[1])]]
            plt.scatter(self.points_cmos[-1][0], self.points_cmos[-1][1], color='#FF0000', marker='x')
            plt.text(self.points_cmos[-1][0]+7, self.points_cmos[-1][1]+7, len(self.points_cmos), color='red')
        self.points_cmos = np.asarray(self.points_cmos)

    def add_to_selected_points(self, event):
        ''' Add CMOS coordinate of clicked point to list of selected points
        '''
        axes = self.image_plot.axes
        if (event.inaxes != axes) or (self.clickbuttons.buttons['Add Point'].active is False):
            return
        self.points_cmos = np.array(list(self.points_cmos)+[[event.xdata, event.ydata]])
        scatterplot = axes.scatter(self.points_cmos[:, 0], self.points_cmos[:, 1], color='#FF0000', marker='x')
        axes.draw_artist(scatterplot)
        text = axes.text(event.xdata+7, event.ydata+7, len(self.points_cmos), color='red')
        axes.draw_artist(text)
        canvas = self.image_plot.figure.canvas
        canvas.blit(axes.bbox)
        
    def delete_last_selected_point(self):
        ''' Deletes the last point from the list of selected points
        '''
        if len(self.points_cmos) == 0:
            return
        self.points_cmos = self.points_cmos[:-1]

        axes = self.image_plot.axes
        axes.cla()
        canvas = self.image_plot.figure.canvas
        self.image_plot = axes.imshow(self.calib_image)
        axes.draw_artist(self.image_plot)
        pointindex = 0
        for point in self.points_cmos:
            scatterplot = axes.scatter(point[0], point[1], color='#FF0000', marker='x')
            axes.draw_artist(scatterplot)
            text = axes.text(point[0]+7, point[1]+7, pointindex+1, color='red')
            axes.draw_artist(text)
            pointindex += 1
        canvas.blit(axes.bbox)
    
    def save_list(self):
        with open(os.path.join(self.geomfolder, 'calib_image_points.dat'), "w") as f:
            for point in self.points_cmos:
                f.writelines(f"{point[0]}\t{point[1]}\n")

class ClickButtonList:
    def __init__(self, axes, button_dict):
        num_of_buttons = len(button_dict.keys())
        buttons = axes.barh(np.arange(num_of_buttons), np.ones(num_of_buttons), height=0.5)
        self.buttons = []
        index = 0
        self.buttons = OrderedDict()
        for key in button_dict:
            self.buttons[key] = Button(buttons[index], key,
                                       single_event=button_dict[key][0],
                                       connected_function=button_dict[key][1])
            index += 1
        self.active_button_id = None

class Button:
    def __init__(self, button, button_name, single_event=True, connected_function=None):
        self.button = button
        self.button_name = button_name
        self.connected_function = connected_function
        self.active = False
        self.single_event = single_event
        self.connect()
        plt.pause(0.1)
        self.draw_button(state='off')

    def connect(self):
        self.cidpress = self.button.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        if self.single_event is True:
            self.cidpress = self.button.figure.canvas.mpl_connect(
                'button_release_event', self.on_release)

    def on_press(self, event):
        if event.inaxes != self.button.axes:
            return
        contains, attrd = self.button.contains(event)
        if contains is False:
            self.active = False
            self.draw_button(state='off')
        else:
            if self.active is False:
                self.active = True
                self.draw_button(state='on')
            else:
                self.active = False
                self.draw_button(state='off')

    def on_release(self, event):
        if self.active is False:
            return
        self.connected_function()
        self.active = False
        self.draw_button(state='off')

    def draw_button(self, state='off'):
        if state == 'off':
            color = '#1f77b4'
        else:
            color = '#800000'
        self.button.set_color(color)
        canvas = self.button.figure.canvas
        axes = self.button.axes
        bbox = axes.get_window_extent().transformed(self.button.figure.dpi_scale_trans.inverted())
        width = int(bbox.width * self.button.figure.dpi/10)
        text = self.button.axes.text(self.button.get_x()+self.button.get_width()/2,
                                     self.button.get_y()+self.button.get_height()/2,
                                     self.button_name, ha='center', va='center', fontsize=width, color='white')
        axes.draw_artist(self.button)
        axes.draw_artist(text)
        axes.axis('off')
        canvas.blit(axes.bbox)

    
# merge1821apdcam()
calibrate = CMOS2021()
# CMOS2021.start()
calibrate.plot_points()
