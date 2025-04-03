#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:10:03 2020

@author: mvecsei
"""

import matplotlib.pyplot as plt
import os
import datetime
import numpy as np

import flap
import flap_w7x_abes
import flap_w7x_webapi as webapi
webapi.register()

class CXRSPlotter():
    def __init__(self, shotID, canvas=None):
        self.shotID=shotID
        self.curr_sample = 0
        self.canvas = canvas
        self.channelrange=1
        self.plot_cxrs()
        

    def get_data(self):
        datapath = '/data'
        flap_w7x_abes.register()

        self.cxrs = flap.get_data('W7X_ABES_CXRS', exp_id=self.shotID,
                                  name="QSI_CXRS")

        self.channelrange = self.cxrs.data.shape[1]
    
    def plot_cxrs(self, setto=None, moveby=None):
    
        if hasattr(self, "cxrs") is False:
            self.get_data()
        
        if setto is not None:
            self.curr_sample = int((setto-1)/100.0*self.channelrange)
        
        if moveby is not None:
            self.curr_sample += moveby
            
        if hasattr(self, "fig") is False:
            if self.canvas is None:
                self.fig, self.ax = plt.subplots()
            else:
                self.fig = self.canvas.fig
                self.ax = self.canvas.fig.add_subplot(1,1,1)
              
        self.ax.cla()    
        self.ax.set_title(f"{self.shotID}  background subtracted CXRS at "+\
                         f"{int(self.cxrs.get_coordinate_object('Device R').values[self.curr_sample]*1000)/1000.0}m")
        beam_on = self.cxrs.data[::2,self.curr_sample,:]
        beam_off = self.cxrs.data[1::2,self.curr_sample,:]
        signal = beam_on[:beam_off.shape[0],:]-beam_off
        signal -= np.min(np.abs(signal))
        norm_signal = signal.transpose()
        newmat = (np.ones(norm_signal.shape)*np.arange(norm_signal.shape[1])*np.max(norm_signal)/4).transpose()
        norm_signal = norm_signal+newmat.transpose()
        # norm_signal =+ np.arange(norm_signal.shape[1])*np.max(norm_signal)*np.ones(norm_signal.shape)
        self.ax.plot(self.cxrs.get_coordinate_object("Wavelength").values[self.curr_sample,:],norm_signal, alpha = 0.3, color= "tab:blue")
        self.ax.set_yticks(newmat[:,0],
                           (self.cxrs.get_coordinate_object("Time").values[1::2]*100).astype(int)/100)
        self.ax.set_ylabel("Time [s]")
        self.ax.set_xlabel("Wavelength [nm]")

        
        
        if self.canvas is None:
            plt.show()
        else:
            self.fig.canvas.draw_idle()
            self.fig.canvas.flush_events()
        

if __name__ == "__main__":
    shotID = '20240926.027'
    plot_signal = True
    if plot_signal == True:
        plotter = CXRSPlotter(shotID)
        plotter.plot_cxrs()
        plotter.plot_cxrs(setto=10)
    
    plot_patchpanel = False
    if plot_patchpanel is True:
        config, spect_config, patchp_config = flap_w7x_abes.cxrs_util.read_fibre_config(exp_id=shotID)
        
        
        flap_w7x_abes.cxrs_util.plot_fibre_config(config)
        flap_w7x_abes.cxrs_util.plot_patchpanel_spectrometer_config(spect_config)
        
        flap_w7x_abes.cxrs_util.plot_patchpanel_optical_channel_config(spect_config, patchp_config)
    

    get_spatcal = False
    if get_spatcal is True:    
        a = flap_w7x_abes.ShotSpatCalCXRS(shotID)
        a.generate_shotdata(options={'Plot':True, 'Overwrite':True})
    
    read_spatcal = False
    if read_spatcal is True:
        a = flap_w7x_abes.ShotSpatCalCXRS(shotID)
        a.read()
        
        #printing out the results
        data = {}
        index = 0
        for channel in a.data["Channel name"]:
            # print(f"{channel} \t {a.data['Device x'][index]} \t {a.data['Device y'][index]}")
            data[channel] = [a.data['Device x'][index], a.data['Device y'][index], np.sqrt(a.data['Device x'][index]**2+a.data['Device y'][index]**2)]
            index += 1
            
        channels = sorted(data.keys())
        for channel in channels:
            print(f"{channel}\t{data[channel]}")
            
    get_cxrs_data = False
    if get_cxrs_data is True:
        shotID= '20240926.027'
        flap.config.read()
        flap_w7x_abes.register()
        
        #test webapi
        try:
            d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              exp_id=shotID,
                              options={'Scale Time': True,
                                       'Cache Data': False,
                                       'Check Time Equidistant': True},
                              object_name='QSI CXRS data',
                              coordinates={'Time': [2, 3]})
        except Exception as e:
            raise ValueError("something is wrong with the flap_w7x_webapi:" + str(e))
            d = flap.get_data('W7X_ABES_CXRS', name='CXRS_TEST',
                                      exp_id=shotID,
                                      object_name=f"CXRS_{shotID}")

        timestep = 1.2 #s
        d_slice = d.slice_data(slicing={"Time":timestep})
        
        plt.figure(figsize=[8,12])
        plt.suptitle(f"QSI CXRS data for {shotID}")
        for index, device_r in enumerate(d.get_coordinate_object("Device R").values):
            plt.subplot(len(d.get_coordinate_object("Device R").values),1,index+1)
            d_slice.slice_data(slicing={"Device R":device_r}).plot(axes="Wavelength")
            plt.ylabel(f"R={int(device_r*100)/100.0}m")
            plt.title("")
            if index < len(d.get_coordinate_object("Device R").values)-1:
                plt.xlabel("")
                plt.xticks([])
            plt.yticks([])
            plt.tight_layout()