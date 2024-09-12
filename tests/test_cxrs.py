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


if __name__ == '__main__':

    shotID = '20240314.032'
    
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
        shotID= '20230330.028'
        flap.config.read()
        flap_w7x_abes.register()
        
        #test webapi
        try:
            d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              exp_id=shotID,
                              options={'Scale Time': True,
                                       'Cache Data': False},
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
            plt.tight_layout()