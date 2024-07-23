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

    shotID = '20230314.032'
    # shotID = '20240314.032'

    # ID = [1696580349815]
    # d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
    #                           options={'Scale Time': True,
    #                                   'Cache Data': False,
    #                                   'Time Start':datetime.datetime(2023, 10, 6, 8, 00, 00),
    #                                   'Time Stop': datetime.datetime(2023, 10, 7, 10, 25, 00)},
    #                           object_name="ROIs")
    # ROIs = np.arange(54)+1
    # plt.pcolormesh(d.data[-1,:])
    # plt.pcolormesh(d.data[-2,:], alpha=0.5, cmap = "gray")

    # plt.title(ROIs[d.data.shape[0]-2])
    # plt.pause(5)
    # plt.close()
    
    d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              options={'Scale Time': True,
                                      'Cache Data': False,
                                      'Timestamp': 1697459453679000000},
                              object_name="ROIs")
    asdas
    
    ID = [1696580349815]
    d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              options={'Scale Time': True,
                                      'Cache Data': False,
                                      'Time Start':datetime.datetime(2023, 10, 6, 19, 00, 00),
                                      'Time Stop': datetime.datetime(2023, 10, 7, 10, 30, 00)},
                              object_name="ROIs")
    
    # plt.imshow(d.data[-1,:])
    # plt.pcolormesh(d.data[-2,:], alpha=0.5, cmap = "gray")
    fibers = dict()
    # fibers["1"] = [65,101]
    # fibers["4"] = [110,150]
    # fibers["7"] = [152,196]
    # fibers["10"] = [202,247]
    # fibers["13"] = [249,296]
    # fibers["16"] = [300,338]
    # fibers["19"] = [340,380]
    # fibers["22"] = [382,420]
    # fibers["25"] = [422,465]
    # fibers["28"] = [480,535]
    # fibers["31"] = [540,580]
    # fibers["34"] = [591,638]
    # fibers["37"] = [640,680]
    # fibers["40"] = [690,730]
    # fibers["43"] = [742,780]
    # fibers["46"] = [782,820]
    # fibers["2"] = [67,115]
    # fibers["5"] = [120,170]
    # fibers["8"] = [172,210]
    # fibers["11"] = [215,256]
    # fibers["14"] = [260,305]
    # fibers["17"] = [310,350]
    # fibers["20"] = [352,395]
    # fibers["23"] = [400,440]
    # fibers["26"] = [445,500]
    # fibers["29"] = [510,555]
    # fibers["32"] = [560,600]
    # fibers["35"] = [610,655]
    # fibers["38"] = [660,700]
    # fibers["41"] = [710,750]
    # fibers["49"] = [835,882]
    # fibers["52"] = [890, 930]
    # fibers["9"] = [193, 229]
    # fibers["12"] = [230,277]
    # fibers["15"] = [280,319]
    # fibers["18"] = [320,362]
    # fibers["21"] = [370,410]
    # fibers["24"] = [412,450]
    # fibers["27"] = [451,512]
    # fibers["30"] = [513,580]
    # fibers["33"] = [581,623]
    # fibers["36"] = [624,660]
    # fibers["39"] = [661,712]
    # fibers["42"] = [713,756]
    # fibers["44"] = [757,792]
    # fibers["47"] = [793,847]
    # fibers["50"] = [848,900]
    # fibers["54"] = [901,958]
    fibers["3"] = [94,135]
    fibers["6"] = [136,178]
    fibers["45"] = [769,810]
    fibers["48"] = [810,860]
    fibers["51"] = [861,908]
    fibers["53"] = [909,941]
    
    index = 1
    for key in fibers:
        plt.subplot(6,1,index)
        data_to_plot = np.sum(d.data[-1,fibers[key][0]:fibers[key][1],:],axis=0)
        plt.plot(data_to_plot-np.min(data_to_plot),label=key)
        index += 1
        plt.legend()


    # plt.title(ROIs[d.data.shape[0]-2])
    plt.legend()
    # plt.pause(10)
    # plt.close()
    
    
    
    # flap_w7x_abes.register()
    # d=flap.get_data('W7X_ABES_CXRS',exp_id=shotID, name="QSI")
    

    config, spect_config, patchp_config = flap_w7x_abes.cxrs_util.read_fibre_config(exp_id=shotID)
    asda

    flap_w7x_abes.cxrs_util.plot_fibre_config(config)
    flap_w7x_abes.cxrs_util.plot_patchpanel_spectrometer_config(spect_config)
    
    flap_w7x_abes.cxrs_util.plot_patchpanel_optical_channel_config(spect_config, patchp_config)


    
    a = flap_w7x_abes.ShotSpatCalCXRS(shotID)
    a.generate_shotdata(options={'Plot':True, 'Overwrite':True})
    
    data = dict()
    index = 0
    for channel in a.data["Channel name"]:
        # print(f"{channel} \t {a.data['Device x'][index]} \t {a.data['Device y'][index]}")
        data[channel] = [a.data['Device x'][index], a.data['Device y'][index], np.sqrt(a.data['Device x'][index]**2+a.data['Device y'][index]**2)]
        index += 1
        
    channels = sorted(data.keys())
    for channel in channels:
        print(f"{channel}\t{data[channel]}")

