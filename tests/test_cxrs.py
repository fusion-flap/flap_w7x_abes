#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:10:03 2020

@author: mvecsei
"""

import matplotlib.pyplot as plt
import os

import flap
import flap_w7x_abes


if __name__ == '__main__':

    # shotID = '20230314.032'
    shotID = '20240314.032'

    # flap_w7x_abes.register()
    # d=flap.get_data('W7X_ABES_CXRS',exp_id=shotID, name="QSI")
    
    config, spect_config, patchp_config = flap_w7x_abes.cxrs_util.read_fibre_config(exp_id=shotID)
    flap_w7x_abes.cxrs_util.plot_fibre_config(config)
    
    flap_w7x_abes.cxrs_util.plot_patchpanel_spectrometer_config(spect_config)
    
    flap_w7x_abes.cxrs_util.plot_patchpanel_optical_channel_config(spect_config, patchp_config)


    
    a = flap_w7x_abes.ShotSpatCalCXRS(shotID)
    a.generate_shotdata(options={'Plot':True, 'Overwrite':True})
    
    data = dict()
    index = 0
    for channel in a.data["Channel name"]:
        # print(f"{channel} \t {a.data['Device x'][index]} \t {a.data['Device y'][index]}")
        data[channel] = [a.data['Device x'][index], a.data['Device y'][index]]
        index += 1
        
    channels = sorted(data.keys())
    for channel in channels:
        print(f"{channel}\t{data[channel]}")