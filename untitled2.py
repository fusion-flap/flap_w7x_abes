#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 14:28:09 2023

@author: apdcam
"""

from  flap_w7x_abes.cmos_as_dataobject import *
import flap_w7x_abes
import flap
from matplotlib import pyplot as plt

flap_w7x_abes.register()

if __name__ == "__main__":
    shotID = "20230307.047"
    config_file = "/repos/flap_modules/flap_spade/flap_defaults.cfg"
    flap.config.read(file_name=config_file)
    cmos = w7x_abes_cmos_get_data(exp_id=shotID)
    # o={'State':{'Chop': 0, 'Defl': 0}}
    # d_beam_on=flap.get_data('W7X_ABES',
    #                         exp_id=shotID,
    #                         name='Chopper_time',
    #                         coordinates={'Time':[np.min(cmos.coordinate('Time')[0]), np.max(cmos.coordinate('Time')[0])]},
    #                         options=o,\
    #                         object_name='Beam_on',
    #                         )
    # times = d_beam_on.coordinate('Time')[1]/2+d_beam_on.coordinate('Time')[2]/2
    # d = cmos.slice_data(slicing={'Time':times})
    print("ok")
    cmos_on = cmos.get_chopstate()
    cmos_off = cmos.get_chopstate(chop=1)


    for timeindex in range(cmos_on.data.shape[0]):
    # for timeindex in range(109,119):
        plt.clf()
        plt.imshow(cmos_on.data[timeindex,:,:]-cmos_off.data[timeindex,:,:])
        plt.show()
        plt.pause(0.1)
