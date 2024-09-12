#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:30:10 2023

@author: apdcam
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.fft import fft, fftfreq

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def cleanup_chopped_data(data):
    data.shape = [data.shape[0]]
    # data.data = np.concatenate(data.data)
    data.del_coordinate("Rel. Time in int(Sample)")
    data.del_coordinate("Start Time in int(Sample)")
    data.del_coordinate("Rel. Sample in int(Sample)")
    data.del_coordinate("Start Sample in int(Sample)")
    data.del_coordinate("Interval(Sample)")
    data.del_coordinate("Interval(Sample) sample index")
    for coord in data.coordinate_names():
        data.get_coordinate_object(coord).dimension_list = [0]
        coord_obj = data.get_coordinate_object(coord)
        # coord_obj.shape = [coord_obj.shape[0]]
        # coord_obj.values = np.concatenate(coord_obj.values)
        coord_obj.start = [coord_obj.values[0]]
        coord_obj.step = [coord_obj.values[1]-coord_obj.values[0]]
        coord_obj.mode.equidistant = True
        coord_obj.values = None
    return data

def dataobj_fft(dataobject):
    N = len(dataobject.data)
    T = dataobject.coordinate("Time")[0][1]-dataobject.coordinate("Time")[0][0]
    xf = fftfreq(N, T)[:N//2]
    yf = 2.0/N*np.abs(fft(dataobject.data)[:N//2])
    return xf, yf
    

def test_W7X_data():
    exp_id="20181018.003"
    time_window_single = [4.25,4.38]
    time_window = [3.25,5]
    time_window_read = [time_window[0], time_window[1]+0.2]
    d_beam_on=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time', coordinates={'Time':time_window_read},
                         options={'State':{'Chop': 0, 'Defl': 0}, 'Start delay': 500, 'End delay': 0}, object_name='Beam_on')
    d_beam_off=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time', coordinates={'Time':time_window_read},
                         options={'State':{'Chop': 1, 'Defl': 0}, 'Start delay': 500, 'End delay': 0}, object_name='Beam_off')
    fig = plt.figure()
    for channel in range(40):       
        data = flap.get_data('W7X_ABES',
                             name=f'ABES-{channel+1}',
                             exp_id=exp_id,
                             coordinates={'Time':time_window_read},
                             object_name='ABES-DATA',
                             options={'Amplitude calibration': False,
                                      'Amplitude calib. path': "/data/cal/",
                                      'Scaling': "Volt",
                                      "Start delay": 500,
                                      "End delay": 0})
        data_on = data.slice_data(slicing={'Sample':d_beam_on})
        plt.subplot(8,5,channel+1)
        yf = 0
        for intervals in np.unique(data_on.coordinate("Interval(Sample)")[0]):
            data_on_curr = data_on.slice_data(slicing={"Interval(Sample)":intervals})
            flap_w7x_abes.add_absolute_time_sample(data_on_curr)
            data_on_curr = cleanup_chopped_data(data_on_curr)
            xf,yf_curr = dataobj_fft(data_on_curr)
            yf += yf_curr
        plt.plot(xf, yf*np.conj(yf))
        
        data_off = data.slice_data(slicing={'Sample':d_beam_off})
        yf = 0
        for intervals in np.unique(data_on.coordinate("Interval(Sample)")[0]):
            data_off_curr = data_off.slice_data(slicing={"Interval(Sample)":intervals})
            flap_w7x_abes.add_absolute_time_sample(data_off_curr)
            data_off_curr = cleanup_chopped_data(data_off_curr)
            xf,yf_curr = dataobj_fft(data_off_curr)
            yf += yf_curr
        plt.plot(xf, yf*np.conj(yf),label=f"ABES-{channel+1}", alpha=0.8)
        plt.yscale("log")
        plt.xscale("log")
        plt.legend(handlelength=0)
        plt.ylim([1e-8, 1e0])
        if channel%5 != 0:
            plt.yticks([])
        else:
            plt.yticks([1e-7, 1e-5, 1e-3, 1e-1])
    fig.text(0.5,0.04, "Frequency [Hz]", ha="center")
    fig.text(0.04,0.5, "Power Spectrum [a.u.]", ha="center", rotation="vertical")
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)

# Reading configuration file in the test directory
thisdir = os.path.dirname(os.path.realpath(__file__))
fn = os.path.join(thisdir,"w7x_config.cfg")
flap.config.read(file_name=fn)

test_W7X_data()