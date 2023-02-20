#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:35:45 2022

@author: apdcam
"""
import flap
import flap_apdcam
import os
import matplotlib.pyplot as plt
import numpy as np
import math

flap_apdcam.register()

def show_all(shot,time=[0,1]):
    datapath='/data'
    d=flap.get_data('APDCAM',
                    name='ADC*',
                    coordinates={'Time':time},
                    options={'Datapath':os.path.join(datapath,shot),
                             'Camera type':'APDCAM-10G_FC',
                             'Scaling':'Volt'
                             }
                    )
    d = d.filter_data(coordinate='Time',options={'Type':'Int','Tau':1e-3})
    d = d.slice_data(slicing={'Time':flap.Intervals(0,3e-4,step=3e-4)},
                     summing={'Interval(Time) sample index':'Mean'})   
    d1 = d.slice_data(slicing={'ADC Channel':flap.Intervals(1,8,step=8)})

    plt.figure()
    d1.plot(plot_type='grid xy',axes=['Interval(ADC Channel)','Interval(ADC Channel) sample index','Time'],
            options={'Y range':[0,2]})    
    plt.suptitle(shot)

def show_all_spectra(shot,timerange=None):

    datapath='/data'
    if (timerange is None):
        tr = None
    else:
        tr= {'Time':timerange}
    d=flap.get_data('APDCAM',
                    name='ADC*',
                    coordinates=tr,
                    options={'Datapath':os.path.join(datapath,shot),
                             'Camera type':'APDCAM-10G_FC',
                             'Scaling':'Volt'
                             }
                    )
    d1 = d.slice_data(slicing={'ADC Channel':flap.Intervals(1,8,step=8)})
    print('Calculating power spectra...')
    d1p = d1.apsd(coordinate='Time',options={'Resolution':1e1,
                                      'Log':True,
                                      'Range':[1e3,1e6]
                                      }
                  )
    plt.rcParams['axes.labelsize'] = 6
    plt.rcParams['axes.titlesize'] = 6
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    plt.rcParams['axes.titlepad'] = -10
    
    d1p.plot(plot_type='grid xy',axes=['Interval(ADC Channel)','Interval(ADC Channel) sample index','Frequency'],
                  options={'Log x':True,'Log y':True,'Y range':[2e-6,1e-4]}
                  )
    plt.suptitle(shot)
    
# plt.close('all')    
# show_all_spectra('T20221214.001')
# show_all_spectra('T20221214.002')
# show_all_spectra('T20221214.004')
# show_all_spectra('T20221214.003')
