#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:38:14 2024

@author: zoletnik
"""

import matplotlib.pyplot as plt

import flap
import flap_w7x_abes
flap_w7x_abes.register()


def plot_chopper(exp_ID, signal='ABES-15', timerange=None, resample=""):
    plt.close('all')

    d_beam_on = flap.get_data('W7X_ABES',
                              exp_id=exp_ID,
                              name='Chopper_time',
                              options={'State': {'Chop': 0, 'Defl': 0},
                                       'Start': 0, 'End': 0}
                              )
    d_beam_off = flap.get_data('W7X_ABES',
                               exp_id=exp_ID,
                               name='Chopper_time',
                               options={'State': {'Chop': 1, 'Defl': 0},
                                        'Start': 0, 'End': 0}
                               )
    chopper_mode = d_beam_on.info['Chopper mode']
    on1, on2, on3 = d_beam_on.coordinate('Time')
    off1, off2, off3 = d_beam_off.coordinate('Time')
    beam_on_time = on3[1]-on2[1]
    beam_off_time = off3[1]-off2[1]
    period_time = beam_on_time + beam_off_time
    if (period_time > 3e-3):
        chop_str = "{:3.0f}-{:3.0f}[ms]".format(
            beam_on_time * 1e3, beam_off_time * 1e3)
    else:
        chop_str = "{:3.0f}-{:3.0f}[us]".format(
            beam_on_time * 1e6, beam_off_time * 1e6)
    txt = '{:s} ... Chopper:{:6s}({:s})'.format(exp_ID, chopper_mode, chop_str)
    print(txt)

    options = {}
    if (period_time > 0.01):
        on_start = 1000
        on_end = -1000
        off_start = 1000
        off_end = -1000
        options['Resample'] = 1e4
    else:
        on_start = 0
        on_end = 0
        off_start = 0
        off_end = 0
        options['Resample'] = None

    if (resample != ""):
        options['Resample'] = resample
    flap_w7x_abes.test_chopper_timing(exp_id=exp_ID,
                                      timerange=timerange,
                                      signal=signal,
                                      resample=None,
                                      x_axis='Time',
                                      start_shift=on_start,
                                      end_shift=on_end
                                      )

# plt.close('all')
# plot_chopper('20230315.025',signal='ABES-24',timerange=[1,1.001],resample=None)