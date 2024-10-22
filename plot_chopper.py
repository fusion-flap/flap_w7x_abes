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


def plot_chopper(exp_ID, signal='ABES-15', timerange=None, resample="", datapath=None):
    """
    Plot a beam signal, mark chopper times and overplot the chopper-period averaged signal.

    Parameters
    ----------
    exp_ID : string
        The experiment ID.
    datapath : string or None, optional
        The data path. If None the one in flap_defaults.cfg will be used. The default is None.
    signal : string, optional
        The signal name. The default is 'ABES-15'.
    timerange : List of two floats or None, optional
        Time range to process. The default is None.
    resample : float of "", optional
        A float number indecates a resample frequency for reading the data. This us useful for flow chopper. 
        The default is "".
 
    Returns
    -------
    None.

    """
    
    plt.close('all')

    chopper_mode,beam_on_time,beam_off_time,period_time,d_beam_on,d_beam_off = flap_w7x_abes.chopper_parameters(exp_ID,datapath=datapath)

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
                                      resample=options['Resample'],
                                      x_axis='Time',
                                      start_shift=on_start,
                                      end_shift=on_end
                                      )

# plt.close('all')
# plot_chopper('20230315.025',signal='ABES-24',timerange=[1,1.001],resample=None)