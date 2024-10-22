#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:00:01 2024

@author: zoletnik
"""
import flap
import flap_w7x_abes

flap_w7x_abes.register()

def chopper_parameters(exp_ID,datapath=None):
    """
    Returns the beam chopper parameters

    Parameters
    ----------
    exp_ID : string
        The experiment ID.
    datapath : string or None, optional
        The data path. If None the one in flap_defaults.cfg will be used. The default is None.

    Returns
    -------
    chopper_mode : text
        The chopper mode: 'Camera','Timed'.
    beam_on_time : float
        The time [s] for which the beam is on.
    beam_off_time : float
        The time [s] for which the beam is off.
    period_time : float
        The period time [s] of the chopping.
    d_beam_on : DataObject
        The data object for the beam on times during the whole discharge.
    d_beam_off : DataObject
        The data object for the beam off times during the whole discharge.
    """
    
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
    beam_on_time = on3[0]-on2[0]
    beam_off_time = off3[0]-off2[0]
    period_time = beam_on_time + beam_off_time
    if (period_time > 3e-3):
        chop_str = "{:3.0f}-{:3.0f}[ms]".format(
            beam_on_time * 1e3, beam_off_time * 1e3)
    else:
        chop_str = "{:3.0f}-{:3.0f}[us]".format(
            beam_on_time * 1e6, beam_off_time * 1e6)
    txt = '{:s} ... Chopper:{:6s}({:s})'.format(exp_ID, chopper_mode, chop_str)

    return chopper_mode,beam_on_time,beam_off_time,period_time,d_beam_on,d_beam_off