# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:02:34 2024

@author: Zoletnik
"""

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def plot_power(exp_ID,timerange=None,signals=['ABES-10','ABES-15','ABES-19', 'ABES-23'],
               datapath=None,background=False,fastchop=True,
               frange=None,fres=None,flog=None):
    """
    Plots the power spectrum of background corrected signals or background signals.

    Parameters
    ----------
    exp_ID : string
        The experiment ID.
    timerange : list of two floats or None, optional
        The time range. If None the whole data is used. The default is None.
    signals : list or string, optional
        The list of signals to process, ABES-xx. The default is ['ABES-10','ABES-15','ABES-19', 'ABES-23'].
    datapath : string or None, optional
        The data path. If None the one in flap_defaults.cfg will be used. The default is None.
    background : boolean, optional
        If True use background signal. The default is False.
    fastchop : boolean or None, optional
        If True average signals in a chopper time. 
        If False calculates spectra from.
        If None uses fastchop if the chopper period time is shorter thatn 100 microsec.
        The default is True.
    frange : list, optional
        Frequency range [Hz]. If None use setting in flap_defaults.cfg. The default is None.
    fres : float, optional
        frequency resolution [Hz].If None use setting in flap_defaults.cfg. The default is None.
    flog : boolean, optional
        If True use logarithmic frequency resolution. If None use setting in flap_defaults.cfg.
        The default is None.

    Raises
    ------
    NotImplementedError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    options = {}
    if (fres is not None):
        options['Resolution'] = fres
    if (frange is not None):
        options['Range'] = frange
    if (flog is not None):
        options['Logarithmic'] = flog
   
    d_beam_on = flap.get_data('W7X_ABES',
                              exp_id=exp_ID,
                              name='Chopper_time',
                              coordinate = {'Time':timerange},
                              options={'State': {'Chop': 0, 'Defl': 0},
                                       'Start': 0, 'End': 0}
                              )
    d_beam_off = flap.get_data('W7X_ABES',
                               exp_id=exp_ID,
                               name='Chopper_time',
                               coordinate = {'Time':timerange},
                               options={'State': {'Chop': 1, 'Defl': 0},
                                        'Start': 0, 'End': 0}
                               )
    chopper_mode = d_beam_on.info['Chopper mode']
    on1, on2, on3 = d_beam_on.coordinate('Time')
    off1, off2, off3 = d_beam_off.coordinate('Time')
    beam_on_time = on3[1]-on2[1]
    beam_off_time = off3[1]-off2[1]
    period_time = beam_on_time + beam_off_time
    if (fastchop is None):
        if (period_time < 1e-4):
            _fastchop = True
        else:
            _fastchop = False
    else:
        _fastchop = fastchop
        
    if (_fastchop):
        if (period_time > 1e-4):
            print("Slow chopper! Not suitable for calculating spectra from background corrected signals.")
        d = flap_w7x_abes.chopped_signals(exp_ID,timerange=timerange,signals=signals,datapath=datapath,background=background)      
        p = d.apsd(coordinate='Time',options=options)
    else:
        if (background):
            d_chop = d_beam_on 
        else:
            d_chop = d_beam_off
        p = d.apsd(coordinate='Time',intervals={'Time':d_chop},options=options)
 
    options = {'Log x':True,'Log y':True}
    p.plot(axes='Frequency',options=options)

        
