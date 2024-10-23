# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:02:34 2024

@author: Zoletnik
"""
import copy
import numpy as np
import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def plot_power(exp_ID,timerange=None,signals=['ABES-10','ABES-15','ABES-19', 'ABES-23'],
               datapath=None,resample="",signal_type='raw',beam_on=None,
               frange=None,fres=None,flog=None,interval_n=None,
               plot_type='xy',options=None,plot_id=None):
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
    resample : float of "", optional
        A float number indecates a resample frequency for reading the data. This us useful for flow chopper. 
        The default is "".
    signal_type : string, optional
        "raw": Now correction for background. 
        "beam": Background corrected beam signal.
        "background": The background signal.
        The default is 'raw'.
    beam_on : boolean or None, optional
        If True calculate power spectra in each chopper period when beam is on and average.
        If False calculates power spectra in each chopper period when beam is on and average.
        If None use the whole time interval for power calculation.
        The default is None.
    frange : list, optional
        Frequency range [Hz]. If None use setting in flap_defaults.cfg. The default is None.
    fres : float, optional
        frequency resolution [Hz].If None use setting in flap_defaults.cfg. The default is None.
    flog : boolean, optional
        If True use logarithmic frequency resolution. If None use default setting in flap_defaults.cfg.
        The default is None.
    interval_n: int or None
        If True use this value for the "Interval_n" options on apsd.
        If None use default setting in flap_defaults.cfg.
    options: dict
        Plot options, see flap.plot
    plot_id: flap.plotID
        The plot ID for overplotting a flap plot.

    Raises
    ------
    NotImplementedError
        DESCRIPTION.

    Returns
    -------
    plot_id:
        The plot id returned by flap.plot()

    """
    apsd_options = {}
    if (fres is not None):
        apsd_options['Resolution'] = fres
    if (frange is not None):
        apsd_options['Range'] = frange
    if (flog is not None):
        apsd_options['Logarithmic'] = flog
    if (interval_n is not None):
        apsd_options['Interval_n'] = interval_n

    if (beam_on is not None):
        if (signal_type !='raw'):
            raise ValueError("Cannot use chopper for interval selection for chopper-averaged signals.")
        
    chopper_mode,beam_on_time,beam_off_time,period_time,d_beam_on,d_beam_off = flap_w7x_abes.chopper_parameters(exp_ID,datapath=datapath)

    d = flap_w7x_abes.get_clean_abes(exp_ID,
                                     signals=signals,
                                     datapath=datapath,
                                     resample=resample,
                                     signal_type=signal_type,
                                     timerange=timerange
                                     )

    # d_beam_on = flap.get_data('W7X_ABES',
    #                           exp_id=exp_ID,
    #                           name='Chopper_time',
    #                           coordinates = {'Time':timerange},
    #                           options={'State': {'Chop': 0, 'Defl': 0},
    #                                    'Start': 0, 'End': 0}
    #                           )
    # d_beam_off = flap.get_data('W7X_ABES',
    #                            exp_id=exp_ID,
    #                            name='Chopper_time',
    #                            coordinates = {'Time':timerange},
    #                            options={'State': {'Chop': 1, 'Defl': 0},
    #                                     'Start': 0, 'End': 0}
    #                            )
    # chopper_mode = d_beam_on.info['Chopper mode']
    # on1, on2, on3 = d_beam_on.coordinate('Time')
    # off1, off2, off3 = d_beam_off.coordinate('Time')
    # beam_on_time = on3[0]-on2[0]
    # beam_off_time = off3[0]-off2[0]
    # period_time = beam_on_time + beam_off_time
        
    if (beam_on is None):
        p = d.apsd(coordinate='Time',options=apsd_options)
    else:
        if (beam_on):
            d_chop = d_beam_on
        else:
            d_chop = d_beam_off
        p = d.apsd(coordinate='Time',intervals={'Time':d_chop},options=apsd_options)
 
    if (options is not None):
        _plot_options = copy.deepcopy(options)
    else:
        _plot_options = {}
        
    try:
        _plot_options['Log x']
    except KeyError:
        if (p.get_coordinate_object('Frequency').mode.equidistant):
            _plot_options['Log x'] = False
        else:
            _plot_options['Log x'] = True
    try:
        _plot_options['Log y']
    except KeyError:
        _plot_options['Log y'] = True
            
    axes = ['Frequency']
    if (len(p.shape) == 1):
        if ((plot_type != 'xy') and (plot_type != 'scatter')):
            raise ValueError("For one signal only xy or scatter plot is possible.")
        p.plot(axes=axes,options=_plot_options)
    else:       
        if (plot_type == 'xy'):
            _plot_type = 'multi xy'
        else:
            _plot_type= plot_type
        if (plot_type == 'grid xy'):
            chn = np.array(p.coordinate('Channel',options={'Change only':True})[0]).size
            if (chn < 3):
                raise ValueError('Two small number of channels for grid xy plot. Use multi xy.')
            if (chn >= 25):
                column_number = 5
            elif (chn > 15):
                column_number = 3
            else:
                column_number = 2
            p = p.slice_data(slicing={'Channel':flap.Intervals(1,column_number,step=column_number)})
            axes=['Interval(Channel)','Interval(Channel) sample index','Frequency']
    plot_id1 = p.plot(plot_type=_plot_type,axes=axes,options=_plot_options,plot_id=plot_id)
    if (_plot_type == 'multi xy'):
        plt.legend(np.array(d.coordinate('Signal name',options={'Change only':True})[0]).flatten())
    return plot_id1
          
