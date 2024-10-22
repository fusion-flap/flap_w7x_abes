# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 15:03:53 2024

Functions to handle background corrected beam and background signals.
@author: Zoletnik
"""
import copy
import numpy as np
import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def get_clean_abes(exp_ID,signals='ABES-*',datapath=None,resample="",signal_type='raw',timerange=None):
    """
    Calculate background corrected beam signals or background signals. Can handle camera and timed chopping as well
    with all frequencies. 

    Parameters
    ----------
    exp_ID : string
        The experiment ID.
    datapath : string or None, optional
        The data path. If None the one in flap_defaults.cfg will be used. The default is None.
    signals : string or list of strings, optional
        The signals to process. The default is 'ABES-*'.
    resample : float of "", optional
        A float number indecates a resample frequency for reading the data. This us useful for flow chopper. 
        The default is "".
    signal_type : string, optional
        "raw": Now correction for background. 
        "beam": Background corrected beam signal.
        "background": The background signal.
        The default is 'raw'.
    timerange : List of two floats or None, optional
        Time range to process. The default is None.

    Returns
    -------
    d : flap.DataObject
        The processed data.
    """
    
    chopper_mode,beam_on_time,beam_off_time,period_time,d_beam_on,d_beam_off = flap_w7x_abes.chopper_parameters(exp_ID,datapath=datapath)
    options = {}
    if (period_time > 0.01):
        options['Resample'] = 1e4
    else:
        options['Resample'] = None
        if (resample != ""):
            print("Resample set for fast chopping. Signal processign might be incorrect.")
    if (resample != ""):
        options['Resample'] = resample
        
    d=flap.get_data('W7X_ABES',
                    exp_id = exp_ID,
                    name = signals,
                    options = options,
                    coordinates = {'Time': timerange}
                    )
    if (signal_type == 'beam'):
        d_on = d.slice_data(slicing={'Time':d_beam_on},
                            summing={'Rel. Sample in int(Time)':'Mean'},
                            options={'Regenerate':True}
                            )
        d_off = d.slice_data(slicing={'Time':d_beam_off},
                             summing={'Rel. Sample in int(Time)':'Mean'},
                             options={'Regenerate':True}
                             )
        d_off = d_off.slice_data(slicing={'Time':d_on},
                                 options={'Interpolation':'Linear'}
                                 )
        d_on.data = d_on.data - d_off.data
        d = d_on
    elif (signal_type == 'background'):
        d = d.slice_data(slicing={'Time':d_beam_off},
                         summing={'Rel. Sample in int(Time)':'Mean'},
                         options={'Regenerate':True}
                         )
    return d
      

def plot_clean_abes(exp_ID,signals='ABES-*',datapath=None,resample="",signal_type='raw',plot_type='xy',timerange=None,options={}):
    """
    Plot background corrected beam signals or background signals. Can handle camera and timed chopping as well
    with all frequencies. 
 
    Parameters
    ----------
    exp_ID : string
        The experiment ID.
    datapath : string or None, optional
        The data path. If None the one in flap_defaults.cfg will be used. The default is None.
    signals : string or list of strings, optional
        The signals to process. The default is 'ABES-*'.
    resample : float of "", optional
        A float number indecates a resample frequency for reading the data. This us useful for flow chopper. 
        The default is "".
    signal_type : string, optional
        "raw": Now correction for background. 
        "beam": Background corrected beam signal.
        "background": The background signal.
        The default is 'raw'.
    timerange : List of two floats or None, optional
        Time range to process. The default is None.
    plot_type: string, optional
        The plot type for the flap.plot method.
 
    Returns
    -------
    d : flap.DataObject
        The processed data.

    """
     
    d = get_clean_abes(exp_ID,
                       signals=signals,
                       datapath=datapath,
                       resample=resample,
                       signal_type=signal_type,
                       timerange=timerange
                       )
    axes = ['Time']
    _options = copy.deepcopy(options)
    if (len(d.shape) == 1):
        if ((plot_type == 'xy') and (plot_type != 'scatter')):
            raise ValueError("For one signal only xy or scatter plot is possible.")
        d.plot(axes=axes,options=_options)
    else:       
        if (plot_type == 'xy'):
            _plot_type = 'multi xy'
        else:
            _plot_type= plot_type
        if (plot_type == 'grid xy'):
            chn = np.array(d.coordinate('Channel',options={'Change only':True})[0]).size
            if (chn < 3):
                raise ValueError('Two small number of channels for grid xy plot. Use multi xy.')
            if (chn >= 25):
                column_number = 5
            elif (chn > 15):
                column_number = 3
            else:
                column_number = 2
            d = d.slice_data(slicing={'Channel':flap.Intervals(1,column_number,step=column_number)})
            axes=['Interval(Channel)','Interval(Channel) sample index','Time']
            try:
                _options['Y range']
            except KeyError:
                _options['Y range'] = [0,np.nanmax(d.data)]
                

        d.plot(plot_type=_plot_type,axes=axes,options=_options)
        if (_plot_type == 'multi xy'):
            plt.legend(np.array(d.coordinate('Signal name',options={'Change only':True})[0]).flatten())
            
        

 