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
               plot_type='xy',options=None,plot_id=None,crosspower_amplitude=False,crosspower_phase=False,
               beam_on_start_delay=None,beam_on_end_delay=None,beam_off_start_delay=None,beam_off_end_delay=None):
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
    crosspower_ampltude: boolean
        If True calculates crosspower between consecutive signals and plots its amplitude. 
        For each signal in <signals> the crosspower will be calculated with the channel following 
        it towards the core plasma.
        Signal ABES-40 will be dropped.
    crosspower_phase: boolean
        Same is crosspower_ampltude, but plots the phase.
    beam_on_start_delay : float or None
        The start delay [microsec] to use for the beam on time relative to the one calculated from the settings.
        If None use the one determined by the data read program which may be non zero for measurements when 
        the camera was run on external timing.
    beam_on_end_delay : float or None
        The end delay [microsec] to use for the beam on time relative to the one calculated from the settings.
        If None use the one determined by the data read program which may be non zero for measurements when 
        the camera was run on external timing.
    beam_off_start_delay : float or None
        The start delay [microsec] to use for the beam off time relative to the one calculated from the settings.
        If None use the one determined by the data read program which may be non zero for measurements when 
        the camera was run on external timing.
    beam_off_end_delay : float or None
        The end delay [microsec] to use for the beam off time relative to the one calculated from the settings.
        If None use the one determined by the data read program which may be non zero for measurements when 
        the camera was run on external timing.
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
    if (crosspower_amplitude and crosspower_phase):
        raise ValueError("Only one of crosspower_amplitude and crosspower_phase can be set.")
    
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

    # Handling regular expressions in signal names
    _signals = copy.deepcopy(signals)
    if (type(_signals) is not list):
        _signals = [_signals]
    try:
        _signals, signal_index = flap.select_signals(d_beam_on.info['signal_list'],_signals)
    except ValueError as e:
        raise e

    # Here we have a signal list suitable for this ABES measurement
    
    if (not (crosspower_amplitude or crosspower_phase)):
        d = flap_w7x_abes.get_clean_abes(exp_ID,
                                         signals=_signals,
                                         datapath=datapath,
                                         resample=resample,
                                         signal_type=signal_type,
                                         timerange=timerange,
                                         beam_on_start_delay=beam_on_start_delay,
                                         beam_on_end_delay=beam_on_end_delay,
                                         beam_off_start_delay=beam_off_start_delay,
                                         beam_off_end_delay=beam_off_end_delay
                                         )
        if (beam_on is None):
            p = d.apsd(coordinate='Time',options=apsd_options)
        else:
            if (beam_on):
                d_chop = d_beam_on
            else:
                d_chop = d_beam_off
            p = d.apsd(coordinate='Time',intervals={'Time':d_chop},options=apsd_options)
    else:
        # Checking that only ABES-xx names are present
        while True:
            for i in range(len(_signals)):
                if (_signals[i][:5] != 'ABES-'):
                    del _signals[i]
                    break
            else:
                break
        if (len(_signals) == 0):
            raise ValueError("Corsspower can be calculated only for ABES-xx signals.")
                
        # Removing ABES-40 as there is no successive signal
        try:
            del _signals[_signals.index('ABES-40')]
        except ValueError:
            pass
        _signals_load = copy.deepcopy(_signals)
        for signal_act in _signals:
            signal_next = signal_act[:5] + str(int(signal_act[5:]) + 1)
            try:
                _signals_load.index(signal_next)
            except ValueError:
                _signals_load.append(signal_next)        
        d = flap_w7x_abes.get_clean_abes(exp_ID,
                                         signals=_signals_load,
                                         datapath=datapath,
                                         resample=resample,
                                         signal_type=signal_type,
                                         timerange=timerange,
                                         beam_on_start_delay=beam_on_start_delay,
                                         beam_on_end_delay=beam_on_end_delay,
                                         beam_off_start_delay=beam_off_start_delay,
                                         beam_off_end_delay=beam_off_end_delay
                                         )
        cpsd_names = []
        cpsd_channels = []
        for i_ch in range(len(_signals)):
            s1 = _signals[i_ch]
            s2 = s1[:5] + str(int(s1[5:]) + 1)
            cpsd_names.append(s1+'--'+s2)
            d1 = d.slice_data(slicing={'Signal name':s1})
            d2 = d.slice_data(slicing={'Signal name':s2})
            cpsd_channels.append(d1.coordinate('Channel',options={'Change only':True})[0][0])
            if (beam_on is None):
                p = d1.cpsd(ref=d2,coordinate='Time',options=apsd_options)
            else:
                if (beam_on):
                    d_chop = d_beam_on
                else:
                    d_chop = d_beam_off
                p = d1.cpsd(ref=d2,coordinate='Time',intervals={'Time':d_chop},options=apsd_options)
            if (i_ch == 0):
                frequency_coord = copy.deepcopy(p.coordinate('Frequency')[0])
                cpsd_data = np.ndarray((frequency_coord.size,len(_signals)),dtype=float)
                if ((p.error is not None) and crosspower_amplitude):
                    cpsd_error = np.ndarray((frequency_coord.size,len(_signals)),dtype=float)
                else:
                    cpsd_error = None
            if (crosspower_amplitude):
                cpsd_data[:,i_ch] = np.absolute(p.data)
                if (p.error is not None):
                    cpsd_error[:,i_ch] = p.error
            else:
                cpsd_data[:,i_ch] = np.angle(p.data)
 
        coord_list = copy.deepcopy(p.coordinates)
        # Removing coordinates and changing the dimension list of the Frequency
        while (True):
            for i in range(len(coord_list)):
                if (coord_list[i].unit.name == 'Frequency'):
                    coord_list[i].dimension_list = [0]
                    continue
                if (coord_list[i].unit.name == 'Signal name (Ref)'):
                    del coord_list[i]
                    break
                if (coord_list[i].unit.name == 'Channel (Ref)'):
                    del coord_list[i]
                    break
                if (coord_list[i].unit.name == 'Channel'):
                    del coord_list[i]
                    break
                if (coord_list[i].unit.name == 'Signal name'):
                    del coord_list[i]
                    break          
            else:
                break
        coord_signal = flap.Coordinate(name='Signal name',
                                       unit='n.a.',
                                       mode=flap.CoordinateMode(equidistant=False),
                                       values=cpsd_names,
                                       dimension_list=[1],
                                       shape=len(cpsd_names)
                                       )
        coord_channel = flap.Coordinate(name='Channel',
                                       unit='n.a.',
                                       mode=flap.CoordinateMode(equidistant=False),
                                       values=cpsd_channels,
                                       dimension_list=[1],
                                       shape=len(cpsd_channels)
                                       )
        coord_list.append(coord_signal)
        coord_list.append(coord_channel)
        if (len(_signals) == 1):
            cpsd_data = cpsd_data[:,0]
            if (cpsd_error is not None):
                cpsd_error = cpsd_error[:,0]
            coord_channel.dimension_list = []
            coord_signal.dimension_list = []
        p = flap.DataObject(cpsd_data,
                            error=cpsd_error,
                            data_unit=flap.Unit(name='Crosspower',unit='a.u'),
                            coordinates=coord_list,
                            exp_id=exp_ID,
                            data_title='Crosspower of consecutive channels',
                            data_source=p.data_source
                            )
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
        _plot_type = plot_type
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
          
