#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:35:45 2022

@author: apdcam
"""
import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def show_all_abes(exp_id,time=[0,5],yrange=[-0.1,0.2],):
    """
    Plot the 40 ABES signals on the beam in a 5x8 pot matrix.

    Parameters
    ----------
    exp_id : string
        The experiment ID.
    time : list of two numbers, optional
        The time range in seconds. If None, will plot all data of the shot. 
        The default is [0,5].

    Returns
    -------
    None.

    """

    d=flap.get_data('W7X_ABES',
                    exp_id=exp_id,
                    name='ABES*',
                    coordinates={'Time':time},
                    options={'Resample':1e3}
                    )
    d = d.slice_data(slicing={'Channel':flap.Intervals(1,5,step=5)})

    d.plot(plot_type='grid xy',axes=['Interval(Channel)','Interval(Channel) sample index','Time'],
            options={'Y range':yrange})    
    plt.suptitle(exp_id)

def show_all_abes_spectra(exp_id,time=[0,1],fres=20,frange=[100,1E6],log_fres=True,beam_on=True):
    """
    Plot the power spectra of all 40 ABES channels.

    Parameters
    ----------
    exp_id : string
        The experiment ID.
    time : list of two numbers, optional
        The time range in seconds. If None, will plot all data of the shot. 
        The default is [0,5].
    fres : float, optional
        Frequency resolution in Hz. The default is 20. If 
    frange : list of two numbers, optional
        The frequency range in Hz. The default is [100,1e6]
    log_fres : bolean, optioinal
        If True the frequency resolution will change proportionally with frequency.
    beam_on : boolean, optional
        If True will process beam-on inervals, otherwise beam-off (chopped)     
        
    Returns
    -------
    None.

    """

    # Loading data
    d=flap.get_data('W7X_ABES',
                    exp_id=exp_id,
                    name='ABES*',
                    coordinates={'Time':time}
                    )
    # Finding processing intervals
    if (beam_on):
        d_chop=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0},'Start':1000,'End':-1000}
                             )
    else:
        d_chop=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0},'Start':1000,'End':-1000}
                             )
    # Calculating power spectra in processing intervals
    d_ps = d.apsd(coordinate='Time',
                  intervals={'Time':d_chop},
                  options={'Interval':1, 'Range':frange,'Res':fres,'Log': log_fres}
                 )
    # Arranging signals in a 5x8 matrix in preparation for plotting
    d_ps = d_ps.slice_data(slicing={'Channel':flap.Intervals(1,5,step=5)})

    d_ps.plot(plot_type='grid xy',
              axes=['Interval(Channel)','Interval(Channel) sample index','Frequency'],
              options={'Log y': True,'Log x': True}
              )    
    plt.suptitle(exp_id)

 