# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 14:40:26 2022

@author: Zoletnik
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

try:
    from .bori_monitor_file_handling import *
except ImportError:
    from flap_w7x_abes.bori_monitor_file_handling import *
    
def plot_beamdata(startdate=None,starttime=None,endtime=None,enddate=None,datapath='',start_datetime=None,end_datetime=None,figure=None,
                  R_emit=81,R_ext=73, last_minutes=None):
    """
    Plot beam data.

    Parameters
    ----------
    startdate : string
        The start date of processing, YYYYMMDD
    starttime : string, optional
        The start time of processing, HHMM The default is None, whiche means 0000
    endtime : string, optional
        The end time of processing, HHMM. The default is None, which means 2400
    enddate : string, optional
        The end date of processing, YYYYMMDD. The default is None which means use the same date as startdate.
    datapath : string, optional
        The data path. The default is ''.
    start_datetime: numpy datetime64 object
        Alternative to specify the start time. If this is set startdate and stattime is not used.
   end_datetime: numpy datetime64 object
       Alternative to specify the end time. If this is set enddate and endtime is not used.
    figure : Matplotlib figure, optional
        The figure to use. The default is None, meaning make new figure.
    R_emit : float, optional
        The resistor on the emitter power supply [MOhm]. The default is 81.
        If None, will calculate from currents.
    R_ext : float, optional
        The resistor on the emitter power supply [MOhm]. The default is 73.
        If None, will calculate from currents.
    last_minutes : float
        Plot only the last minutes indicated by this argument.

    Returns
    -------
    None.

    """
    
    labelsize = 10
    linewidth = 1
    plt.rcParams['lines.linewidth'] = linewidth
    plt.rcParams['axes.linewidth'] = linewidth
    plt.rcParams['axes.labelsize'] = labelsize 
    plt.rcParams['axes.titlesize'] = labelsize 
    plt.rcParams['xtick.labelsize'] =  labelsize 
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = linewidth
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.minor.width'] = linewidth
    plt.rcParams['xtick.minor.size'] = 2
    plt.rcParams['ytick.labelsize'] = labelsize 
    plt.rcParams['ytick.major.width'] = linewidth
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.minor.width'] = linewidth
    plt.rcParams['ytick.minor.size'] = 2
    plt.rcParams['legend.fontsize'] = labelsize 
#    plt.rcParams['suptitle.fontsize'] = 10
    
    data_names = ['Emit Current A','HV Em Meas Voltage','HV Ex Meas Voltage','HV Em Meas Current','HV Ex Meas Current',
                  'TC Oven Top','TC Oven Bottom','TC Torus Side Cone','TC Emit Side Cone','FC1 in','FC2 in','FC Polarity',
                  'FC1 Resistor Current mA','FC2 Resistor Current mA','VG HighVac1','VG HighVac2',
                  'Neut Shut Closed']
    t,d,u = read_data(data_names=data_names,startdate=startdate,starttime=starttime,endtime=endtime,datapath=datapath)
    d_dict = {}
    for key,data in zip(data_names,d):
        d_dict[key] = data 
    
    if ((t[-1] -t[0]) / np.timedelta64(1,'s') > 300):
        time = (t - t[0]) / np.timedelta64(1,'m')
        time_unit = 'min'
    else:
        time = (t - t[0]) / np.timedelta64(1,'s')
        time_unit = 's'
        
    if (figure is None):
        plt.close('all')
        figure= plt.figure(figsize=(25,17))
    else:
        plt.figure(figure)
    gs = gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=plt.gca(),hspace=0.5,wspace=0.5)
    plt.suptitle('Start time: {:s}, R_emit={:3.0f}MOhm  R_ext={:3.0f}MOhm'.format(str(t[0]),R_emit,R_ext),fontsize=16)


    if (last_minutes is not None):
        if (time_unit == 'min'):
            xlim = [time[np.nonzero(time > time[-1] - last_minutes)[0][0]], time[-1] ]
        else:
            print('time')
            xlim = [time[np.nonzero(time > time[-1] - last_minutes)[0][0] * 60], time[-1]]
    else:
        xlim = None
    
    ax = plt.subplot(gs[0,0:2])
    plt.plot(time,d_dict['Emit Current A'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[A]')
    plt.title('Emitter heating current')
    if (xlim is not None):
        plt.xlim(*xlim)
    
    plt.subplot(gs[0,2:4],sharex=ax)
    plt.plot(time,d_dict['HV Em Meas Voltage'])
    plt.plot(time,d_dict['HV Ex Meas Voltage'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[kV]')
    plt.title('Accelerator voltages')
    plt.legend(['$HV_{emit}$','$HV_{ext}$'])
    
    plt.subplot(gs[1,2:4],sharex=ax)
    plt.plot(time,d[3])
    plt.plot(time,d[4])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[mA]')
    plt.title('HV PS currents')
    plt.legend(['$HV_{emit}$','$HV_{ext}$'])
    trange = plt.xlim()
    
    nonzero_current = np.logical_and(d[3] > 0.05,d[4] > 0.05)
    ind = np.nonzero(np.logical_and(d_dict['HV Em Meas Voltage'] < d_dict['HV Ex Meas Voltage'],nonzero_current))[0]
    if (len(ind) != 0):
        if (R_emit is None):
            R_emit = np.mean(d_dict['HV Em Meas Voltage'][ind] / d_dict['HV Em Meas Current'][ind])
        if (R_ext is None):
            R_ext = np.mean(d_dict['HV Ex Meas Voltage'][ind] / d_dict['HV Ex Meas Current'][ind])
        print("R_emit={:3.0f}MOhm  R_ext={:3.0f}MOhm".format(R_emit,R_ext))
        I_emit_ext = d_dict['HV Ex Meas Voltage'] / R_ext -  d[4] 
        I_beam = d_dict['HV Em Meas Current'] - d_dict['HV Ex Meas Voltage'] / R_emit - I_emit_ext
        plt.subplot(gs[2,2:4],sharex=ax)
        plt.plot(time,I_beam)
        plt.plot(time,I_emit_ext)
        plt.xlabel('Time [{:s}]'.format(time_unit))
        plt.ylabel('[mA]')
        plt.title('Beam current')
        plt.legend(['Beam current','E-ext current'])
        
    plt.subplot(gs[1,0:2],sharex=ax)
    legend = []
    ind = np.nonzero(np.logical_and(d_dict['TC Oven Bottom'] > 20, d_dict['TC Oven Bottom'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Oven Bottom'],20,300))
        legend.append('Oven')
    ind = np.nonzero(np.logical_and(d_dict['TC Oven Top'] > 20, d_dict['TC Oven Top'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Oven Top'],20,300))
        legend.append('Top')
    ind = np.nonzero(np.logical_and(d_dict['TC Torus Side Cone'] > 20, d_dict['TC Torus Side Cone'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Torus Side Cone'],20,300))
        legend.append('Cone torus side')
    ind = np.nonzero(np.logical_and(d_dict['TC Emit Side Cone'] > 20, d_dict['TC Emit Side Cone'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Emit Side Cone'],20,300))
        legend.append('Cone emit side')
    plt.legend(legend)    
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[C]')
    plt.title('Neutralizer temperatures')
    plt.xlim(*trange)

    ind_shut = np.nonzero(d_dict['Neut Shut Closed'] != 0)[0]
    fc2_act = d_dict['FC2 in']
    fc2_act[ind_shut] = 0                         
    ind1 = np.nonzero(d_dict['FC1 in'] == 1)[0]
    ind2 = np.nonzero(fc2_act == 1)[0]   
    if ((len(ind1) != 0) or (len(ind2 != 0))):
        legend = []
        plt.subplot(gs[3,2:4],sharex=ax)
        if (len(ind1) != 0):
            plt.plot(time[ind1],d_dict['FC1 Resistor Current mA'][ind1])
            legend.append('FC1 Resistor Current mA')
        if (len(ind2) != 0):
            plt.plot(time[ind2],d_dict['FC2 Resistor Current mA'][ind2])
            legend.append('FC2 Resistor Current mA')
        plt.legend(legend)    
        plt.xlabel('Time [{:s}]'.format(time_unit))
        plt.ylabel('[mA]')
        plt.title('Faraday cup currents')
        plt.xlim(*trange)
        

    plt.subplot(gs[2,0:2],sharex=ax)
    plt.plot(time,d_dict['VG HighVac1'])
    plt.plot(time,d_dict['VG HighVac2'])
    plt.legend(['High vac 1','High vac 2'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[mbar]')
    plt.title('Vacuum pressure')
    plt.xlim(*trange)
    plt.yscale('log')
    
    plt.subplot(gs[3,0:2],sharex=ax)
    I_extra_emit = d_dict['HV Em Meas Voltage'] / R_emit - d_dict['HV Em Meas Current']
    I_extra_ext = d_dict['HV Ex Meas Voltage'] / R_ext - d_dict['HV Ex Meas Current']
    plt.plot(time,I_extra_emit)
    plt.plot(time,-I_extra_ext)
    plt.legend(['Into Emit. PS','From Ext. PS'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[mA]')
    plt.title('Extra currents')
    plt.xlim(*trange)

if __name__ == "__main__":
    datetoget=sys.argv[1]
    plot_beamdata(startdate=datetoget,datapath='/data',last_minutes=20)
    plt.savefig("/home/apdcam/Measurement/beam_20min.png")
    #plot_beamdata(startdate='20240918',datapath='z:/Data/W7-X_ABES/Beam',figure=fig)
    # fig = plt.figure(figsize=(25,17))
    # plot_beamdata(startdate='20230223',datapath='c:/Users/Zoletnik/Root/tmp',figure=fig)    
