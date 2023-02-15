#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:42:34 2023

@author: refydani
"""

from tdms_data_reading import read_tdms_data
import matplotlib.pyplot as plt

def plot_tdms_data(channels=['HV Em Meas Voltage','HV Em Meas Current','HV Ex Meas Voltage','HV Ex Meas Current','TC Top'],shot=None,group_name=None):
    plt.close('all')
    fig6, ax = plt.subplots(len(channels),sharex=True, sharey=False)
    fig6.suptitle("W7-X ABES beam parameters, shot#"+shot)
    for i in range(len(channels)):
        d=read_tdms_data(channels[i],shot=shot,group_name=group_name)
        ax[i].plot(d['time'],d['data'])
        ax[i].set_title(channels[i])
        ax[i].set_ylabel(d['unit'])
        if i == len(channels)-1:
            ax[i].set_xlabel("time")
    # plt.tight_layout()
    plt.show()

if __name__ == '__main__':  
    shot='20230214.029'
    group_name=None
    # channels=['HV Em Meas Voltage','HV Em Meas Current','HV Ex Meas Voltage','HV Ex Meas Current','VG HighVac1','FC2 Resistor Current mA','- FC Current mA']
    channels=['HV Em Meas Voltage','HV Em Meas Current','HV Ex Meas Voltage','HV Ex Meas Current','- FC Voltage V','FC1 Resistor Current mA','- FC Current mA']
    chdata=plot_tdms_data(channels=channels,shot=shot,group_name=group_name)