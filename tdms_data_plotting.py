#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:42:34 2023

@author: refydani
"""

from tdms_data_reading import read_tdms_data
import matplotlib.pyplot as plt

def plot_tdms_data(channels=None,shot=None):
    plt.close('all')
    fig6, ax = plt.subplots(len(channels),sharex=True, sharey=False)
    fig6.suptitle("W7-X ABES beam parameters, shot#"+shot)
    for i in range(len(channels)):
        time,data,unit=read_tdms_data(channels[i],shot=shot)
        ax[i].plot(time,data)
        ax[i].set_title(channels[i])
        ax[i].set_ylabel(unit)
        if i == len(channels):
            ax[i].set_xlabel("time")
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':  
    shot='T20230210.003'
    channels=['HV Em Meas Voltage','HV Em Meas Current','HV Ex Meas Voltage','HV Ex Meas Current','TC Oven']
    chdata=plot_tdms_data(channels=channels,shot=shot)