#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:45:43 2023

@author: refydani
"""

import os
import re
import matplotlib.pyplot as plt
import numpy as np
import glob

from nptdms import TdmsFile
from utc_offset import UTC_offset

def read_tdms_data(channel,shot=None,search_dir='/data/W7-X/APDCAM/'):
    print_avail=0
    search_dir = '/data/W7-X/APDCAM/'
    dl=glob.glob(search_dir+'*'+shot+'*')
    if len(dl) > 0:
        for search_dir in dl:
            fl=glob.glob(search_dir+'/*'+shot+'*.tdms')
            if len(fl) == 0: 
                raise ValueError("Shot directory not found in '{:s}'".format(search_dir))
    else:
        raise ValueError("Shot directory not found in '{:s}'".format(search_dir))
    data=[]
    time=[]
    file_list=[]
    group_list=[]
    channel_list=[]
    for fn in fl:
        fn = os.path.join(search_dir,fn)
        file_list.append(fn)
        with TdmsFile.open(fn) as tdms_file:
            # print(tdms_file.groups())
            for group in tdms_file.groups():
                group_list.append(group.name)
                # print(tdms_file[group.name].channels())
                for ch in tdms_file[group.name].channels():
                    channel_list.append(ch.name)
                    # print(ch.name)
                    if re.match(channel,ch.name):
                        tdms_data = tdms_file[group.name][channel]
                        unit=tdms_data.properties['Unit']
                        # print(tdms_data[:])
                        # break
                    else:
                        if re.match('TimeStamp',ch.name):
                            tdms_time = tdms_file[group.name]['TimeStamp']
                        else:
                            continue
                try: 
                    data.append(tdms_data[:])
                    time.append(tdms_time[:])
                except NameError:
                    print("Channel '{:s}' is not found in file {:s}.".format(str(channel),fn))
                    # return
                # data.append(tdms_data[:])
    time=np.array(np.concatenate(time).flat)
    data=np.array(np.concatenate(data).flat)
    arr1inds = time.argsort()
    time = time[arr1inds[::-1]]
    data = data[arr1inds[::-1]]
    time = time[~np.isnan(data)]
    data = data[~np.isnan(data)]
    if print_avail:
        file_list=[*set(file_list)]
        group_list=[*set(group_list)]
        channel_list=[*set(channel_list)]
        print('files available:')
        print(file_list)
        print('groups available:')
        print(group_list)
        print('channels available:')
        print(channel_list)
    return time,data,unit
    
if __name__ == '__main__':  
    shot='T20230207.002'
    time, data,unit=read_tdms_data('HV Em Meas Voltage',shot=shot)
    plt.plot(time,data)
    plt.ylabel(unit)
    plt.xlabel('time')
    plt.show()