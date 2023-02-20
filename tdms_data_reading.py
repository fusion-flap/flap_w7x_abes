#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 16:38:29 2023

@author: refydani
"""

import os
import re
import numpy as np
import glob
import pickle
import tqdm
import matplotlib.pyplot as plt
from nptdms import TdmsFile


def loadfile(filename):
    with open(filename,'rb') as f:  
        return pickle.load(f)    
    
def savefile(filename,variables):
    with open(filename, 'wb') as f:  
        pickle.dump(variables, f)  
def read_one_channel(channel,shot=None,search_dir='/data/W7-X/APDCAM/'):
    dl=glob.glob(search_dir+'*'+shot+'*')
    if len(dl) > 0:
        for search_dir in dl:
            fl=glob.glob(search_dir+'/*'+shot+'*.tdms')
            if len(fl) == 0: 
                raise ValueError("TDMS file not found in '{:s}'".format(search_dir))
    else:
        raise ValueError("Shot ("+shot+") directory not found in '{:s}'".format(search_dir))
    tdms_data=[]
    for fn in fl:
        fn = os.path.join(search_dir,fn)
        with TdmsFile.open(fn) as tdms_file:
            for group in tdms_file.groups():
                for ch in tdms_file[group.name].channels():
                    # if re.match(channel,ch.name):
                    if channel == ch.name:
                        # print(group.name,ch.name)
                        chdata=tdms_file[group.name][ch.name]
                        data =chdata[:]
                        try:
                            unit=chdata.properties['Unit']
                        except:
                            unit=''
                        # unit=chdata.properties['Unit']
                        tdms_data.append({'shot':shot,'group':group.name,'channel':ch.name,'data':data,'unit':unit})
    if not tdms_data:
        raise ValueError("Channel '{:s}' is not found in file {:s}.".format(str(channel),fn))
    return tdms_data

def read_tdms_data(channel,shot=None,search_dir='/data/W7-X/APDCAM/',save_dir='/data/W7-X/Beam/tdms_processed/',group_name=None):
    # search_dir = '/data/W7-X/APDCAM/'
    file_save =save_dir+shot+'_tdms_processed.pkl'
    if group_name == None:
        group_name='.'
    else:
        group_name='.*'+group_name+'.*'
    tdms_data=[]
    try:
        tdms_data=loadfile(file_save)
    except:
        tdms_data=read_one_channel(channel,shot=shot,search_dir=search_dir)
        savefile(file_save,tdms_data)   
    if len(list(filter(lambda d: d['channel'] == channel, tdms_data))) == 0:
        tdms_data.extend(read_one_channel(channel,shot=shot,search_dir=search_dir))
        savefile(file_save,tdms_data) 
    if len(list(filter(lambda d: d['channel'] == 'TimeStamp', tdms_data))) == 0:
        tdms_data.extend(read_one_channel('TimeStamp',shot=shot,search_dir=search_dir))
        savefile(file_save,tdms_data) 
    data=[]
    time=[]
    for d in list(filter(lambda d: d['channel'] == channel, tdms_data)):
        if re.match(group_name,d['group']):
            data.append(d['data'])
            unit=d['unit']
    for d in list(filter(lambda d: d['channel'] == 'TimeStamp', tdms_data)):
        if re.match(group_name,d['group']):
            time.append(d['data'])
    time=np.array(np.concatenate(time).flat)
    data=np.array(np.concatenate(data).flat)
    arr1inds = time.argsort()
    time = time[arr1inds]
    data = data[arr1inds]
    time = time[~np.isnan(data)]
    data = data[~np.isnan(data)]
    # return time,data,unit
    return {'channel':channel,'data':data,'time':time,'unit':unit}
    
def read_tdms_data_slow(channel,shot=None,search_dir='/data/W7-X/APDCAM/',group_name=None):
    print_avail=1
    search_dir = '/data/W7-X/APDCAM/'
    if group_name == None:
        group_name='.'
    else:
        group_name='.*'+group_name+'.*'
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
            for group in tdms_file.groups():
                group_list.append(group.name)
                if re.match(group_name,group.name):
                    for ch in tdms_file[group.name].channels():
                        channel_list.append(ch.name)
                        if re.match(channel,ch.name):
                            tdms_data = tdms_file[group.name][channel]
                            unit=tdms_data.properties['Unit']
                            if ~print_avail: break
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
        for channel in channel_list:
            print(channel)
    return {'channel':channel,'data':data,'time':time,'unit':unit}

if __name__ == '__main__':  
    shot='20230216.069'
    search_dir=r'C:/Users/refyd/Documents/BES/W7X/data/'
    save_dir=r'C:/Users/refyd/Documents/BES/W7X/tdms_processed/'
    channel='HV Em Meas Current'
    # d=read_tdms_data(channel,shot=shot,search_dir=search_dir,save_dir=save_dir)
    d=read_tdms_data(channel,shot=shot)
  
    plt.plot(d['time'],d['data'])
    plt.ylabel(d['unit'])
    plt.xlabel('time')
    plt.show()