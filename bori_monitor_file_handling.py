# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 21:15:27 2021

@author: Zoletnik
"""

import matplotlib.pyplot as plt
import numpy as np
try:
    try:
        from .nptdms_mod import TdmsFile
    except ImportError:
        from nptdms_mod import TdmsFile
except ModuleNotFoundError:
        print("nptdms module not available, beam monitor file handling cannot be done.")

import datetime
import os
import copy

try:
    from .utc_offset import UTC_offset
except ImportError:
    from utc_offset import UTC_offset


def find_files(startdate=None,starttime='0000',start_datetime=None,
               enddate=None,endtime='2359',end_datetime=None,
               datapath='',UTC_offset_minutes=None,verbose=True
               ):
    """
    Finds the files in which data for the indicated time interval (local time) can be found.

    Parameters
    ----------
    startdate : string, optional
        The start date YYYYMMDD. If not set the actual daye is assumed.
    starttime : string, optional
        The start time in local time, format: HHMM. The default is '0000'. The time is local time
    start_datetime: numpy datetime64 object
        Alternative to specify the start time. If this is set startdate and stattime is not used.
    enddate: string, optional
        The end date YYYYMMDD. If not set the same day is assumed as the start date.
    endtime: string, optional
        The end time in local time, format: HHMM. The default is '2459'. The time is local time.
    end_datetime: numpy datetime64 object
        Alternative to specify the end time. If this is set enddate and endtime is not used.
    datapath : string
        The data path.
    UTC_offset_minutes : int, optional
        The time offset of the local time realtive to UTC. 
        If None it will be determined from the date.
    verbose: boolean
        Print messages about progress of processing.

    Raises
    ------
    ValueError
        No suitable file found.

    Returns
    -------
    Returns the filenames, start and end times sorted by increasing starttime.
    file_list: list of strings
        The list of suitable filenames (with full path)
    file_starttime_list: list of np.datetime64 objects
        The start times of the files in local times.
    file_endtime_list: list of np.datetime64 objects
        The end times of the files in local times.
        

    """
    
    if (start_datetime is None):
        if (startdate is None):
            _startdate = str(datetime.date.today().year) + str(datetime.date.today().month) + str(datetime.date.today().day)
        else:
            _startdate = startdate
        _start_datetime = np.datetime64(_startdate[:4]+'-'+_startdate[4:6]+'-'+_startdate[6:8]+'T'+
                                        starttime[:2]+':'+starttime[2:4]
                                        )
    else:
        _start_datetime = start_datetime

    if (end_datetime is None):        
        if (enddate is None):
            _enddate = _startdate
        else:
            _enddate = enddate
        _end_datetime = np.datetime64(_enddate[:4]+'-'+_enddate[4:6]+'-'+_enddate[6:8]+'T'+
                                      endtime[:2]+':'+endtime[2:4]
                                      )
    else:
        _end_datetime = end_datetime
        
    date_str =str(_start_datetime)
    _date = date_str[:4] + date_str[5:7] + date_str[8:10]
    _UTC_offset_minutes = UTC_offset(UTC_offset_minutes=UTC_offset_minutes,date=_date)   
    file_starttime_list = []
    file_endtime_list = []
    fname_list = []
    # Starting search 13 h before the desired time as one log file is about 11 hour long
    act_datetime = _start_datetime - np.timedelta64(13,'h')
    dirname = datapath
    try:
        l = os.listdir(dirname)
    except FileNotFoundError:
        raise ValueError("No logfiles found, bad datapath?")
            
    for f in l:
        ind = f.find('.tdms')
        if (ind < 0):
            continue
        if (f[ind:] != '.tdms'):
            continue
        if (f[:5] == 'bori_'):
            try:
               file_start_datetime = np.datetime64(f[5:9]+'-'+f[9:11]+'-'+f[11:13]+'T'+f[14:16]+':'+f[16:18]+':'+f[18:20]) 
            except ValueError:
                continue
        else:
            try:
               file_start_datetime = np.datetime64(f[0:4]+'-'+f[4:6]+'-'+f[6:8]+'T'+f[9:11]+':'+f[11:13]+':'+f[13:15]) 
            except ValueError:
                continue            
        if (file_start_datetime > _end_datetime):
            continue
        if (file_start_datetime < _end_datetime - np.timedelta64(2,'D')):
            continue
        fname = os.path.join(dirname,f)
        if (verbose):
            print('Checking {:s}'.format(fname),flush=True)
        try:
            with TdmsFile.open(fname) as tdms_file:                    
#                tdms_version =  tdms_file.properties['Version']
                try:
                    ch_t = tdms_file['MonitorData']['TimeStamp']
                    file_start_datetime_from_file = ch_t[0] + np.timedelta64(_UTC_offset_minutes, 'm')
                    file_end_datetime_from_file = ch_t[-1] + np.timedelta64(_UTC_offset_minutes, 'm')
                    MAXDIFF = 500 # maximum time difference between time in file name and contents [s]
                    for i in range(11):
                        if (abs(file_start_datetime_from_file - file_start_datetime) > np.timedelta64(MAXDIFF,'s')):  
                            file_start_datetime_from_file = ch_t[i + 1] + np.timedelta64(_UTC_offset_minutes, 'm')
                            print("Warning: timestamp [{:d}] and time from file name are significantly different. (file: {:s})".format(i,fname),flush=True)
                        else:
                            break
                    if (abs(file_start_datetime_from_file - file_start_datetime) > np.timedelta64(MAXDIFF,'s')):  
                        print("Warning: Cannot find reasonable start time in file (file: {:s})".format(fname),flush=True)
                        print("Omitting this file.",flush=True)
                        continue                 
                    if (file_start_datetime_from_file > _end_datetime):
                        continue
                    for i in range(10):
                        if (file_end_datetime_from_file < file_start_datetime_from_file):
                            file_end_datetime_from_file = ch_t[-(2 + i)] + np.timedelta64(_UTC_offset_minutes, 'm')
                            print("Warning: timestamp [{:d}] in file is bad. (file: {:s})".format(-(i + 1),fname),flush=True)
                        else:
                            break
                except (IndexError, IOError, ValueError):
                    print("Broken file? {:s}".format(fname))
                    continue
        except ValueError:
            print("Broken file: {:s}".format(fname))
            continue
        if (file_end_datetime_from_file < file_start_datetime_from_file):  
            print("Warning: Cannot find reasonable end time in file (file: {:s})".format(fname),flush=True)
            print("Omitting this file.",flush=True)
            continue  
        if (file_end_datetime_from_file < _start_datetime):
            continue
        fname_list.append(fname)
        file_starttime_list.append(file_start_datetime)
        file_endtime_list.append(file_end_datetime_from_file)
    if (len(fname_list) == 0):
        raise ValueError('No data found.')
    rel_starttime = (np.array(file_starttime_list) - np.amin(np.array(file_starttime_list))).astype(float)
    ind = np.argsort(rel_starttime)
    fname_sorted = []
    starttime_sorted = []
    endtime_sorted = []
    for i in ind.tolist():
        fname_sorted.append(fname_list[i])
        starttime_sorted.append(file_starttime_list[i])
        endtime_sorted.append(file_endtime_list[i])
    return fname_sorted,starttime_sorted,endtime_sorted
            
def read_data(data_names=None,startdate=None,starttime='0000',start_datetime=None,
              enddate=None,endtime='2359',end_datetime=None,
               datapath='',UTC_offset_minutes=None,verbose=True):
    """
    Reads multiple data from the TDMS files specified by the start and end dates/times.

    Parameters
    ----------
    data_names : string or list of strings
        The channel names to read. The default is None.
    startdate : string, optional
        The start date YYYYMMDD. If not set the actual daye is assumed.
    starttime : string, optional
        The start time HHMM. The default is '0000'. The time is local time
    start_datetime: numpy datetime64 object
        Alternative to specify the start time. If this is set startdate and stattime is not used.
    enddate: string, optional
        The end date YYYYMMDD. If not set the same day is assumed as the start date.
    endtime: string, optional
        The end time HHMM. The default is '2459'. The time is local time.
    end_datetime: numpy datetime64 object
        Alternative to specify the end time. If this is set enddate and endtime is not used.
    datapath : string, optional
        The data path. The default is for the beam computer.
    UTC_offset_minutes : int, optional
        The time offset of the local time realtive to UTC. 
        If None it will be determined from the date.
    verbose: boolean
        Print messages about progress of processing.

    Raises
    ------
    ValueError
        No suitable file found.
    IOError
        Error deadeing data.

    Returns
    -------
    time : numpy datetime64 array
        The timstamps for the measurements. 
        Use (time-time[0])/np.timedelta64(1,'s') to convert to time relative to first time in seconds.
    data : List of numpy float arrays
        The data.
    data_unit : List of string(s)
        List of unit(s) corresponding to data

    """
    if (starttime is None):
        _starttime = '0000'
    else:
        _starttime= starttime
    if (endtime is None):
        _endtime = '2359'
    else:
        _endtime= endtime
    if (len(_starttime) != 4):
        raise ValueError('starttime should be given as HHMM')
    if (len(_endtime) != 4):
        raise ValueError('endtime should be given as HHMM')
    if (data_names is None):
        raise ValueError('Data names should be set.')
    if (type(data_names) is not list):
        _data_names = [data_names]
    else:
        _data_names = data_names
    if (start_datetime is None):
        if (startdate is None):
            _startdate = str(datetime.date.today().year) + str(datetime.date.today().month) + str(datetime.date.today().day)
        else:
            if (len(startdate) != 8):
                raise ValueError("date should be given as YYYMMDD")
            _startdate = startdate
        _start_datetime = np.datetime64(_startdate[:4]+'-'+_startdate[4:6]+'-'+_startdate[6:8]+'T'+
                                        _starttime[:2]+':'+_starttime[2:4]
                                        )
    else:
        _start_datetime = start_datetime

    if (end_datetime is None):        
        if (enddate is None):
            _enddate = _startdate
        else:
            if (len(enddate) != 8):
                raise ValueError("date should be given as YYYMMDD")
            _enddate = enddate
        _end_datetime = np.datetime64(_enddate[:4]+'-'+_enddate[4:6]+'-'+_enddate[6:8]+'T'+
                                      _endtime[:2]+':'+_endtime[2:4]
                                      )
    else:
        _end_datetime = end_datetime

    date_str =str(_start_datetime)
    _date = date_str[:4] + date_str[5:7] + date_str[8:10]
    _UTC_offset_minutes = UTC_offset(UTC_offset_minutes=UTC_offset_minutes,date=_date)   

    fnames, starttimes, endtimes = find_files(start_datetime=_start_datetime,end_datetime=_end_datetime,
                                              datapath=datapath,UTC_offset_minutes=_UTC_offset_minutes,verbose=verbose
                                              )

    time = np.ndarray(0,dtype=np.datetime64)
    data = [np.ndarray(0,dtype=float)]*len(_data_names)
    for fn in fnames:
        if (verbose):
            print('Processing {:s}'.format(fn),flush=True)
        with TdmsFile.open(fn) as tdms_file:
            try:
                tdms_version =  tdms_file.properties['Version']
            except KeyError:
                tdms_version = 1

            ch_t = tdms_file['MonitorData']['TimeStamp']
            try:
                t = ch_t[:] + np.timedelta64(_UTC_offset_minutes, 'm')
                ind = np.nonzero(np.logical_and(t >= _start_datetime,
                                                t <= _end_datetime 
                                                )
                                 )[0]
            except (IndexError, IOError, ValueError):
                print("Error reading from file: {:s}".format(fn))
                continue      
            if (len(ind) == 0):
                print("Cannot find data in time interval in file: {:s}".format(fn))
                continue
            data_unit=[]
            read_error = False
            data_save = copy.deepcopy(data)
            for i,signame in enumerate(_data_names):
                if (verbose):
                    print("Reading '{:s}'".format(signame),flush=True)
                ch = tdms_file['MonitorData'][signame] 
                d = ch[:]
                try:
                    data[i] = np.concatenate((data[i],d[ind].astype(float)))
                except IndexError:
                    read_error = True
                    break
            if (read_error):
                data = data_save
                if (verbose):
                    print("Broken file '{:s}'".format(fn),flush=True)           
            else:
                time = np.concatenate((time,t[ind]))
    if (len(time) < 2):
        raise ValueError("NO data found.")
    dt = np.diff((time - time[0]) / np.timedelta64(1,'s'))
    ind_bad = np.nonzero(dt <= 0)[0]
    if (len(ind_bad) != 0):
        print('Removed {:d} bad time points.'.format(len(ind_bad)),flush=True)
        ind_good = np.nonzero(dt > 0)[0]
        if (len(ind_good) == 0):
            raise ValueError("No good time points found.")
        time = np.concatenate((np.array([time[0]]),time[ind_good + 1]))
        for i in range(len(data)):
            data[i] = np.concatenate((np.array([data[i][0]]),data[i][ind_good + 1]))
    
    return time,data,data_unit

def read_channels(startdate,datapath,page='MnitorData'):
    channels = []
    file = find_files(startdate,datapath=datapath)
    with TdmsFile.open(file[0][0]) as tdms_file:
            ch_t = tdms_file[page]
    channel_list=ch_t.channels()
    for c in channel_list:
        channels.append(c.name)
    return channels

def channel_list(file,page='MonitorData'):
    """
    Return the list of channels in the TDMS file.

    Parameters
    ----------
    file : string
        The full file name.

    Returns
    -------
    channels : list of strings
        The channel list.

    """
    channels = []
    with TdmsFile.open(file) as tdms_file:
        ch_t = tdms_file[page]
    channel_list=ch_t.channels()
    for c in channel_list:
        channels.append(c.name)
    return channels
    
def page_list(file):
   with TdmsFile.open(file) as tdms_file:
       groups = tdms_file.groups()
       group_names = []
       for g in groups:
           group_names.append(g.name)
   return group_names
    
    
# def read_exp_tdms(exp_id,datapath='/data/W7X/APDCAM'):
#     dirname = os.path.join(datapath,exp_id)
#     files = os.listdir(dirname)
#     filelist = []
#     phase_list = ['Neut-Active','HV-Raise','HV-Beam']
#     phase_data = [[],[],[]]
#     for f in files:
#         if f[-5:] == '.tdms':
#             d = {'Filename':os.path.join(dirname,f)}
#             d['Pages'] = page_list(os.path.join(dirname,f))
#             for p in d['Pages']:
#                 d[p] = {'Channels':channel_list(os.path.join(dirname,f),page=p)}
#                 with TdmsFile.open(os.path.join(dirname,f)) as tdms_file:
#                     d[p]['Start time'] = tdms_file[p]['TimeStamp'][0] 
#                     d[p]['End time'] = tdms_file[p]['TimeStamp'][-1]             
#                 for i,ph in enumerate(phase_list):
#                     if (p[:len(ph)] == ph):
#                         pdi = {'Filename' : os.path.join(dirname,f),
#                                'Start time' : d[p]['Start time'],
#                                'End time' : d[p]['End time'] 
#                                } 
#                         for ii in range()
#                         phase_data[i].append(
#                         phase_data[i]
                        
                        
#              filelist.append(d)
#     print(phase_list)
    
            
