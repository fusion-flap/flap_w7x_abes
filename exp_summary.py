# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 21:16:22 2024

@author: Zoletnik
"""
import numpy as np
import os
import re
import time
import matplotlib.pyplot as plt
import pickle

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def exp_summary(exp_ID,timerange=None,datapath=None,channels=range(10,26),test=False):
    """
    Return a single line summary of the experiment: expID, max APDCAM signal (uncalibrated), time range, channel range

    Parameters
    ----------
    exp_ID : string
        The exp ID, YYYMMDD.xxx
    timerange : list, optional
        The processing timerange. The default is None.
    datapath : string, optional
        The datapath. The default is None.

    Returns
    -------
    txt : string
        The experiment summary.
    data : dict
        The names are the data names: exp_ID, Choopper mode, Beam on, Beam off, Mean signal, ...
        The values are the data.
    """
    
    channels_str = ['ABES-'+str(a) for a in channels]
    data = {}
    txt = exp_ID
    

    print('Processing ' + exp_ID,flush=True)

    options = { 'Scaling':'Volt',
                'Amplitude calibration' : False
                }
    if (datapath is not None):
        options['Datapath'] = datapath
    
    try:
        options['State'] ={'Chop': 0, 'Defl': 0}
        options['Start'] = 0
        options['End'] = 0
        
        d_beam_on=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options=options
                                 )
        options['State'] ={'Chop': 1, 'Defl': 0}
        options['Start'] = 0
        options['End'] = 0
        d_beam_off=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options=options
                                 )           

        chopper_mode = d_beam_on.info['Chopper mode']
        on1,on2,on3 = d_beam_on.coordinate('Time')
        off1,off2,off3 = d_beam_off.coordinate('Time')
        beam_on_time = on3[1]-on2[1]
        beam_off_time = off3[1]-off2[1]
        period_time = beam_on_time + beam_off_time
        if (period_time > 3e-3):
            chop_str = "{:3.0f}-{:3.0f}[ms]".format(beam_on_time * 1e3, beam_off_time * 1e3)
        else:
            chop_str = "{:3.0f}-{:3.0f}[us]".format(beam_on_time * 1e6, beam_off_time * 1e6) 
        txt += ' ... Chopper: {:6s}({:s})'.format(chopper_mode,chop_str)
        data['exp_ID'] = exp_ID
        data['Chopper mode'] = chopper_mode
        data['Beam on time'] = beam_on_time
        data['Beam off time'] = beam_off_time
     
        if (period_time > 0.01):
            on_start = 1000
            on_end = -1000
            off_start = 1000
            off_end = -1000
            options['Resample'] = 1e4
        else:
            on_start = 0
            on_end = 0
            off_start = 0
            off_end = 0
            options['Resample'] = None
            
        d_beam_on=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options={'State':{'Chop': 0, 'Defl': 0},'Start':on_start,'End':on_end}
                                 )
        d_beam_off=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options={'State':{'Chop': 1, 'Defl': 0},'Start':off_start,'End':off_end}
                                 ) 
          
        sig = np.zeros(len(channels))
        for i in range(len(channels)):
            print('  Processing '+channels_str[i],flush=True)
            d=flap.get_data('W7X_ABES',
                            exp_id = exp_ID,
                            name = channels_str[i],
                            options = options,
                            coordinates = {'Time': timerange}
                            )
            data['Measurement start'] = float(d.info['Config']['APDCAM_starttime'])
            data['Measurement end'] = float(d.info['Config']['APDCAM_sampletime'] * d.info['Config']['APDCAM_samplenumber']) + data['Measurement start']
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
            if (test):
                plt.figure()
                d_on.plot(axes='Time')
                d_off.plot(axes='Time')
            d_on_data = d_on.data
            d_off_data = d_off.data
            ind = np.nonzero(np.logical_and(np.isfinite(d_off_data),
                                            np.isfinite(d_on_data)
                                            )
                             )[0]
            d_on_data = d_on_data[ind]
            d_off_data = d_off_data[ind]
            d = d_on_data - d_off_data
            if (i == 0):
                sig = np.zeros((len(d),len(channels)))
            sig[:,i] = d
        
        max_loc = np.unravel_index(np.argmax(sig), sig.shape)
        d_max = sig[max_loc]
        mean_signal = np.median(np.mean(sig,axis=1))
        
        # txt += ' ... Max:{:4.0f}[mV] '.format(d_max * 1000)
        txt += "... Mean signal: {:3d}[mV]".format(int(mean_signal*1000))
        txt += " ... Max: {:4d}[mV] at {:6s}/{:7.3f}s".format(int(d_max*1000),channels_str[max_loc[1]],round(d_on.get_coordinate_object("Time").values[max_loc[0]],3))
        data['Mean signal'] = mean_signal * 1000
        data['Max signal'] = d_max * 1000
        data['Max signal channel'] = channels_str[max_loc[1]]
        data['Max signal time'] = d_on.get_coordinate_object("Time").values[max_loc[0]]
        timescale = d_on.coordinate('Time')[0][ind]
        s = np.sum(sig,axis=1)
        if (np.max(s) <= 0):
            txt += '... Time range: ---'
            data['Good signal start'] = np.nan
            data['Good signal start'] = np.nan
        else:
            ind = np.nonzero(s >= np.max(s) * 0.1)[0]
            txt += ' ... Time range: ({:6.2f}-{:6.2f})[s]'.format(timescale[ind[0]], timescale[ind[-1]])
            data['Good signal start'] = timescale[ind[0]]
            data['Good signal start'] = timescale[ind[-1]]                                               
    except Exception as e:
        txt += ' --- {:s} ---'.format(str(e))
        data['exp_ID'] = exp_ID
        data['Chopper mode'] = ''
        data['Beam on time'] = np.nan
        data['Beam off time'] = np.nan
        data['Measurement start'] = np.nan
        data['Measurement end'] = np.nan
        data['Mean signal'] = np.nan
        data['Max signal'] = np.nan
        data['Max signal channel'] = np.nan
        data['Max signal time'] = np.nan 
        data['Good signal start'] = np.nan
        data['Good signal start'] = np.nan

    return txt,data
        
         
def exp_summaries(exp_ids,datapath=None,timerange=None,file='exp_summaries.txt',datafile='exp_summaries.dat'):  
    """
    Generate an experiment summary of a series of experiments.

    Parameters
    ----------
    exp_ids : str
        Experiment IDs. *, ? supported as in Unix regexp. Like 2018101?.*
    datapath : str, optional
        The data path. The default is None.
    timerange : list, optional
         The processing timerange. The default is None.
    file: str
        File name where to write the list    
    datafile: str
        The name of the output data file. This is a pickle file contaiining a dictionary. 
        The names are the data names: exp_ID, Choopper mode, Beam on, Beam off, Mean signal, ...
        The values are lists.
  
    Returns
    -------
    txts : list of texts
        The list of summary texts.
    data : dict
        Dictinoary with the data

    """
    
    def sort_by_name(entry):
        return entry.name
        
    if (datapath is not None):
        options = {'Datapath':datapath}
    else:
        options = {}
    default_options = {'Datapath':datapath}
    _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES')
    dp = _options['Datapath']
    regexp = exp_ids.replace('.','\.') 
    regexp = regexp.replace('*','.*') 
    txts = []
    data = {}
    id_list = list(os.scandir(dp))
    id_list.sort(key=sort_by_name)
    exp_list = []
    for idi in id_list:
        if (idi.is_dir()):
            if ((len(idi.name) == 12) and (idi.name[8] == '.') and (re.search(regexp,idi.name) is not None)):
                exp_list.append(idi.name)
    exp_list.sort()
#    print(exp_list)
    with open(file,"wt") as f:
        for exp in exp_list:
            txt,d = exp_summary(exp,datapath=dp,timerange=timerange)
            time.sleep(0.1)
            f.writelines(txt + '\n')
            f.flush()
            for k in d.keys():
                try:
                    data[k].append(d[k])  
                except KeyError:
                    data[k] = []
                    data[k].append(d[k]) 
            for k in data.keys():
                if (len(data[k]) != len(data[list(data.keys())[0]])):
                    raise ValueError("exp_summary did not return all keys.")
                                            
    if (datafile is not None):
        with open(datafile,"wb") as f:
            pickle.dump(data,f)
       
    return txts,data
        
#exp_summary('20230328.028')
