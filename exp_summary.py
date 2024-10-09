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
    timerange : listE, optional
        The processing timerange. The default is None.
    datapath : string, optional
        The datapath. The default is None.

    Returns
    -------
    txt : string
        The experiment summary.

    """
    
    channels_str = ['ABES-'+str(a) for a in channels]
    txt = exp_ID
    print('Processing ' + exp_ID,flush=True)

    options = { 'Scaling':'Volt',
                'Amplitude calibration' : False
                }
    if (datapath is not None):
        options['Datapath'] = datapath
    
    try:
        d_beam_on=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options={'State':{'Chop': 0, 'Defl': 0},'Start':0,'End':0}
                                 )
        d_beam_off=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options={'State':{'Chop': 1, 'Defl': 0},'Start':0,'End':0}
                                 )           

        chopper_mode = d_beam_on.info['Chopper mode']
        on1,on2,on3 = d_beam_on.coordinate('Time')
        off1,off2,off3 = d_beam_on.coordinate('Time')
        beam_on_time = on3[1]-on2[1]
        beam_off_time = off3[1]-off2[1]
        period_time = beam_on_time + beam_off_time
        if (period_time > 3e-3):
            chop_str = "{:3.0f}-{:3.0f}[ms]".format(beam_on_time * 1e3, beam_off_time * 1e3)
        else:
            chop_str = "{:3.0f}-{:3.0f}[us]".format(beam_on_time * 1e6, beam_off_time * 1e6) 
        txt += ' ... Chopper:{:6s}({:s})'.format(chopper_mode,chop_str)
     
        if (period_time > 0.01):
            on_start = 1000
            on_end = -1000
            off_start = 1000
            off_end = -1000
            options['Resample'] = 1e4
        else:
            on_start = 0
            on_end = -3
            off_start = 0
            off_end = -3
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
            
        d_max = np.max(sig)
        txt += ' ... Max:{:4.0f}[mV] '.format(d_max * 1000)
        timescale = d_on.coordinate('Time')[0][ind]
        s = np.sum(sig,axis=1)
        if (np.max(s) <= 0):
            txt += '... Time range: ---'
        else:
            ind = np.nonzero(s >= np.max(s) * 0.1)[0]
            txt += ' ... Time range:({:6.2f}-{:6.2f})[s]'.format(timescale[ind[0]], timescale[ind[-1]])
    except Exception as e:
        txt += ' --- {:s} ---'.format(str(e))
    return txt
        
         
def exp_summaries(exp_ids,datapath=None,timerange=None,file='exp_summaries.txt'):  
    """
    Generate an experiment summary of a series of experiments.

    Parameters
    ----------
    exp_ids : string
        Experiment IDs. *, ? supported as in Unix regexp. Like 2018101?.*
    datapath : string, optional
        The data path. The default is None.
    timerange : listE, optional
         The processing timerange. The default is None.
    file: string
        File name where to write the list       
  
    Returns
    -------
    txts : list of texts
        The list of summary texts.

    """
    
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
    id_list = os.scandir(dp)
    exp_list = []
    for idi in id_list:
        if (idi.is_dir()):
            if ((len(idi.name) == 12) and (idi.name[8] == '.') and (re.search(regexp,idi.name) is not None)):
                exp_list.append(idi.name)
    exp_list.sort()
#    print(exp_list)
    with open(file,"wt") as f:
        for exp in exp_list:
            txt = exp_summary(exp,datapath=dp,timerange=timerange)
            time.sleep(2)
            f.writelines(txt + '\n')
            f.flush()
    return txts
        
#exp_summary('20230328.028')
