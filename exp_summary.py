# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 21:16:22 2024

@author: Zoletnik
"""
import numpy as np
import os
import re

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def exp_summary(exp_ID,timerange=None,datapath=None):
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
    

    options = {'Resample':1e3,
               'Scaling':'Volt',
               'Amplitude calibration' : False
               }
    if (datapath is not None):
        options['Datapath'] = datapath

    txt = exp_ID
    print('Processing ' + exp_ID)
    try:
        d=flap.get_data('W7X_ABES',
                        exp_id = exp_ID,
                        name ='ABES-[1-40]',
                        options = options,
                        coordinates = {'Time': timerange}
                        )
        
        d_max = d.slice_data(summing={'Time':'Max','Channel':'Max'})
        txt += ' ... Max:{:4.0f}[mV] '.format(d_max.data * 1000)
        
        d_sum = d.slice_data(summing={'Channel': 'Mean'})
        timescale = d_sum.coordinate('Time')[0]
        t_on= timescale[np.nonzero(d_sum.data > np.max(d_sum.data) * 0.1)[0]]
        txt += ' ... Signal timerange:( {:6.2f}-{:6.2f})[s]'.format(t_on[0],t_on[-1])
        try:
            d=flap.get_data('W7X_ABES',
                            exp_id = exp_ID,
                            name ='ABES-[1-40]',
                            options = {'Resample':1e3,
                                       'Scaling':'Volt',
                                       'Amplitude calibration' : True
                                       },
                            coordinates = {'Time': timerange}
                            )
        
            d_dist = d.slice_data(slicing={'Time':flap.Intervals(t_on[0],t_on[-1])},
                                  summing={'Time':'Mean'}
                                  )
            chscale = d_dist.coordinate('Channel',options={'Change':True})[0]
            ch_on= chscale[np.nonzero(d_dist.data > np.max(d_dist.data) * 0.1)[0]]
            txt += ' ... Channel range:({:6.2f}-{:6.2f})'.format(ch_on[0],ch_on[-1])
        except Exception:
            txt += ' ... Channel range:  (no calib)  '
    except Exception:
        txt += ' ... No data'
    return txt
         
def exp_summaries(exp_ids,datapath=None,timerange=None):  
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
    for idi in id_list:
        if (idi.is_dir()):
            if ((len(idi.name) == 12) and (idi.name[8] == '.') and (re.search(regexp,idi.name) is not None)):
                txts.append(exp_summary(idi.name,datapath=dp,timerange=timerange))
    return txts
        
