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
import numpy as np
import copy

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def w7x_summary(exp_ID,datapath='https://w7x-logbook.ipp-hgw.mpg.de/api//log/history/XP/'):
    """ Reads some information about the experiment. Needs to be run on the W7-X servers

    Parameters
    ----------
    exp_ID : string
        The exp ID, YYYMMDD.xxx
    datapath : str, optional
        The path were to find the log data. The default is 'https://w7x-logbook.ipp-hgw.mpg.de/api//log/history/XP/'.

    Returns
    -------
    dict:
        The disctionary with the data
        'Config': Main magnetic configuration

    """
    
    import requests, json

    path = datapath+'XP_'+exp_ID[:9]+"{:d}".format(int(exp_ID[9:]))
#    print(path)
    try:
        res = requests.get(path)
    except Exception as e:
        print("Not on W7-X machine?")
        raise RuntimeError("Cannot connect to W7-X logbook server.")

    if res.status_code == 200:
        result = {}
        result['exp_ID'] = exp_ID
        r = res.json()
        try:
            for c in r['hits']['hits'][-1]['_source']['tags']:
                if (c['description'] == 'Configuration main coils'):
                    result['Config'] = c['Main field']
        except IndexError:
            pass
        try:
            result['Config']
        except KeyError:
            result['Config'] = '???'
        return result
    
def summary_line(ABES_data=None, W7X_data=None):
    """ Returns a summary line created from ABES and W7X data. At least one should be present.
    

    Parameters
    ----------
    ABES_data : dict, optional
        ABES data returned by exp_summary. The default is None.
    W7X_data : dict, optional
        W7X data returned by w7x_summary. The default is None.

    Raises
    ------
    ValueError
        No data found or exp_ID in two data do not agree.

    Returns
    -------
    The text line.

    """
    
    
    if ((ABES_data is not None) and (W7X_data is not None)):
        if (ABES_data['exp_ID'] != W7X_data['exp_ID']):
            raise ValueError('Different exp_ID in ABES and W7X data.')
    if (ABES_data is not None):
        txt = ABES_data['exp_ID']
    elif (W7X_data is not None):
        txt = W7X_data['exp_ID']
    else:
        raise ValueError('No data.')
    if (W7X_data is not None):
        txt += "{:s}".format('(' + W7X_data['Config'] + ')')
        txt += '.' * (25 - (len(W7X_data['Config']) + 2))
    if (ABES_data is not None):
        if ((ABES_data['Chopper mode'] == '') or (ABES_data['Beam on time'] is None) or (ABES_data['Beam off time'] is None)):
            try:
                txt += ' --- {:s} ---'.format(ABES_data['error'])
            except KeyError:
                txt += ' --- {:s} ---'.format('Error')
        else:        
            period_time = ABES_data['Beam on time'] + ABES_data['Beam off time']
            if (period_time > 3e-3):
                chop_str = "{:3.0f}-{:3.0f}[ms]".format(ABES_data['Beam on time'] * 1e3, ABES_data['Beam off time'] * 1e3)
            else:
                chop_str = "{:3.0f}-{:3.0f}[us]".format(ABES_data['Beam on time'] * 1e6, ABES_data['Beam off time'] * 1e6) 
            txt += ' Chopper: {:6s}({:s})'.format(ABES_data['Chopper mode'],chop_str)
            txt += "... Mean signal: {:3d}[mV]".format(int(ABES_data['Mean signal']))
            txt += " ... Max: {:4d}[mV] at {:6s}/{:7.3f}s".format(int(ABES_data['Max signal']),ABES_data['Max signal channel'],round(ABES_data['Max signal time'],3))
            timetext = '... Time range: ---'
            if 'Good signal start' in ABES_data.keys():
                if (ABES_data['Good signal start'] != np.nan):                
                    timetext = ' ... Time range: ({:6.2f}-{:6.2f})[s]'.format(ABES_data['Good signal start'], ABES_data['Good signal end'])
            
            txt += timetext
        # if we're in extended mode:
        if "APD voltage" in ABES_data.keys():
            txt += f" ... APD voltage: {(ABES_data['APD voltage'])}"
            txt += f" ... MPPC voltage: {(ABES_data['MPPC voltage'])}"
            txt += f" ... Clock source: {(ABES_data['Clock source'])}"
            if "Camera" not in txt:
                txt += f" ... Chopper period: {(ABES_data['Chopper period'])}"
    return txt
 
def exp_summary(exp_ID,timerange=None,datapath=None,channels=range(10,30),test=False, extended=False):
    """
    Return a single line summary of the experiment and a dictionaty with data.

    Parameters
    ----------
    exp_ID : string
        The exp ID, YYYMMDD.xxx
    timerange : list, optional
        The processing timerange. The default is None.
    datapath : string, optional
        The datapath. The default is None.
    extended : bool, optional
        Whether to show some additional information  about the APDCAM configuration
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
#    txt = exp_ID
    

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
#        txt += ' ... Chopper: {:6s}({:s})'.format(chopper_mode,chop_str)
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
        
        beam_on_state = {'State':{'Chop': 0, 'Defl': 0},'Start':on_start,'End':on_end}
        d_beam_on=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options=beam_on_state
                                 )
        beam_off_state = {'State':{'Chop': 1, 'Defl': 0},'Start':on_start,'End':on_end}
        d_beam_off=flap.get_data('W7X_ABES',
                                 exp_id=exp_ID,
                                 name='Chopper_time',
                                 options=beam_off_state
                                 ) 
          
        sig = np.zeros(len(channels))
        maxraw = 0
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
                plt.show()
                plt.pause(0.1)
            d_on_data = d_on.data
            d_off_data = d_off.data
            ind = np.nonzero(np.logical_and(np.isfinite(d_off_data),
                                            np.isfinite(d_on_data)
                                            )
                             )[0]
            d_on_data = d_on_data[ind]
            d_off_data = d_off_data[ind]
            d_clean = d_on_data - d_off_data
            maxraw = max([maxraw,
                          max(d_on_data),
                          max(d_off_data)])
            if (i == 0):
                sig = np.zeros((len(d_clean),len(channels)))
            sig[:,i] = d_clean
        
        max_loc = np.unravel_index(np.argmax(sig), sig.shape)
        d_max = sig[max_loc]
        mean_signal = np.mean(np.mean(sig,axis=1))
        
#        txt += "... Mean signal: {:3d}[mV]".format(int(mean_signal*1000))
#        txt += " ... Max: {:4d}[mV] at {:6s}/{:7.3f}s".format(int(d_max*1000),channels_str[max_loc[1]],round(d_on.get_coordinate_object("Time").values[max_loc[0]],3))
        data['Mean signal'] = mean_signal * 1000
        data['Max signal'] = d_max * 1000
        data['Max signal channel'] = channels_str[max_loc[1]]
        data['Max signal time'] = d_on.get_coordinate_object("Time").values[max_loc[0]]
        timescale = d_on.coordinate('Time')[0][ind]
        s = np.sum(sig,axis=1)
        if (np.max(s) <= 0):
#            txt += '... Time range: ---'
            data['Good signal start'] = np.nan
            data['Good signal end'] = np.nan
        else:
            ind = np.nonzero(s >= np.max(s) * 0.1)[0]
            data['Good signal start'] = timescale[ind[0]]
            data['Good signal end'] = timescale[ind[-1]]  

        if extended:
            data["APD voltage"] = d.info['Config']['APDCAM_bias1']
            data["MPPC voltage"] = d.info['Config']['APDCAM_bias2']
            data["Clock source"] = d.info['Config']['APDCAM_clock_source']
            if chopper_mode != "Camera":
                data["Chopper period"] = d.info['Config']['Chopper period']
            data["Max raw signal"] = maxraw
#            txt += ' ... Time range: ({:6.2f}-{:6.2f})[s]'.format(timescale[ind[0]], timescale[ind[-1]])
                                             

    except Exception as e:
#        txt += ' --- {:s} ---'.format(str(e))
        data['error'] = str(e)
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
        data['Good signal end'] = np.nan
        if extended:
            data["APD voltage"] = np.nan
            data["MPPC voltage"] = np.nan
            data["Clock source"] = np.nan
            data["Chopper period"] = np.nan

    txt = summary_line(ABES_data=data, W7X_data=None)
    return txt,data
        
         
def exp_summaries(exp_ids,datapath=None,timerange=None,file='exp_summaries.txt',datafile='exp_summaries.dat',abes=True,w7x=False):  
    """
    Generate an experiment summary of a series of experiments.

    Parameters
    ----------
    exp_ids : str
        Experiment IDs. *, ? supported as in Unix regexp. Like 2018101?.*. Will not be used if
        abes is False.
    datapath : str, optional
        The data path. The default is None.
    timerange : list, optional
         The processing timerange. The default is None.
    file: str
        File name where to write the list. This will be generated from the updated data.    
    datafile: str
        The name of the output data file. This is a pickle file contaiining a dictionary. 
        The names are the data names: exp_ID, Choopper mode, Beam on, Beam off, Mean signal, ...
        The values are lists. Additionally it may contain W7-X data.
    abes: bool, optional, The default is True
        If True will read ABES data and write a new datafile.
    w7x: bool, optional, the default is False
        If True will read W7-X data and add the information to the experiments in the 
        datafile or to the ABES data being read, if abes is True. 
  
    Returns
    -------
    txts : list of texts
        The list of summary texts.
    data : dict
        Dictinoary with the data

    """
    
    def sort_by_name(entry):
        return entry.name
        
    txts = []
    if (abes):
        if (datapath is not None):
            options = {'Datapath':datapath}
        else:
            options = {}
        default_options = {'Datapath':datapath}
        _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES')
        dp = _options['Datapath']
        regexp = exp_ids.replace('.','\\.') 
        regexp = regexp.replace('*','.*') 
        data = {}
        id_list = list(os.scandir(dp))
        id_list.sort(key=sort_by_name)
        exp_list = []
        for idi in id_list:
            if (idi.is_dir()):
                if ((len(idi.name) == 12) and (idi.name[8] == '.') and (re.search(regexp,idi.name) is not None)):
                    exp_list.append(idi.name)
        exp_list.sort()
    else:
        with open(datafile,"rb") as f:
            data_orig = pickle.load(f)
            data = copy.deepcopy(data_orig)
        exp_list = data['exp_ID']
    with open(file,"wt") as f:
        for i_exp,exp in enumerate(exp_list):
            print(exp)
            if (abes):
                txt,data_abes = exp_summary(exp,datapath=dp,timerange=timerange)
            else:
                data_abes = {}
                for k in data_orig.keys():
                    try:
                        data_abes[k] = data_orig[k][i_exp]
                    except IndexError:
                        aaa=0
            d = copy.deepcopy(data_abes)
            if (w7x):
                data_w7x = w7x_summary(exp)
                d.update(data_w7x) 
            else:
                data_w7x = None
            txt = summary_line(ABES_data=data_abes, W7X_data=data_w7x)
            txts.append(txt)
            f.writelines(txt + '\n')
            f.flush()
            if (abes):
                for k in d.keys():
                    try:
                        data[k].append(d[k])
                    except KeyError:
                        data[k] = [d[k]]
            else:
                for k in data_w7x.keys():
                    if (i_exp == 0):
                        try:
                            data[k]
                        except KeyError:
                            data[k] = [None] * len(data['exp_ID'])                 
                    data[k][i_exp] = data_w7x[k]  
                                                            
    if (datafile is not None):
        with open(datafile,"wb") as f:
            pickle.dump(data,f)
       
    return txts,data

def find_experiment(search_dict,datafile='exp_summaries.dat',list_keys='True'):
    """Find an experiment in the datafile fulfilling certain criteria.

    Parameters
    ----------
    search_dict : dict
        Dictionary with search criteria. Keys should be identical to keys in datafile.
        If a value for a key:
            is a single element then experiments will be searched where the value is equal to this.
            is a list then experiments will be serached where the value is between (equal) to the two values in the list. 
    datafile : str, optional, deafult is 'exp_summaries.dat'
        The name of the datafile written by exp_summaries. 
    list_keys : bool
        Lists the keys which can can be used in the search.

    Raises
    ------
    ValueError
        Imvalid key or single element list in value.

    Returns
    -------
    array-like
        The experiment IDs fulfilling the requirement.

    Use this module like this to find experiments with timed fast chopper with signal above 30 mV::
        
        print(find_experiment({"Chopper mode":'Timed',
                               'Beam on time':[1e-6,1e-5],
                               'Mean signal':[30,3000]
                               },datafile=df)
              )        
        
    """

    with open(datafile,"rb") as f:
        df = pickle.load(f)
    if (list_keys):
        print(list(df.keys()))
    index = np.array(list(range(len(df['exp_ID']))))
    for k in search_dict.keys():
        if (type(search_dict[k]) is list):
            try:
                ind = np.nonzero(np.logical_and(np.array(df[k])[index] >= search_dict[k][0],
                                                np.array(df[k])[index] <= search_dict[k][1]
                                                )
                                 )[0]
                index = index[ind]
            except KeyError:
                txt = ''
                for kk in df.keys():
                    txt += ', ' + kk
                raise ValueError('Key "{:s}" is not present in datafile. \n Valid keys are:{:s}'.format(k,txt))
            except IndexError:
                raise ValueError('Error searching for range of data in "{:s}"'.format(k))
        else:
            try:
                index = index[np.nonzero((np.array(df[k])[index] == search_dict[k]))[0]]
            except KeyError:
                txt = ''
                for kk in df.keys():
                    txt += ', ' + kk
                raise ValueError('Key "{:s}" is not present in datafile. \n Valid keys are:{:s}'.format(k,txt))
    return np.array(df['exp_ID'])[index]
        
if __name__ == '__main__':
    pass
    # out = exp_summary('20250312.099')
    # df = 'c:/Users/Zoletnik/OneDrive - energia.mta.hu/Megosztott dokumentumok - FPL/Projects/Experiments/W7-X/ABES/op21/Log/exp_summaries_2023.dat'
    # print(find_experiment({"Chopper mode":'Timed',
    #                         'Beam on time':[1e-6,1e-5],
    #                         'Mean signal':[40,300]
    #                         },datafile=df)
    #       )
