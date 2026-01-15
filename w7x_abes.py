# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:23:49 2018

@author: Zoletnik

This is the flap module for W7-X alkali BES diagnostic
"""

import os
import time
import warnings
from decimal import Decimal
import gc
import psutil
import numpy as np
import copy
import h5py
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial

import flap
from . import spatcal

if (flap.VERBOSE):
    print("Importing flap_w7x_abes")

def abes_get_config(xml):
    retval = {}
    try:
        retval['version'] = xml.head.attrib['Version']
    except:
        retval['version'] = '1.0'
    try:
        retval['ShotID'] = xml.head.attrib['ShotID']
    except KeyError:
        raise ValueError("Invalid config file head format.")
    if (retval['version'] != '1.0'):    
        try:
            retval['Time'] = xml.head.attrib['Time']
            retval['Date'] = xml.head.attrib['Date']
        except KeyError:
            raise ValueError("Invalid config file head format.")
    if (retval['version'] != '1.0'):    
        retval['TriggerTime'] = Decimal(xml.get_element('System','TriggerTime')['Value'])
    else:
        retval['TriggerTime'] = Decimal(xml.get_element('APDCAM','Trigger')['Value'])
    # Micrometer settings for spatial calibration
    if xml.get_element('System','APD_H-Micrometer')['Unit'] == 'mm' and \
       xml.get_element('System','APD_V-Micrometer')['Unit'] == 'mm':
        retval['H-Micrometer'] = float(xml.get_element('System', 'APD_H-Micrometer')['Value'])
        retval['V-Micrometer'] = float(xml.get_element('System', 'APD_V-Micrometer')['Value'])
    else:
        raise ValueError('H-Micrometer and V-Micrometer units should be in mm')
    retval['APDCAM_state'] = int(xml.get_element('APDCAM','State')['Value'])
    if (retval['APDCAM_state'] == 1):
        ADCDiv = Decimal(xml.get_element('APDCAM', 'ADCDiv')['Value'])
        ADCMult = Decimal(xml.get_element('APDCAM', 'ADCMult')['Value'])
        retval['APDCAM_f_ADC'] = Decimal(20e6) * ADCMult / ADCDiv
        samplediv = Decimal(xml.get_element('APDCAM', 'Samplediv')['Value'])
        retval['APDCAM_f_sample'] = retval['APDCAM_f_ADC'] / samplediv
        retval['APDCAM_sampletime'] = Decimal(1.) / retval['APDCAM_f_sample']
        retval['APDCAM_samplenumber'] = int(xml.get_element('APDCAM', 'SampleNumber')['Value'])
        retval['APDCAM_bits'] = int(xml.get_element('APDCAM', 'Bits')['Value'])
        trigger = Decimal(xml.get_element('APDCAM', 'Trigger')['Value'])
        if (trigger < 0):
            trigger = Decimal(0)
        retval['APDCAM_starttime'] = trigger + retval['TriggerTime']
        mask1 = int(xml.get_element('APDCAM', 'ChannelMask1')['Value'],16)
        mask2 = int(xml.get_element('APDCAM', 'ChannelMask2')['Value'],16)
        mask3 = int(xml.get_element('APDCAM', 'ChannelMask3')['Value'],16)
        mask4 = int(xml.get_element('APDCAM', 'ChannelMask4')['Value'],16)
        chmask = mask1 + (mask2 << 32) + (mask3 << 64) + (mask4 << 96)
        retval['APDCAM_chmask'] = chmask
        retval['APDCAM_bias1'] = float(xml.get_element('APDCAM','DetectorBias1')['Value'])
        retval['APDCAM_bias2'] = float(xml.get_element('APDCAM','DetectorBias2')['Value'])
        retval['APDCAM_det_temp'] = float(xml.get_element('APDCAM','DetectorTemp')['Value'])
        if (retval['version'] != '1.0'):  
            try:
                s = int(xml.get_element('APDCAM','SystemClockSource')['Value'])
            except:
                s = 0
        else:
            s = int(xml.get_element('System','SystemClockSource')['Value'])
        if (s == 1):
            retval['APDCAM_clock_source'] = 'External'
        else:
            retval['APDCAM_clock_source'] = 'Internal'

    chopmode = int(xml.get_element('Chopper','Mode')['Value'])
    clk = xml.get_element('Chopper','BaseClockFrequency')
    clk_freq = Decimal(clk['Value'])
    if (clk['Unit'] == 'Hz'):
        pass
    elif (clk['Unit'] == 'kHz'):
        clk_freq = clk_freq * Decimal(1000)
    elif (clk['Unit'] == 'MHz'):
        clk_freq = clk_freq * Decimal(1E6)
    else:
        raise ValueError("Unknown chopper clock frequency unit.")
    retval['Chopper clock'] = clk_freq
    if (chopmode == 1):
        retval['Chopper mode'] = 'Camera'
        retval['CMOS exptime'] = Decimal(xml.get_element('CMOS','Exptime')['Value'])
        retval['CMOS frametime'] = Decimal(xml.get_element('CMOS','Frametime')['Value'])
        retval['CMOS frame number'] = int(xml.get_element('CMOS','FrameNumber')['Value'])
    else:
        retval['Chopper mode'] = 'Timed'
    if (retval['version'] != '1.0'):    
        sch = xml.get_element('Chopper','SchemeFileContents')['Value']
        if (chopmode == 0):
            retval['Chopper period'] = Decimal(xml.get_element('Chopper','PeriodTime')['Value']) \
                                         /Decimal(1000000)
    else:
        pol_enable = int(xml.get_element('Chopper','PolEnable')['Value'])
        tor_enable = int(xml.get_element('Chopper','TorEnable')['Value'])
        if ((tor_enable == 1) and (pol_enable == 1)):
            raise ValueError("Toroidal and Poloidal chopper is enabled at the same time. Cannot handle this.")
        sch = ""
        if (tor_enable == 1):
            if (int(xml.get_element('Chopper','TorStartTimeCLK')['Value'] )!= 0) :
                raise ValueError("Toroidal chopper start time not 0, cannot handle.")
            retval['Chopper period'] = (int(xml.get_element('Chopper','TorOnTimeCLK')['Value']) 
                                       + int(xml.get_element('Chopper','TorOffTimeCLK')['Value'])) / clk_freq
            if (retval['Chopper mode'] == 'Camera'):
                sch = "[General] <NL>Name = Cam_Chop <NL>[Frame 1] <NL>Chopper = 1 <NL>Deflection = 0 <NL>[Frame 2] <NL>Chopper = 0 <NL>Deflection = 0" 
            else:
                period = int(xml.get_element('Chopper','TorOnTimeCLK')['Value']) + int(xml.get_element('Chopper','TorOnTimeCLK')['Value'])
                phase1 = float(xml.get_element('Chopper','TorOnTimeCLK')['Value']) / float(period) * 100
                phase2 = float(xml.get_element('Chopper','TorOffTimeCLK')['Value']) / float(period) * 100
                sch = "[General] <NL>Name=Simple fast <NL>"
                sch += "[Phase 1] <NL>Length = "+str(int(round(phase1)))+" <NL>Chopper = 0 <NL>Deflection = 0 <NL>"
                sch += "[Phase 2] <NL>Length = "+str(int(round(phase2)))+" <NL>Chopper = 1 <NL>Deflection = 0 <NL>"
        if (pol_enable == 1):
            if (int(xml.get_element('Chopper','PolStartTimeCLK')['Value']) != 0) :
                raise ValueError("Poloidal chopper start time not 0, cannot handle.")
            retval['Chopper period'] = (int(xml.get_element('Chopper','PolOnTimeCLK'))['Value'] 
                                       + int(xml.get_element('Chopper','PolOffTimeCLK')['Value'])) / clk_freq
            if (retval['Chopper mode'] == 'Camera'):
                sch = "[General] <NL>Name = Cam_Chop <NL>[Frame 1] <NL>Chopper = 0 <NL>Deflection = 1 <NL>[Frame 2] <NL>Chopper = 0 <NL>Deflection = 0" 
            else:
                period = int(xml.get_element('Chopper','PolOnTimeCLK')['Value']) + int(xml.get_element('Chopper','PolOnTimeCLK')['Value'])
                phase1 = float(xml.get_element('Chopper','PolOnTimeCLK')['Value']) / float(period) * 100
                phase2 = float(xml.get_element('Chopper','PolOffTimeCLK')['Value']) / float(period) *100
                sch = "[General] <NL>Name=Simple fast <NL>"
                sch += "[Phase 1] <NL>Length = "+str(int(round(phase1)))+" <NL>Chopper = 0 <NL>Deflection = 0 <NL>"
                sch += "[Phase 2] <NL>Length = "+str(int(round(phase2)))+" <NL>Chopper = 0 <NL>Deflection = 1 <NL>"
    sch.replace("\n","")
    sch = sch.split('<NL>')
    retval['Chopper scheme'] = sch
    ch_list = []
    fibre_list = []
    det_type_list = []
    adc_list = []
    for i in range(64):
        try:
            ch = xml.get_element('Optics', 'ADC' + str(i + 1))['Value']
            adc_list.append(i + 1)
            ch_arr = ch.split('-')
            if ((ch_arr[0][0] >= '1') and (ch_arr[0][0] <= '9')):
                ch_arr[0] = 'ABES-' + ch_arr[0]
            ch_list.append(ch_arr[0])
            fibre_list.append(ch_arr[1])
            if (len(ch_arr) > 2):
                det_type_list.append(ch_arr[2])
            else:
                det_type_list.append(' ')
        except Exception:
            pass

    def sortkey(in_string):
        signal_name = in_string[0]
        if (signal_name[0:5] == 'ABES-'):
            if (len(signal_name[5:]) == 1):
                signal_name = 'ABES-0' + signal_name[5:]
        return signal_name


    # Ordering the channel list
    ch_list, adc_list, fibre_list,det_type_list = zip(*sorted(zip(ch_list,
                                                                  adc_list,
                                                                  fibre_list,
                                                                  det_type_list),
                                                      key=sortkey))
    retval['ADC_list'] = list(adc_list)
    retval['signal_list'] = list(ch_list)
    retval['fibre_list'] = list(fibre_list)
    retval['det_type_list'] = list(det_type_list)

    return retval

def calibrate(data_arr, signal_proc, read_range, exp_id=None, options=None):
    """ Do calibration of the data if it is requested in options
    """
    try:
        calib = options['Amplitude calibration']
    except KeyError:
        calib = False
    if (type(calib) is not bool):
        raise ValueError("Invalid calibration setting. Should be True or False.")
    if (not calib):
        return data_arr, None, None

    try:
        calibration_path = options['Amplitude calib. path']
    except KeyError:
        calibration_path = 'cal'

    if (options['Amplitude calib. file'] is not None):
        calibration_file = options['Amplitude calib. file']
        index_start = [0]
        index_stop = [data_arr.size[1]]
        cal_files = [calibration_file]
    else:
        fn = os.path.join(calibration_path, exp_id + '.cal')
        with open(fn, 'r', encoding='ascii') as infile:
            # These lists will contain the information from the cal file
            index_start = []
            index_stop = []
            cal_files = []
            cal_factors = []
            for line in infile:
                line = line.split()
                if (len(line) == 1):
                    index_start = [0]
                    index_stop = [data_arr.shape[0]]
                    cal_files = [line[0]]
                else:
                    try:
                        cal_timerange = [float(line[1]), float(line[2])]
                    except (IndexError, ValueError):
                        infile.close()
                        raise ValueError("Invalid shot calibration file: "+fn)
                    if ((read_range[0] > cal_timerange[1]) or   \
                        (read_range[1] < cal_timerange[0])):
                        # The validity of this calibration is outside the data range
                        continue
                    if (cal_timerange[0] > read_range[0]) :
                        s = float(cal_timerange[0] - read_range[0])  \
                            / (read_range[1] - read_range[0]) \
                            * data_arr.shape[0]
                    else:
                        s = 0
                    if (cal_timerange[1] < read_range[1]) :
                        e = float(cal_timerange[1] - read_range[0]) \
                            / (read_range[1] - read_range[0])  \
                            * data_arr.shape[0]
                    else:
                        e = data_arr.shape[0]
                    index_start.append(s)
                    index_stop.append(e)
                    cal_files.append(line[0])
    #Sorting the validity intervals according to the start index
    z = zip(index_start,index_stop,cal_files)
    z = list(z)
    z.sort()
    index_start,index_stop,cal_files = zip(*z)
    # Going through the intervals and looking for breaks in the validity
    prev_end = 0
    for i in range(len(index_start)):
        if (index_start[i] > prev_end):
            break
        prev_end = index_stop[i]
        if (prev_end >= data_arr.shape[0]):
            index_start = index_start[0:i+1]
            index_stop = index_stop[0:i+1]
            cal_files = cal_files[0:i+1]
            break
    if (prev_end < data_arr.shape[0]):
        raise ValueError("There is no valid calibration for part of the time inteval.")
    # Reading calibration data
    calfac = []
    calfac_err = []
    cal_channels = []
    for fn in cal_files:
        try:
            fname = os.path.join(calibration_path,fn)
            f = h5py.File(fname, 'r')
        except Exception as e:
            raise OSError("Error reading calibration from {:s}: {:s}".format(fname,str(e)))
        try:
            calfac.append(copy.deepcopy(np.array(f['Calibration_factors'])))
            calfac_err.append(copy.deepcopy(np.array(f['Calibration_factor_errors'])))
            cal_channels.append(copy.deepcopy(list(f['Channels'])))
        except Exception as e:
            raise e
        f.close()

    signal_coord = flap.Coordinate(name='Signal name',
                                   unit='n.a.',
                                   mode=flap.CoordinateMode(equidistant=False),
                                   shape=len(calfac_err[0]),
                                   values=np.array([chan.decode('utf-8').split(' ')[0] for chan in cal_channels[0]]),
                                   dimension_list=[0])
    calfac_err_dataobject = flap.DataObject(data_array=calfac_err[0]/calfac[0],
                                            exp_id=exp_id,
                                            data_shape=calfac_err[0].shape,
                                            coordinates = [signal_coord],
                                            data_title = 'Relative calibration factor error')

    # Doing the calibration
    data_err = np.zeros(data_arr.shape)
    if (data_arr.dtype.kind != 'f'):
            data_arr = float(data_arr)
    for i in range(len(index_start)):
        # Collecting the calibration factors andd errors
        calfac_act = np.empty(len(signal_proc),dtype=float)
        calfac_act_err = np.empty(len(signal_proc),dtype=float)
        for i_ch in range(len(signal_proc)):
            try:
                sig = signal_proc[i_ch].encode("ascii")
                if (len(sig) < 7) and sig:
                    sig += b' '
                i_cal = cal_channels[i].index(sig)
            except IndexError:
                if (signal_proc[i_ch][0:5] == 'ABES-'):
                    raise ValueError("No calibration data for signal "+signal_proc[i_ch])
            if (data_arr.ndim == 2):
                data_arr[index_start[i]:index_stop[i], i_ch] /= calfac[i][i_cal].astype("float")
                data_err[index_start[i]:index_stop[i], i_ch] = data_arr[index_start[i]:index_stop[i], i_ch] / calfac[i][i_cal] * calfac_err[i][i_cal]
            else:
                data_arr[index_start[i]:index_stop[i]] /= calfac[i][i_cal].astype("float")
                data_err[index_start[i]:index_stop[i]] = data_arr[index_start[i]:index_stop[i]] / calfac[i][i_cal] * calfac_err[i][i_cal]
    return data_arr, data_err, calfac_err_dataobject

def calculate_amplitude_calibration(shotID, options={}):

    options_default={'Time window': None,
                     'Sample window': None,
                     'Amplitude calib. path': flap.config.get('Module W7X_ABES', 'Amplitude calib. path'),
                     'Overwrite': False,
                     'Partial intervals': False,
                     'Chop delay': read_chopshift(shotID),
                     'Plot': False}
    options = {**options_default, **options}
    options['Amplitude calibration'] = False
    chop_delay = copy.deepcopy(options['Chop delay'])
    plot_data = options['Plot']
    del options['Chop delay']
    del options['Plot']
    
    filename = os.path.join(options['Amplitude calib. path'], shotID+'.cld')
    if os.path.exists(filename) is True and options['Overwrite'] is False:
            raise ValueError('Calibration data ('+filename+") already exists. To overwrite, set options['Overwrite']=True.")
    del options['Overwrite']

    if (options['Time window'] is not None):
        del options['Sample window']
        data_range = copy.deepcopy(options['Time window'])
        del options['Time window']
        light_profile = w7x_abes_get_data(exp_id=shotID, data_name='ABES-*', options=options,
                                        coordinates=[flap.Coordinate(name='Time',c_range=data_range)])
    elif (options['Sample window'] is not None):
        del options['Time window']
        data_range = copy.deepcopy(options['Sample window'])
        del options['Sample window']
        light_profile = w7x_abes_get_data(exp_id=shotID, data_name='ABES-*', options=options,
                                        coordinates=[flap.Coordinate(name='Sample',c_range=data_range)])
    else:
        raise ValueError('A "Sample window" or a "Time window" vector of length 2 should be given in'+
                         ' the options dictionary.')

    # Removing the background
    beam_on = proc_chopsignals(dataobject=light_profile,
                       on_options={'Start delay': chop_delay[0], 'End delay': chop_delay[1], 'Partial intervals':options['Partial intervals']},
                       off_options={'Start delay': chop_delay[0], 'End delay': chop_delay[1],  'Partial intervals':options['Partial intervals']},
                       options=options, test=False)

    del light_profile
    beam_on = beam_on.slice_data(summing={'Time': 'Mean'})
    if plot_data is True:
        plt.subplot(2,1,1)
        beam_on.plot(axes='Channel')
        plt.subplot(2,1,2)
        plt.plot(beam_on.error/beam_on.data)
    
    # Getting and formatting the relevant data
    channel_names = beam_on.get_coordinate_object('Signal name').values
    calfac = beam_on.data
#    calfac = [cf.encode('ascii') for cf in calfac]
    calfac_err = beam_on.error
#    calfac_err = [cfe.encode('ascii') for cfe in calfac_err]
    sorted_cal = np.asarray(sorted(zip(channel_names, calfac, calfac_err)))
    calfac = np.asarray(sorted_cal[:,1]).astype(float)
    calfac_err = sorted_cal[:,2].astype(float)
    channel_names = sorted_cal[:,0]
    channel_names = [(ch_name+' ').encode('ascii') if len(ch_name)<7 else ch_name.encode('ascii') for ch_name in channel_names]
    
    # Saving the calibration data
    with h5py.File(filename+'new', 'w', libver='earliest') as f:
        f.create_dataset('Calibration_factors', data=calfac, dtype=float)
        f.create_dataset('Calibration_factor_errors', data=calfac_err, dtype=float)
        f.create_dataset('Channels', data=channel_names)
        if os.path.exists(filename) is True:
            os.remove(filename)
    
    os.rename(filename+'new', filename)

def process_chopper_setup(config):
    sch = config['Chopper scheme']
    if (config['Chopper mode'] == 'Camera'):
        chop = []
        defl = []
        frame_number = 0
        line_counter = 2
        while True:
            for il in range(line_counter,len(sch)):
                if (sch[il].strip() == '[Frame '+str(len(chop)+1)+']'):
                    break
            else:
                break
            start_n = len(chop)
            for il1 in range(il+1,len(sch)):
                if (len(sch[il1].strip()) == 0):
                    continue
                if (sch[il1].strip()[0] == '['):
                    break
                par  = sch[il1].split('=')
                if (par[0].strip() == 'Chopper'):
                    chop.append(int(par[1].strip()))
                if (par[0].strip() == 'Deflection'):
                    defl.append(int(par[1].strip()))
            if ((len(chop) == start_n) or (len(chop) != len(defl))):
                raise ValueError("Invalid chopper file contents.")
            line_counter = il1
        length = [100./len(chop)]*len(chop)
    else:
        chop = []
        defl = []
        length = []
        frame_number = 0
        line_counter = 2
        while True:
            for il in range(line_counter,len(sch)):
                if (sch[il].strip() == '[Phase '+str(len(chop)+1)+']'):
                    break
            else:
                break
            start_n = len(chop)
            for il1 in range(il+1,len(sch)):
                if (len(sch[il1].strip()) == 0):
                    continue
                if (sch[il1].strip()[0] == '['):
                    break
                par  = sch[il1].split('=')
                if (par[0].strip() == 'Chopper'):
                    chop.append(int(par[1].strip()))
                if (par[0].strip() == 'Length'):
                    length.append(float(par[1].strip()))
                if (par[0].strip() == 'Deflection'):
                    defl.append(int(par[1].strip()))
            if ((len(chop) == start_n) or (len(chop) != len(defl)) or (len(chop) != len(length))):
                raise ValueError("Invalid chopper file contents.")
            line_counter = il1
    return length, chop, defl

def chopper_timing_data_object(config, options, read_samplerange=None):
    """ Determine the chopper timing and return a data object with Time and Sample coordinates

        options:
        'Phase' : Phase/frame to process in the chopper period (0, ...)
        'State'  :  Dictionary with {'Chop': 0 or 1, 'Defl': 0 or 1}
        'Start delay' : Start delay relative to the phase start in microsec or None. If None use seom standard value depending on settings.
        'End delay' : Delay at the end of the phase in microsec. Positive means later times. If None use seom standard value depending on settings.
    """
    try:
        sch = config['Chopper scheme']
    except KeyError:
        raise("No chopper information in configuration file.")

    state_set = False
    if (options['State'] is not None):
        state = options['State']
        state_set = True
    else:
        if (options['Phase'] is not None):
            phase = options['Phase']
            if (phase is not ...):
                phase = int(phase)
        else:
            raise ValueError("Either State or Phase should be set to determine chopper timing.")
    try:
        length, chop, defl = process_chopper_setup(config)
    except Exception as e:
        raise e

    if (state_set):
        for i in range(len(chop)):
            try:
                if ((chop[i] == state['Chop']) and (defl[i] == state['Defl'])):
                    phase = i
                    break
            except KeyError:
                raise ValueError("Invalid chopper state description. Valid: {'Chop': 0 or 1, 'Defl': 0 or 1}")
        try:
            phase
        except NameError:
            raise ValueError("No such chopper state.")
    else:
        if ((phase != ...) and (phase >= length(chop))):
            raise ValueError("Invalid chopper phase number.")

    if (phase is ...):
        start_phase = 0
        end_phase = len(chop) - 1
    else:
        start_phase = phase
        end_phase = phase
    chop_clk_in_sample = config['APDCAM_f_sample'] / config['Chopper clock']
    if (chop_clk_in_sample >= 1):
        if (abs(chop_clk_in_sample - round(chop_clk_in_sample)) > 1e-6):
            print("Chopper clock period time is not compatible with sample period.")
        else:
            chop_clk_in_sample = round(chop_clk_in_sample)
    else:
        if (abs(Decimal(1)/chop_clk_in_sample - round(Decimal(1)/chop_clk_in_sample)) > 1e-6):
            print("Chopper clock period time is not compatible with sample period (sampletime longer).")            


    # These are values [microsec] determined from the signals
    switch_time = 7/Decimal(2000000)
    switch_time_sample = round(switch_time / config['APDCAM_sampletime'])
    if (config['APDCAM_f_ADC'] == Decimal(20e6)):
        if (config['APDCAM_f_sample'] == Decimal(2e6)):
            instrument_delay = -9/Decimal(1000000)
            # instrument_delay = -14/Decimal(1000000)
            instrument_delay = -28/Decimal(1000000)
            # instrument_delay = -30/Decimal(1000000)
            # instrument_delazy= -17/3-1/3*period time[microsec]
        elif (config['APDCAM_f_sample'] == Decimal(1e6)):
            instrument_delay = -6/Decimal(1000000)
        else:
            warnings.warn("Unknown sample delay relative to chopper. Check chopper timing.")
            instrument_delay = -9/Decimal(1000000)
    else:
        warnings.warn("Unknown sample delay relative to chopper. Check chopper timing.")
        instrument_delay = -9/Decimal(1000000)
       
        
    clock_unit = Decimal(1.) / config['Chopper clock']
    if (config['Chopper mode'] == 'Camera'):
        period_time_clk = round(config['CMOS frametime'] / Decimal(1000) \
                                * config['Chopper clock']) * len(chop)
        # Calculating the chopper start and stop times in chopper clock units relative to the the APDCAM start time
        if (chop[start_phase] == 1):
            start_clk = round(( config['CMOS frametime'] / Decimal(1000) * start_phase
                              -(config['CMOS frametime'] - config['CMOS exptime']) / 2 / Decimal(1000) \
                              - 2 / Decimal(100000) ) / clock_unit \
                              )
            stop_clk =  round((config['CMOS frametime'] / Decimal(1000) * (end_phase + 1)\
                              -(config['CMOS frametime'] - config['CMOS exptime']) / 2 / Decimal(1000)  \
                              + 2 / Decimal(100000) ) / clock_unit \
                              )
        else:
            start_clk = round((config['CMOS frametime'] / Decimal(1000) * start_phase \
                              -(config['CMOS frametime'] - config['CMOS exptime']) / 2 / Decimal(1000)  \
                              + 2 / Decimal(100000) ) / clock_unit  \
                              )
            stop_clk =  round((config['CMOS frametime'] / Decimal(1000) * (end_phase + 1) \
                              -(config['CMOS frametime'] - config['CMOS exptime']) / 2 /  Decimal(1000)  \
                              - 2 / Decimal(100000)) / clock_unit  \
                              )
    else:
        period_time_clk = int(round(config['Chopper period'] / clock_unit))
        st = 0
        for i in range(start_phase):
            st += length[i]
        start_clk = int(round(st / 100 * period_time_clk))
        st = 0
        for i in range(end_phase + 1):
            st += length[i]
        stop_clk = int(round(st / 100 * period_time_clk))

    if (options['Start delay'] is None):
        if (config['APDCAM_clock_source'] == 'External'):
            if (config['Chopper period'] > 0.01):
                start_delay_sample = int(round(1e-3 / float(config['APDCAM_sampletime'])))
                print("APDCAM clock source is external, chopper drifts relative to signal. Shifting chopper start by 1 ms.")
            else:
                raise ValueError("APDCAM clock source is external, chopper drifts relative to signal. Set 'start delay'.")
        else:
            start_delay_sample = 0
    else:
        try:
            start_delay = float(options['Start delay'])
            start_delay_sample = int(round(start_delay *1e-6 / float(config['APDCAM_sampletime'])))            
        except ValueError:
            raise ValueError("Invalid chopper start delay.")

    if (options['End delay'] is None):
        if (config['APDCAM_clock_source'] == 'External'):
            if (config['Chopper period'] > 0.01):
                stop_delay_sample = int(round(-1e-3 / float(config['APDCAM_sampletime'])))
                print("APDCAM clock source is external, chopper drifts relative to signal. Shifting chopper end forward by 1 ms.")
            else:
                raise ValueError("APDCAM clock source is external, chopper drifts relative to signal. Set 'stop delay'.")
        else:
            stop_delay_sample = 0
    else:
        try:
            stop_delay = float(options['End delay'])
            stop_delay_sample = int(round(stop_delay * 1e-6 / float(config['APDCAM_sampletime'])))
        except ValueError:
            raise ValueError("Invalid chopper end delay.")

    if (read_samplerange is not None):
        start_sample = round(start_clk * chop_clk_in_sample)
        period_time_sample = round(period_time_clk * chop_clk_in_sample)
        if (period_time_sample < 2):
            raise ValueError("Less than one sample in chopper period.")
        if (start_sample < read_samplerange[0]):
            start_ind = int((read_samplerange[0] - start_sample)
                            / period_time_sample + 1)
            start_clk += period_time_clk * start_ind
            stop_clk += period_time_clk * start_ind
        stop_sample = stop_clk * chop_clk_in_sample
        n_period =  int(float((read_samplerange[1] - stop_sample) / period_time_sample) + 0.5)
        if (n_period <= 0):
            raise ValueError("No chopper intervals in time (sample) range.")
    else:
        n_period = int(round(config['APDCAM_samplenumber'] /  (period_time_clk * chop_clk_in_sample)))

    mode = flap.CoordinateMode(equidistant=True, range_symmetric=False)
    c_sample = copy.deepcopy(flap.Coordinate(name='Sample',
                                             unit='n.a.',
                                             mode=mode,
                                             start=round(start_clk * chop_clk_in_sample
                                                    + instrument_delay/config['APDCAM_sampletime'])
                                                    + start_delay_sample + switch_time_sample,
                                             step=round(period_time_clk*chop_clk_in_sample) ,
                                             value_ranges=[0,round((stop_clk-start_clk) * chop_clk_in_sample)
                                                           + stop_delay_sample - start_delay_sample - switch_time_sample],
                                           dimension_list=[0]
                                           ))
    c_time = copy.deepcopy(flap.Coordinate(name='Time',
                                           unit='Second',
                                           mode=mode,
                                           start = c_sample.start * float(config['APDCAM_sampletime'])+ float(config['APDCAM_starttime']),
                                           step=period_time_clk*clock_unit,
                                           value_ranges=[0, c_sample.value_ranges[1] * float(config['APDCAM_sampletime'])],
                                           dimension_list=[0]
                                           ))


    d = copy.deepcopy(flap.DataObject(data_shape=[n_period],
                                      data_unit=flap.Unit(name='Interval', unit='n.a.'),
                                      coordinates=[c_sample, c_time],
                                      data_source = 'W7X_ABES',
                                      exp_id = config['ShotID'],
                                      info=config
                                      )
                          )
    d.check()


    return d


def w7x_abes_get_data(exp_id=None, data_name=None, no_data=False, options=None, coordinates=None, data_source=None):
    """ Data read function for the W7-X Alkali BES diagnostic
    data_name: ABES-xxx, Rx, Lx, ... (string) depending on configuration file
               Chopper_times : To read the chopper state (will return intervals in Time and Sample ccoordinate)

    exp_id: Experiment ID, YYYYMMDD.xxx
    Unix style regular expressions are allowed with * and []
                       Can also be a list of data names, eg. ['ABES-1','ABES-4']
    coordinates: List of flap.Coordinate() or a single flap.Coordinate
                 Defines read ranges. The following coordinates are interpreted:
                     'Sample': The read samples
                     'Time': The read times
                     Only a single equidistant range is interpreted in c_range.
    options:
        'Scaling':  string
           'Digit'
           'Volt'
        'Offset timerange': list
           Time range for offset subtraction
        'Datapath': string
           Data path (string)
        'Amplitude calibration': bool
           True/False do/don't do amplitude calibration of the data
        'Amplitude calib. path': string
           Calibration directory name
        'Amplitude calib. file': string
           Calibration cld file name.
        'Partial intervals' : bool
           True: Keep partial chopper intervals extending through the start/end of timerange
           False: Keep only full chopper intervals
        'Resample' : float
           Resample data to this sample frequency
        'Spatial calibration': bool
           True: Attempt spatial calibration, add device coordinates
            
        For further options see Chopper_times see chopper_timing_data()

    """
    if (exp_id is None):
        raise ValueError('exp_id should be set for W7X ABES.')
    default_options = {'Datapath': 'data',
                       'Scaling':'Digit',
                       'Offset timerange': None,
                       'Amplitude calibration': False,
                       'Amplitude calib. path': 'cal',
                       'Amplitude calib. file': None,
                       'Phase' : None,
                       'State' : None,
                       'Start delay': None,
                       'End delay': None,
                       'Spatial calibration': False,
                       'Spatial calib. path': 'spatcal',
                       'Partial intervals': False,
                       'Resample' : None
                       }
    _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES')
    try:
        datapath_base = _options['Datapath']
    except (KeyError, TypeError):
        datapath_base = 'data'
    if (type(exp_id) is not str):
        raise ValueError("exp_id should be a string of format yyyymmdd.xxx")
    datapath = os.path.join(datapath_base,exp_id)
    xmlfile = os.path.join(datapath, exp_id + '_config.xml')
    xml = flap.FlapXml()
    if os.path.exists(xmlfile) is False:
        raise ValueError('XML file '+xmlfile+' does not exist')
    try:
        while os.access(xmlfile, os.R_OK) is False:
            time.sleep(0.1)
        xml.read_file(xmlfile)
    except Exception:
        raise IOError("Error reading XML file:" + xmlfile)
    try:
        if (xml.head.tag != 'ShotSettings'):
            raise ValueError("No ShotSettings entry found in XML file " + xmlfile + ".")
        if (xml.head.attrib['Experiment'] != "W7-X A-BES"):
            raise ValueError(xmlfile + " is not a W7-X ABES config file.")
    except Exception:
        raise ValueError("File format error in " + xmlfile)
    try:
        config = abes_get_config(xml)
    except Exception as e:
        raise e

    # Ensuring that the data name is a list
    if type(data_name) is not list:
        chspec = [data_name]
    else:
        chspec = data_name

    # Finding read_range (timerange) and read_samplerange

    read_range = None
    read_samplerange = None
    if (coordinates is not None):
        if (type(coordinates) is not list):
             _coordinates = [coordinates]
        else:
            _coordinates = coordinates
        for coord in _coordinates:
            if (type(coord) is not flap.Coordinate):
                raise TypeError("Coordinate description should be flap.Coordinate.")
            if ((coord is None) or (coord.c_range is None)):
                continue
            if (coord.unit.name == 'Time'):
                if (coord.mode.equidistant):
                    read_range = [float(coord.c_range[0]),float(coord.c_range[1])]
                    if (read_range[1] <= read_range[0]):
                        raise ValueError("Invalid read timerange.")
                else:
                    raise NotImplementedError("Non-equidistant Time axis is not implemented yet.")
                break
            if coord.unit.name == 'Sample':
                if (coord.mode.equidistant):
                    read_samplerange = coord.c_range
                    if (read_samplerange[1] <= read_samplerange[0]):
                        raise ValueError("Invalid read samplerange.")

                else:
                    raise \
                        NotImplementedError("Non-equidistant Sample axis is not implemented yet.")
                break
    if ((read_range is None) and (read_samplerange is None)):
        read_samplerange = np.array([0,config['APDCAM_samplenumber']])
    if (read_range is not None):
        read_range = np.array(read_range)
    if (read_samplerange is None):
        read_samplerange = np.rint((read_range - float(config['APDCAM_starttime']))
                                   / float(config['APDCAM_sampletime'])
                                   )
    else:
        read_samplerange = np.array(read_samplerange)
    if ((read_samplerange[1] < 0) or (read_samplerange[0] >= config['APDCAM_samplenumber'])):
        raise ValueError("No data in time range.")
    if (read_samplerange[0] < 0):
        read_samplerange[0] = 0
    if (read_samplerange[1] >= config['APDCAM_samplenumber']):
        read_samplerange[1] = config['APDCAM_samplenumber']
    read_range = float(config['APDCAM_starttime']) \
                       + read_samplerange * float(config['APDCAM_sampletime'])
    if (_options['Resample'] is not None):
        if (_options['Resample'] > 1 / config['APDCAM_sampletime']):
            raise ValueError("Resampling frequency should be below the original sample frequency.")
        resample_binsize = int(round((1 / _options['Resample']) / float(config['APDCAM_sampletime'])))
    try:
        ch_chop = chspec.index('Chopper_time')
        chopper_signal = True
    except ValueError:
        chopper_signal = False
    if (chopper_signal):
        if (len(chspec) != 1):
            raise ValueError("Chopper_time should be read separately from data.")
        # Reading chopper timing
        try:
            chop = chopper_timing_data_object(config, _options, read_samplerange=read_samplerange)
        except Exception as e:
            raise e
        return chop

    # Finding the desired channels
    signal_list = config['signal_list']
    ADC_list = config['ADC_list']
    try:
        signal_proc, signal_index = flap.select_signals(signal_list,chspec)
    except ValueError as e:
        raise e
    ADC_proc = []
    for i in signal_index:
        ADC_proc.append(ADC_list[i])


    scale_to_volts = False
    dtype = np.int16
    data_unit = flap.Unit(name='Signal',unit='Digit')
    if _options is not None:
        try:
            if (_options['Scaling'] == 'Volt'):
                scale_to_volts = True
                dtype = float
                data_unit = flap.Unit(name='Signal',unit='Volt')
        except (NameError, KeyError):
            pass

    try:
        offset_timerange = _options['Offset timerange']
    except (NameError, KeyError):
        offset_timerange = None

    try:
        if (offset_timerange is not None):
            if (type(offset_timerange) is not list):
                raise ValueError("Invalid Offset timerange. Should be list or string.")
            if ((len(offset_timerange) != 2) or (offset_timerange[0] >= offset_timerange[1])) :
                raise ValueError("Invalid Offset timerange.")
            offset_samplerange = np.rint((np.array(offset_timerange) - float(config['APDCAM_starttime']))
                                       / float(config['APDCAM_sampletime']))
            if ((offset_samplerange[0] < 0) or (offset_samplerange[1] >= config['APDCAM_samplenumber'])):
                raise ValueError("Offset timerange is out of measurement time.")
            offset_data = np.empty(len(ADC_proc), dtype='int16')
            for i_ch in range(len(ADC_proc)):
                fn = os.path.join(datapath, "Channel_{:03d}.dat".format(ADC_proc[i_ch] - 1))
                with open(fn,"rb") as f:
                    try:
                        f.seek(int(offset_samplerange[0]) * 2, os.SEEK_SET)
                        d = np.fromfile(f, dtype=np.int16, count=int(offset_samplerange[1]-offset_samplerange[0])+1)
                    except Exception:
                        raise IOError("Error reading from file: " + fn)
                offset_data[i_ch] = np.int16(np.mean(d))
            if (scale_to_volts):
                offset_data = ((2 ** config['APDCAM_bits'] - 1) - offset_data) \
                            / (2. ** config['APDCAM_bits'] - 1) * 2
            else:
                offset_data = (2 ** config['APDCAM_bits'] - 1) - offset_data
    except ValueError as e:
        warnings.warn(str(e))
        offset_timerange = None

    ndata_read = int(read_samplerange[1] - read_samplerange[0] + 1)
    if (_options['Resample'] is not None):
        ndata_out = int(ndata_read / resample_binsize) 
        ndata_read = ndata_out * resample_binsize
    else:
        ndata_out = ndata_read

    if (no_data is False):
        if ndata_out * len(ADC_proc) * 32 > psutil.virtual_memory().available:
            raise MemoryError("Not enough memory for reading data")

        if (len(ADC_proc) != 1):
            if (_options['Resample'] is not None):
                data_arr = np.empty((ndata_out, len(ADC_proc)), dtype=float)
                data_err = np.empty((ndata_out, len(ADC_proc)), dtype=float)
            else:
                data_arr = np.empty((ndata_out, len(ADC_proc)), dtype=dtype)
                data_err = None
        for i in range(len(ADC_proc)):
            fn = os.path.join(datapath, "Channel_{:03d}.dat".format(ADC_proc[i] - 1))
            with open(fn,"rb") as f:
                try:
                    f.seek(int(read_samplerange[0]) * 2, os.SEEK_SET)
                except Exception:
                    raise IOError("Error reading from file: " + fn)
    
                try:
                    d = np.fromfile(f, dtype=np.int16, count=ndata_read)
                except Exception:
                    raise IOError("Error reading from file: " + fn)
                if (len(d) != ndata_read):
                    raise(IOError("Truncated file, could not read enough data."))
                if (scale_to_volts):
                    d = ((2 ** config['APDCAM_bits'] - 1) - d) \
                                / (2. ** config['APDCAM_bits'] - 1) * 2
                else:
                    d = (2 ** config['APDCAM_bits'] - 1) - d
                if (offset_timerange is not None):
                        d -= offset_data[i]
            if (_options['Resample'] is not None):
                d = d.astype(float)
                d_resample = np.zeros(ndata_out,dtype=float)
                d_error = np.zeros(ndata_out,dtype=float)
                if (ndata_out > resample_binsize):
                    for i_slice in range(0,resample_binsize):
                        d_resample += d[slice(i_slice,len(d),resample_binsize)]
                        d_error += d[slice(i_slice,len(d),resample_binsize)] ** 2
                else:
                    for i_resamp in range(0,len(d_resample)):
                        d_resample[i_resamp] = np.sum(d[i_resamp * resample_binsize : (i_resamp + 1) * resample_binsize])
                        d_error[i_resamp] = np.sum(d[i_resamp * resample_binsize : (i_resamp + 1) * resample_binsize] ** 2)
                d_error = np.sqrt(d_error / resample_binsize - (d_resample / resample_binsize) ** 2)
                d = d_resample / resample_binsize
                if (len(ADC_proc) == 1):
                    data_err = d_error
                else:
                    data_err[:,i] = d_error
            if (len(ADC_proc) == 1):
                data_arr = d
            else:
                data_arr[:,i] = d    
        try:
            data_arr, data_err, calfac_err = calibrate(data_arr, signal_proc, read_range, exp_id=exp_id, options=_options)
        except Exception as e:
            raise e
        data_dim = data_arr.ndim    
    else:
        if (len(ADC_proc) != 1):
            data_dim = 2
        else:
            data_dim = 1

    coord = [None]*data_dim*2
    c_mode = flap.CoordinateMode(equidistant=True)
    if (read_range is None):
        read_range = float(config['APDCAM_sampletime']) + read_samplerange * float(config['APDCAM_sampletime'])
    if (_options['Resample'] is not None):
        tstart = read_range[0] + float(config['APDCAM_sampletime']) * resample_binsize / 2
        tstep = float(config['APDCAM_sampletime']) * resample_binsize
    else:
        tstart = read_range[0]
        tstep = float(config['APDCAM_sampletime'])
    coord[0] = copy.deepcopy(flap.Coordinate(name='Time',
                                             unit='Second',
                                             mode=c_mode,
                                             start=tstart,
                                             step=tstep,
                                             dimension_list=[0])
                             )

    if (_options['Resample'] is not None):
        s_start = read_samplerange[0] + resample_binsize / 2
        s_step = resample_binsize
    else:
        s_start = read_samplerange[0]
        s_step = 1
    coord[1] = copy.deepcopy(flap.Coordinate(name='Sample',
                                             unit='n.a.',
                                             mode=c_mode,
                                             start=s_start,
                                             step=s_step,
                                             dimension_list=[0])
                             )
    if (data_dim > 1):
        ch_proc = []
        for ch in signal_proc:
            if (ch[0:4] != 'ABES'):
                ch_proc = []
                break
            ch_proc.append(int(ch[5:]))
        c_mode = flap.CoordinateMode(equidistant=False)
        if (ch_proc != []):
            coord[2] = copy.deepcopy(flap.Coordinate(name='Channel',
                                                     unit='n.a.',
                                                     mode=c_mode,
                                                     shape=len(ch_proc),
                                                     values=ch_proc,
                                                     dimension_list=[1])
                                 )
        else:
            coord = coord[:2]+[coord[3]]
        coord[-1] = copy.deepcopy(flap.Coordinate(name='Signal name',
                                                 unit='n.a.',
                                                 mode=c_mode,
                                                 shape=len(signal_proc),
                                                 values=signal_proc,
                                                 dimension_list=[1])
                                 )

    data_title = "W7-X ABES data"
    if (data_arr.ndim == 1):
        data_title += " (" + signal_proc[0] + ")"

    d = flap.DataObject(data_array=data_arr,
                        error = data_err,
                        data_unit=data_unit,
                        coordinates=coord,
                        exp_id=exp_id,
                        data_title=data_title,
                        info={'Options':_options, 'Calibration factor error': calfac_err,'Config':config},
                        data_source="W7X_ABES")
    if _options['Spatial calibration'] is True:
        # Getting the spatial calibration
        d = add_coordinate(d, ['Device R', 'Beam axis'],
                           options={"Shot spatcal dir": flap.config.get("Module W7X_ABES","Spatial calibration directory")})
    return d


def add_coordinate(data_object,
                   coordinates,
                   exp_id=None,
                   options=None):
    '''
    Routine for adding spatial data to W7-X ABES measurements
    data_object - the object to which the coordinate should be added
    coordinates - a list of coordinate names to be added
                 available: "Device x", "Device y", "Device z", "Device R", "Device Z", "Beam axis"
    options - a dictionary of options
              available: 'spatcal_dir' - the location of calibration data
    '''
    default_options = {'Shot spatcal dir': os.path.join(os.path.dirname(os.path.abspath(__file__)),'spatcal'),
                       "Channels":''
                       }
    _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES')

    if exp_id is None:
        exp_spatcal = spatcal.ShotSpatCal(data_object.exp_id, options=_options)
    else:
        exp_spatcal = spatcal.ShotSpatCal(exp_id, options=_options)
    try:
        exp_spatcal.read(options=_options)
    except Exception as e:
        print(e)

    # getting the dimension of the channel coordinate, this should be the same as the spatial coordinate
    data_coord_list = np.array([coord.unit.name for coord in data_object.coordinates])
    if 'Signal name' in data_coord_list:
        channel_coordinate_dim = data_object.get_coordinate_object('Signal name').dimension_list[0]
        channel_names = data_object.get_coordinate_object('Signal name').values
    else:
        channel_coordinate_dim = data_object.get_coordinate_object('Channel').dimension_list[0]
        channel_names = data_object.get_coordinate_object('Channel').values
    for coord_name in coordinates:
        coord_object = exp_spatcal.create_coordinate_object(channel_coordinate_dim, coord_name,
                                                         channel_names=channel_names)
        data_object.add_coordinate_object(coord_object)

    return data_object

def regenerate_time_sample(d):
    """ Regenerate Time and Sample coordinate after chopper slicing aon Sample coordinate
        and averaging for chopper interval.
    """
    try:
        # Trying to get Time coordinate. If not present regenerating it
        d.get_coordinate_object('Time')
    except ValueError:
        ct = d.get_coordinate_object('Start Time in int(Sample)')
        c_shift = d.get_coordinate_object('Rel. Time in int(Sample)')
        if (c_shift.dimension_list != []):
            raise ValueError("Rel. Time in int(Sample) is not constant.")
        if (not ct.mode.equidistant):
            if c_shift.values == None:
                c_shift.values = c_shift.start
            try:
                ct.values += c_shift.values[0]
            except IndexError:
                ct.values += c_shift.values
            #check if new coordinate is equidistant
            if len(ct.dimension_list) == 1:
                steps = ct.values[1:]-ct.values[:-1]
                accuracy = np.max(steps)/np.min(steps)
                if accuracy-1 < 1e-10:
                    ct.start = ct.values[0]
                    ct.step = np.mean(steps)
                    ct.mode.equidistant = True
        else:
            try:
                ct.start += c_shift.values[0]
            except IndexError:
                ct.start += c_shift.values
            except TypeError:
                pass
        ct.unit.name='Time'
        
        d.del_coordinate('Rel. Time in int(Sample)')
    try:
        # Trying to get Sample coordinate. If not present regenerating it
        d.get_coordinate_object('Sample')
    except ValueError:
        ct = d.get_coordinate_object('Start Sample in int(Sample)')
        c_shift = d.get_coordinate_object('Rel. Sample in int(Sample)')
        if (c_shift.dimension_list != []):
            raise ValueError("Rel Sample in int(Sample) is not constant.")
        if (not ct.mode.equidistant):
            if c_shift.values == None:
                c_shift.values = c_shift.start
            try:
                ct.values += c_shift.values[0]
            except IndexError:
                ct.values += c_shift.values
            #check if new coordinate is equidistant
            if len(ct.dimension_list) == 1:
                steps = ct.values[1:]-ct.values[:-1]
                accuracy = np.max(steps)/np.min(steps)
                if accuracy-1 < 1e-10:
                    ct.start = ct.values[0]
                    ct.step = np.mean(steps)
                    ct.mode.equidistant = True
        else:
            try:
                ct.start += c_shift.values[0]
            except IndexError:
                ct.start += c_shift.values
            except TypeError:
                pass
        ct.unit.name='Sample'
        d.del_coordinate('Rel. Sample in int(Sample)')
    try:
        d.del_coordinate('Interval(Sample)')
        d.del_coordinate('Interval(Sample) sample index')
    except ValueError:
        pass

def proc_chopsignals_single(dataobject=None, exp_id=None,timerange=None,
                            samplerange=None, signals='ABES-1',
                            on_options=None, off_options=None, test=None, options={}):
    """ Calculate signals in beam on and beam/off phases of the measurement and
        correct the beam-on phases with the beam-off. The result is "ABES" and "ABES_back" data object
        in the FLAP storage.
        INPUT:
            There are two modes of this function:
                The first one obtains the data directly from a measurement file
                    exp_id: exp_id (no default)
                    signals: List of measurement signals. Default is ABES-[1-40]
                The second one takes a dataobject as an input
                    dataobject: the input dataobject file
            timerange: Time range to process. Default is all times.
            on_options: Options for the  for the get_data function when beam_on is read
            off_options: Options for the get_data function when beam_off is read
            test: Plot test plots if True
            options:
                Average Chopping Period: Whether the data should be averaged for a single chopper time period. Per 
                                         default this is True. If this option is False then all measurement time points
                                         of the beam on state is saved. This could be useful if slow beam chopping is
                                         used and the background signal is reasonably constant relative to the analyzed
                                         process.
        OUTPUT: The background subtracted A-BES data
    """

    options_default = {'Average Chopping Period': True,
                       'Off-axis Correction':False,
                       'Deflection': 0}
    options = {**options_default, **options}

    # Obtaining the chopper data
    if dataobject is not None:
        exp_id = dataobject.exp_id

    o = copy.deepcopy(on_options)
    if 'Datapath' in options.keys():
        o['Datapath'] = options['Datapath']
    if 'W7X_ABES' not in flap.list_data_sources():
        register()
    o.update({'State':{'Chop': 0, 'Defl': options['Deflection']}})

    if timerange is None and samplerange is None:
        if "Sample" not in dataobject.coordinate_names():
            timerange = [np.min(dataobject.coordinate('Time')[0]), np.max(dataobject.coordinate('Time')[0])]
            coord_dict = {'Time':timerange}
        else:
            samplerange = [np.min(dataobject.coordinate('Sample')[0]), np.max(dataobject.coordinate('Sample')[0])]
            coord_dict = {'Sample':samplerange}
    elif timerange is not None:
        coord_dict = {'Time':timerange}
    else:
        coord_dict = {'Sample':samplerange}

    d_beam_on=flap.get_data('W7X_ABES',
                            exp_id=exp_id,
                            name='Chopper_time',
                            coordinates=coord_dict,
                            options=o,\
                            object_name='Beam_on',
                            )
    
    o = copy.deepcopy(off_options)
    if 'Datapath' in options.keys():
        o['Datapath'] = options['Datapath']
    o.update({'State':{'Chop': 1, 'Defl': 0}})
    d_beam_off=flap.get_data('W7X_ABES',
                            exp_id=exp_id,
                            name='Chopper_time',
                            coordinates=coord_dict,
                            options=o,\
                            object_name='Beam_off',
                                )
        
    if (test):
        from matplotlib import pyplot as plt
        import time
        plt.close('all')
        if dataobject is None:
            flap.plot('ABES', axes='Time', plot_options={'marker': 'o'})
        else:
            dataobject.plot(axes='Time', options={'All':True}, plot_options={'marker': 'o'})
#        flap.plot('ABES',axes='Time',plot_type='scatter')
        d_beam_on.plot(plot_type='scatter', axes=['Time', plt.ylim()[1]], options={'Force': True,'All': True})
        d_beam_off.plot(plot_type='scatter', axes=['Time', plt.ylim()[0]], options={'Force': True,'All': True})
        plt.savefig(str(time.time())+'.png')

    # Background subtraction
    if dataobject is None:
        # in this case the flap storage is used for obtaining the data by experiment ID
        flap.get_data('W7X_ABES',
                      exp_id=exp_id,
                      coordinates={'Time':timerange},
                      name=signals,
                      object_name='ABES'
                      )
        if (test):
            plt.close('all')
            flap.plot('ABES', axes='Time', plot_options={'marker': 'o'})
            d_beam_on.plot(plot_type='scatter', axes=['Time', 2], options={'Force': True,'All': True})
            d_beam_off.plot(plot_type='scatter', axes=['Time', 0.1], options={'Force': True,'All': True})
        d = flap.slice_data('ABES',slicing={'Sample':d_beam_on})
        
        if options['Average Chopping Period'] is True:
            d = d.slice_data(summing={'Rel. Sample in int(Sample)':'Mean'})
            regenerate_time_sample(d)
        else:
            add_absolute_time_sample(d)

        flap.add_data_object(d,'ABES_on')
        
        d = flap.slice_data('ABES',slicing={'Sample':d_beam_off})
        d = d.slice_data(summing={'Rel. Sample in int(Sample)':'Mean'})
        regenerate_time_sample(d)
        
        if options["Off-axis correction"] is True:
            raise NotImplementedError("Off-axis correction not implemented if input is not a dataobject")
        
        flap.add_data_object(d,'ABES_off')
        flap.slice_data('ABES_off',slicing={'Time':flap.get_data_object('ABES_on')},options={'Inter':'Linear'},output_name='ABES_back')
        # Ensuring that only those samples are kept which also have a background
        #    flap.slice_data('ABES_on',slicing={'Start Sample in int(Sample)':flap.get_data_object('ABES_off_resampled')},options={'Inter':'Linear'},output_name='ABES_on')
        if (test):
            plt.figure()
            flap.plot('ABES')
            flap.plot('ABES_on',plot_type='scatter')
            flap.plot('ABES_on')
            flap.plot('ABES_off',plot_type='scatter')
            flap.plot('ABES_off')
            flap.plot('ABES_back',plot_type='scatter')
         
        d=flap.get_data_object('ABES_on')
        d_back = flap.get_data_object('ABES_back')
        d.data -= d_back.data
        
        
        #add electrc noise using offset data
        del o['State']
        if 'Signal name' in dataobject.coordinate_names():
            try:
                offset = flap.get_data('W7X_ABES',
                              exp_id=exp_id,
                              coordinates={'Time':[-0.1,0.1]},
                              name=dataobject.get_coordinate_object('Signal name').values[0],
                              object_name='ABES'
                              )
            except OSError:
                #no amplitude calibration present yet
                offset = flap.get_data('W7X_ABES',
                              exp_id=exp_id,
                              coordinates={'Time':[-0.1,0.1]},
                              name=dataobject.get_coordinate_object('Signal name').values[0],
                              object_name='ABES',
                              options = {'Amplitude calibration': False}
                              )
            chop = np.mean(d_beam_on.coordinate("Time")[2]-d_beam_on.coordinate("Time")[1])
            if chop < 1e-4:
                if np.min(offset.coordinate("Time")[0]) < -1e-6 and options['Average Chopping Period'] is True:
                    offset = offset.slice_data(slicing={"Time": flap.Intervals(-0.1,-1e-6)})
                    electric_noise = np.var(offset.data, axis=0)
                    d.error = np.sqrt(d.error**2+electric_noise)
                elif options['Average Chopping Period'] is True:
                    standard_error = np.array([8.52008759e-03, 6.05577881e-03, 5.61784578e-03, 4.01606559e-03,
                           4.80343614e-03, 4.96760886e-03, 4.31857827e-03, 4.61216619e-03,
                           5.07540636e+00, 5.77341359e+00, 5.17886695e+00, 6.67524467e+00,
                           4.10515027e+00, 3.20796480e+00, 4.06675798e+00, 4.07734300e+00,
                           4.09595980e+00, 3.70180586e+00, 4.78968865e+00, 3.76662174e+00,
                           4.02454327e+00, 3.90762905e+00, 5.00180528e+00, 3.25465952e+00,
                           4.32225647e+00, 3.23294377e+00, 5.75256366e+00, 4.56699600e+00,
                           5.02907424e+00, 4.54773953e+00, 6.10916503e+00, 4.75039061e+00,
                           3.76776202e+00, 5.17901390e+00, 4.44659086e+00, 6.81016910e+00,
                           7.95957335e+00, 1.45460916e+01, 1.10220756e+01, 1.85718833e+01])
                    d.error = np.sqrt(d.error**2+standard_error)

#        flap.add_data_object(d,'ABES')

#        # error approximation
#        d_beam_off = flap.get_data_object('ABES_off')
##        beam_off_data = d_beam_off.slice_data(summing={'Rel. Sample in int(Sample)':'Mean'})
##        regenerate_time_sample(beam_off_data)
#        beam_off_data = beam_off_data.slice_data(slicing={'Time':d_beam_off},options={'Inter':'Linear'})
#        background_error = np.average((d_beam_off.data-beam_off_data.data.reshape(np.shape(d_beam_off.data)))**2, axis=0)\
#                           *len(beam_off_data.data)/(len(beam_off_data.data)-1)


        flap.delete_data_object(['ABES_on','ABES_off','Beam_on','Beam_off'],exp_id=exp_id)
        try:
            flap.delete_data_object('ABES')
        except:
            pass
        gc.collect()
        if (test):
            plt.figure()
            d.plot(axes='Time')
            
        return d
    else:

        # in this case the passed dataobject is used and only the chopper data is obtained from file
        #there is sometimes some weird slicing error, probably  due to rounding, that is to be corrected in the following
        while np.min(d_beam_on.coordinate("Sample")[0]) < np.min(dataobject.coordinate("Sample")[0]):
                d_beam_on.get_coordinate_object("Sample").start += d_beam_on.get_coordinate_object("Sample").step[0]
        dataobject_beam_on = dataobject.slice_data(slicing={'Sample': d_beam_on})

        # For dataobject_beam_on.data the 0 dimension is along a constant 'Start Time in int(Time)' and 
        # "Rel. Time in int(Time)" varies. For dimension 1 it is 'Start Time in int(Time)' that varies
        dataobject_beam_on = process_chopped_dataobject(dataobject_beam_on, options=options)

        if d_beam_off.coordinate('Sample')[0][-1] > dataobject.coordinate('Sample')[0][-1]:
            d_beam_off.shape = [int((dataobject.coordinate('Sample')[0][-1]-d_beam_off.coordinate('Sample')[0][0])/
                                   d_beam_off.get_coordinate_object('Sample').step)]

        dataobject_beam_off = dataobject.slice_data(slicing={'Sample': d_beam_off}, options={'Partial intervals':True})

        dataobject_beam_off.error = None
        dataobject_beam_off = process_chopped_dataobject(dataobject_beam_off, options={'Average Chopping Period': True})

        dataobject_background = dataobject_beam_off.slice_data(slicing={'Time': dataobject_beam_on},
                                                               options={'Inter': 'Linear'})
        if test is True:
            from matplotlib import pyplot as plt
            import time
            plt.figure()
            plt.scatter(dataobject_beam_on.get_coordinate_object("Time").values, dataobject_beam_on.data, color ='green')
            plt.scatter(dataobject_beam_off.get_coordinate_object("Time").values, dataobject_beam_off.data, color='red')
            plt.savefig(str(time.time())+'.png')
#            plt.plot(dataobject.coordinate("Time")[0], dataobject.data, color='blue')

        dataobject_beam_on.data -= dataobject_background.data.reshape(np.shape(dataobject_beam_on.data))

        # calculating the error for the beam off part
        dataobject_beam_off_error = copy.deepcopy(dataobject_beam_off)
        dataobject_beam_off_error.data = dataobject_beam_off.error
        background_error = dataobject_beam_off_error.slice_data(slicing={'Time': dataobject_beam_on}, options={'Inter': 'Linear'})
        background_error = background_error.data.reshape(np.shape(dataobject_beam_on.data))
        dataobject_beam_on.error = np.asarray(np.sqrt(dataobject_beam_on.error**2 + background_error**2))

        #adding the calibration factor error
        if 'Calibration factor error' in dataobject_beam_on.info:
            calfac_error = dataobject_beam_on.info['Calibration factor error']
            if 'Signal name' in dataobject.coordinate_names() and calfac_error is not None:
                calfac_curr = calfac_error.slice_data(slicing =
                                                      {'Signal name': dataobject.get_coordinate_object('Signal name').values[0]})
                calfac_curr_errors = dataobject_beam_on.data * calfac_curr.data
                dataobject_beam_on.error = np.asarray(np.sqrt(dataobject_beam_on.error**2 + calfac_curr_errors**2))
            else:
                print('Calibration error not considered')

        #add electrc noise using offset data
        del o['State']
        if 'Signal name' in dataobject.coordinate_names() and options['Average Chopping Period'] is True:
            try:
                offset = flap.get_data('W7X_ABES',
                              exp_id=exp_id,
                              coordinates={'Time':[-0.1,0.1]},
                              name=dataobject.get_coordinate_object('Signal name').values[0],
                              object_name='ABES'
                              )
            except OSError:
                #no amplitude calibration present yet
                offset = flap.get_data('W7X_ABES',
                              exp_id=exp_id,
                              coordinates={'Time':[-0.1,0.1]},
                              name=dataobject.get_coordinate_object('Signal name').values[0],
                              object_name='ABES',
                              options = {'Amplitude calibration': False}
                              )
            chop = np.mean(d_beam_on.coordinate("Time")[2]-d_beam_on.coordinate("Time")[1])
            if chop<1e-4:
                if np.min(offset.coordinate("Time")[0]) < -1e-6:
                    offset_time = [-0.1,-1e-6]
                    offset = offset.slice_data(slicing={"Time": flap.Intervals(-0.1,-1e-6)})
                    electric_noise = np.var(offset.data, axis=0)
                    dataobject_beam_on.error = np.sqrt(dataobject_beam_on.error**2+electric_noise)
                else:
                    standard_error = np.array([8.52008759e-03, 6.05577881e-03, 5.61784578e-03, 4.01606559e-03,
                           4.80343614e-03, 4.96760886e-03, 4.31857827e-03, 4.61216619e-03,
                           5.07540636e+00, 5.77341359e+00, 5.17886695e+00, 6.67524467e+00,
                           4.10515027e+00, 3.20796480e+00, 4.06675798e+00, 4.07734300e+00,
                           4.09595980e+00, 3.70180586e+00, 4.78968865e+00, 3.76662174e+00,
                           4.02454327e+00, 3.90762905e+00, 5.00180528e+00, 3.25465952e+00,
                           4.32225647e+00, 3.23294377e+00, 5.75256366e+00, 4.56699600e+00,
                           5.02907424e+00, 4.54773953e+00, 6.10916503e+00, 4.75039061e+00,
                           3.76776202e+00, 5.17901390e+00, 4.44659086e+00, 6.81016910e+00,
                           7.95957335e+00, 1.45460916e+01, 1.10220756e+01, 1.85718833e+01])
                    if len(dataobject.get_coordinate_object('Signal name').values)<len(standard_error):
                        ch = [int(chname.split("-")[1])-1 for chname in dataobject.get_coordinate_object('Signal name').values]
                        standard_error = standard_error[ch]
                        if len(ch) == 1:
                            standard_error = standard_error[0]
                    # raise ValueError(str(dataobject.get_coordinate_object('Signal name').values))
                    dataobject_beam_on.error = np.sqrt(dataobject_beam_on.error**2+standard_error)


        if test is True:
            plt.plot(dataobject_beam_on.get_coordinate_object("Time").values,
                     dataobject_beam_on.data)
    
        return dataobject_beam_on

def chopped_signals(exp_ID,timerange=None,signals='ABES-[1-40]',datapath=None,background=False):
    """
    
     *****THIS IS OUTDATED! Use get_clean_abes ****
    Processes chopped ABES measurements. Averages the signal in chopping periods and corrects for the background.
    Returns a data object with all the requested channels. The time resolution will be the chopper period time.
    For background correction the background is linearly interplated between chopping periods to find the background
    when the beam is on. If background is requested the timescale will be the background measurement time vector, 
    otherwise the beam-on time vector.

    Parameters
    ----------
    exp_ID : string
        Experiment ID. E.g. 202410109.062
    timerange : list of two floats, optional
        The time range to process. The default is None, all times will be processed.
    signals : string or string array, optional
        The signals to process. The default is 'ABES-[1-40]'.
    datapath : string, optional
        The data path. The default is None, in this case the default from the flap_defaults.cfg 
        config file will be used.
    background : boolean, optional
        If True the background signal is returned. The default is False.

    Returns
    -------
    flap.dataObject
        The data.

    """    
    
    print("Chopped_signals is outdated. Use get_clean_abes")
    options = {}
    if (datapath is not None):
        options['Datapath'] = datapath
   

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
    off1,off2,off3 = d_beam_off.coordinate('Time')
    beam_on_time = on3[1]-on2[1]
    beam_off_time = off3[1]-off2[1]
    period_time = beam_on_time + beam_off_time
    
    on_start = 0
    on_end = 0
    off_start = 0
    off_end = 0
 
    if (d_beam_on.info['APDCAM_clock_source'] == 'External'):
        if (period_time < 0.01):
            warnings.warn("Fast chopping with external clock, background separation might not be correct.")
        else:
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
      
    d=flap.get_data('W7X_ABES',
                    exp_id = exp_ID,
                    name = signals,
                    options = options,
                    coordinates = {'Time': timerange}
                    )
    d_off = d.slice_data(slicing={'Time':d_beam_off},
                         summing={'Rel. Sample in int(Time)':'Mean'},
                         options={'Regenerate':True,'Partial intervals':False}
                         )

    if (background):
        return d_off
    else:
        d_on = d.slice_data(slicing={'Time':d_beam_on},
                            summing={'Rel. Sample in int(Time)':'Mean'},
                            options={'Regenerate':True,'Partial intervals':False}
                            )
        d_off = d_off.slice_data(slicing={'Time':d_on},
                                 options={'Interpolation':'Linear'}
                                 )
        d_on.data = d_on.data - d_off.data
        return d_on
    
    
def proc_chopsignals(dataobject=None, exp_id=None,timerange=None,
                     samplerange=None, signals='ABES-[1-40]', on_options=None,
                     off_options=None, test=None, options={}):
    """ Calculate signals in beam on and beam/off phases of the measurement and
        correct the beam-on phases with the beam-off. Further information in the comments of 
        proc_chopsignals_single
    """
    naming_conventions = ["Channel", "Signal name", "Device R", "Beam axis"]
    channel_naming = []
    for name in naming_conventions:
        if name in dataobject.coordinate_names():
            channel_naming.append(name)
    
#    if len(channel_naming)==0:
    # Added S. Zoletnik 9 Aug 2024
    if (len(dataobject.shape) == 1):
       processed_data = proc_chopsignals_single(dataobject=dataobject, timerange=timerange,
                                                samplerange=samplerange,
                                                test=test, on_options=on_options,
                                                off_options=off_options,  options=options)
    else:
        # The channels are processed in parallel
        partial_proc_func = partial(proc_chopsignals_single, timerange=timerange, samplerange=samplerange,
                                    test=test, on_options=on_options,
                                    off_options=off_options,  options=options)
        channels = dataobject.get_coordinate_object(channel_naming[0]).values
        num_of_channels = len(channels)
        divide_to = 10
        results =[]

        for index in range(divide_to+1):
            with mp.Pool(int((mp.cpu_count()+1)/2)) as pool:
                curr_channels=channels[index*int(1+num_of_channels/divide_to):(index+1)*int(1+num_of_channels/divide_to)]
                channel_data = [dataobject.slice_data(slicing={channel_naming[0]: channel}) for channel
                                                               in curr_channels]
                del curr_channels
                gc.collect()
                results = results + pool.map(partial_proc_func, channel_data)
                print('Multichannel signal processing finished '+str(int(len(results)/num_of_channels*100))+'%')
        for channel_processed_data in results:
            if not ("processed_data" in locals()):
                processed_data = channel_processed_data
            else:
                if len(processed_data.data.shape) == 1:
                    processed_data.data=np.stack((processed_data.data, channel_processed_data.data), axis=1)
                    processed_data.error=np.stack((processed_data.error, channel_processed_data.error), axis=1)
                else:
                    processed_data.data=np.hstack((processed_data.data, np.expand_dims(channel_processed_data.data, axis=1)))
                    processed_data.error=np.hstack((processed_data.error, np.expand_dims(channel_processed_data.error, axis=1)))

        del results
        del channel_data
        gc.collect()

        processed_data.shape = processed_data.data.shape
        # adding the channel coordinates
        channel_dimension = dataobject.get_coordinate_object(channel_naming[0]).dimension_list
        for coordinate in dataobject.coordinates:
            if coordinate.dimension_list == channel_dimension:
                naming_coord = copy.deepcopy(coordinate)
                naming_coord.dimension_list = [1]
                processed_data.add_coordinate_object(naming_coord)
    return processed_data

def process_chopped_dataobject(dataobject, options={}):
    ''' Processes the input channel data dataobject which has been already
    sliced with a chopper dataobject. It flattens the data and adds the proper
    Time coordinate. If requested it only averages the data over the chopping
    time periods.
    dataobject - the channel data sliced with a chopper DataObject
    options 'Average Chopping Period' - boolean,  whether to average the date
                                        over a chopping interval
    '''
    options_default = {'Average Chopping Period': True}
    options = {**options_default, **options}
    if options['Average Chopping Period'] is True:
        #calculating the error of the beam on part
        reltime_coord = dataobject.get_coordinate_object("Rel. Time in int(Sample)")
        reltime_size = dataobject.data.shape[reltime_coord.dimension_list[0]]
        average = np.nanmean(dataobject.data, axis=reltime_coord.dimension_list[0], keepdims=True)
        beam_on_error = np.nanmean((dataobject.data-average)**2, axis = reltime_coord.dimension_list[0])/(reltime_size-1)
        dataobject.get_coordinate_object("Rel. Time in int(Sample)").dimension_list=[]
        dataobject.get_coordinate_object("Rel. Sample in int(Sample)").dimension_list=[]

        #the averaged density profile - #flap is currently not capable of doing this correctly
        # need to set up a time coordinate correctly
        regenerate_time_sample(dataobject)
        times = copy.deepcopy(dataobject.coordinate("Time")[0])[0]
        samples = copy.deepcopy(dataobject.coordinate("Sample")[0])[0]
        dataobject.coordinates = [dataobject.get_coordinate_object("Time"), dataobject.get_coordinate_object("Sample")]

        #taking the average over the Rel.Time in int(Time) coordinate
        dataobject.data = average[0]
        dataobject.shape = dataobject.data.shape
        dataobject.error = np.sqrt(beam_on_error)
    else:
        add_absolute_time_sample(dataobject)
        times = copy.deepcopy(dataobject.coordinate("Time")[0]).flatten()
        samples = copy.deepcopy(dataobject.coordinate("Sample")[0]).flatten()
        dataobject.data = dataobject.data.flatten()
        times = times[np.logical_not(np.isnan(dataobject.data))] # removing nans due to the padding of the slicer
        samples = samples[np.logical_not(np.isnan(dataobject.data))] # removing nans due to the padding of the slicer
        dataobject.data = dataobject.data[np.logical_not(np.isnan(dataobject.data))]
        dataobject.shape = dataobject.data.shape
        dataobject.error = np.zeros(dataobject.data.shape)
        dataobject.coordinates = [dataobject.get_coordinate_object("Time"), dataobject.get_coordinate_object("Sample")]

    dataobject.get_coordinate_object("Time").values = times
    dataobject.get_coordinate_object("Time").shape = times.shape
    dataobject.get_coordinate_object("Time").dimension_list = [0]
    dataobject.get_coordinate_object("Time").mode.equidistant = False
    dataobject.get_coordinate_object("Time").start=None
    dataobject.get_coordinate_object("Time").step=None

    dataobject.get_coordinate_object("Sample").values = samples
    dataobject.get_coordinate_object("Sample").shape = times.shape
    dataobject.get_coordinate_object("Sample").dimension_list = [0]
    dataobject.get_coordinate_object("Sample").mode.equidistant = False
    dataobject.get_coordinate_object("Sample").start=None
    dataobject.get_coordinate_object("Sample").step=None

    #ordering the data along the time coordinate
    time_data = np.asarray(sorted(zip(dataobject.get_coordinate_object("Time").values,
                                      dataobject.get_coordinate_object("Sample").values,
                                      dataobject.data)))

    dataobject.get_coordinate_object("Time").values = time_data[:,0]
    dataobject.get_coordinate_object("Sample").values = time_data[:,1]
    dataobject.data = time_data[:,2]

    return dataobject

def add_absolute_time_sample(dataobject):
    """ Creates a coordinate 'Time' to the input dataobject from proc_chopsignals. This can be used for slicing the
            data
        INPUT:
            d - dataobject, with coordinates Start Time in int(Sample) and Rek. Time in int(Sample)
        OUTPUT: None, the coordinate is added to the original dataobject
    """
    # Finding the coordinates for the dataobject
    coords = [coord.unit.name for coord in dataobject.coordinates]
    coord_index = 0
    dimension_list_time = []
    for coordinate in coords:
        # will need to set for the time index
        if coordinate == 'Start Time in int(Sample)':
            dimension_list_time = list(dataobject.coordinates[coord_index].dimension_list)+dimension_list_time
        elif coordinate == 'Rel. Time in int(Sample)':
            dimension_list_time = list(dataobject.coordinates[coord_index].dimension_list)+dimension_list_time
        coord_index = coord_index+1
    dimension_list_time = list(np.unique(np.asarray(dimension_list_time)))

    name = 'Time'
    time_coord_value = dataobject.coordinate('Start Time in int(Sample)')[0] +\
                       dataobject.coordinate('Rel. Time in int(Sample)')[0]
    time_coord = flap.Coordinate(name=name, unit='Second', values=time_coord_value, shape=np.shape(time_coord_value),
                 mode=flap.CoordinateMode(equidistant=False), dimension_list=dimension_list_time)
    name = 'Sample'
    sample_coord_value = dataobject.coordinate('Start Sample in int(Sample)')[0] +\
                       dataobject.coordinate('Rel. Sample in int(Sample)')[0]
    sample_coord = flap.Coordinate(name=name, unit='n.a.', values=sample_coord_value, shape=np.shape(sample_coord_value),
                 mode=flap.CoordinateMode(equidistant=False), dimension_list=dimension_list_time)
    dataobject.add_coordinate_object(time_coord)
    dataobject.add_coordinate_object(sample_coord)

def get_pearson_matrix(dataobject, options={}):
    """ This function calculates and adds the pearson matrix for the channels of the backgorund corrected
        A-BES dataobject.
        INPUT:
            dataobject: The background corrected A-BES data, should have only a Channel/Signal name dimension 
            and a 'Time' dimension.
            options:
                "Time Window": If defined, it is the time window for which the fluctuations should be analzed to
                calculate the Pearson matrix. If undefined, then the whole time interval of datobject is used
        OUTPUT:
            The Pearson matrix as a dataobject which is the size number_of_channels*number_of_channels
    """
    options_default = {'Time Window': None}
    options = {**options_default, **options}
    if options['Time Window'] is None:
        options['Time Window'] = np.array([np.min(dataobject.coordinate('Time')[0]),
                                           np.max(dataobject.coordinate('Time')[0])])

    dataobject = copy.deepcopy(dataobject) # so that the original one is not accidentally overwritten
    dataobject.error = None

    # Creating the matrix for the calculation of the correlation
    # Getting the channel dimension:
    data_coord_list = np.array([coord.unit.name for coord in dataobject.coordinates])
    if 'Channel' in data_coord_list:
        ch_coord = 'Channel'
        channel_coordinate_dim = dataobject.get_coordinate_object('Channel').dimension_list[0]
        channel_names = np.unique(dataobject.get_coordinate_object('Channel').values)
    else:
        ch_coord = 'Signal name'
        channel_coordinate_dim = dataobject.get_coordinate_object('Signal name').dimension_list[0]
        channel_names = np.unique(dataobject.get_coordinate_object('Signal name').values)

    # Slicing along the channels and flattening the data
    channel_id = 0
    for channel in channel_names:
        channel_data = dataobject.slice_data(slicing={ch_coord: channel}).data
        channel_data = channel_data.flatten()
        if channel_id == 0:
            corrcoeff_input = np.zeros([len(channel_names), len(channel_data)])
        corrcoeff_input[channel_id, :] = channel_data
        channel_id += 1

    # Calculating the correlation matrix
    pearson_matrix = np.corrcoef(corrcoeff_input)

#    channel_coord = flap.Coordinate(name=ch_coord, unit='', values=channel_names, shape=np.shape(channel_names),
#                                    mode=flap.CoordinateMode(equidistant=False), dimension_list=[0, 1])
    coordinates = []
    #adding all spatial coordinates that are in the original dataobject
    for coord in data_coord_list:
        if dataobject.get_coordinate_object(coord).dimension_list[0] == channel_coordinate_dim:
            coordinates = coordinates + [dataobject.get_coordinate_object(coord)]
            other_axis = copy.deepcopy(dataobject.get_coordinate_object(coord))
            other_axis.dimension_list = [1]
            other_axis.unit.name = other_axis.unit.name + ' 2'
            coordinates = coordinates + [other_axis]


    pearson_dataobject = flap.DataObject(data_array=np.asarray(pearson_matrix),
                                         data_unit=flap.Unit(name='Correlation', unit='1'),
                                         coordinates=coordinates, exp_id=dataobject.exp_id,
                                         data_title='Pearson Correlation Matrix', data_shape=pearson_matrix.shape,
                                         error=None)

    return pearson_dataobject

def read_chopshift(shotID):
    "Reads by how much the chopping has to be shifted to align the data"
    location = os.path.dirname(os.path.realpath(__file__))
    location = os.path.join(location, 'w7x_chop_shift')
    with open(location, "r", encoding='utf-8') as f:
        lines = list(f)
        chopshift={}
        for line in lines[1:]:
            data = line.split('\t')
            chopshift[data[0]] = np.asarray([float(data[1]), float(data[-1].split('\n')[0])])
        if shotID in chopshift:
            delay = chopshift[shotID]
        else:
            delay = np.asarray([flap.config.get("Module W7X_ABES","Start delay"),
                                flap.config.get("Module W7X_ABES","End delay")])
    return delay

def write_chopshift(shotID, start, end):
    'Writes the data into the w7x_chop_shift file'
    location = os.path.dirname(os.path.realpath(__file__))
    location = os.path.join(location, 'w7x_chop_shift')
    try:
        os.rename(location, location+'old')
    except FileNotFoundError:
        print(f"{location} not found, probably disk is full. Free up space and press Enter.")
        input()
    separator = '\t'
    found_shot = False
    with open(location, 'w', encoding='utf-8') as fout:
        with open(location+'old', "r", encoding='utf-8') as fin:
            for line in fin:
                data = line.split(separator)
                if data[0] == shotID:
                    data[1] = str(start)
                    data[-1] = str(end)+'\n'
                    found_shot = True
                dataout = separator.join(data)
                fout.write(dataout)
        if found_shot == False:
            data = [shotID, str(start), '',str(end)+'\n']
            dataout = separator.join(data)
            fout.write(dataout)
    os.remove(location+'old')

  
def register(data_source=None):
    flap.register_data_source('W7X_ABES', get_data_func=w7x_abes_get_data, add_coord_func=add_coordinate)
    from .cxrs_main import w7x_abes_cxrs_get_data, cxrs_add_coordinate
    flap.register_data_source('W7X_ABES_CXRS', get_data_func=w7x_abes_cxrs_get_data, add_coord_func=cxrs_add_coordinate)

