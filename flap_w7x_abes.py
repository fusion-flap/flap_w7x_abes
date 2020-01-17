# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:23:49 2018

@author: Zoletnik

This is the flap module for W7-X alkali BES diagnostic
"""

import os.path
from decimal import *
import numpy as np
import copy
import h5py
import matplotlib.pyplot as plt

import flap
from . import spatcal

if (flap.VERBOSE):
    print("Importing flap_w7x-abes")

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
        retval['TriggerTime'] = Decimal(xml.get_element('System','SystemTrigger')['Value'])
    retval['APDCAM_state'] = int(xml.get_element('APDCAM','State')['Value'])
    if (retval['APDCAM_state'] is 1):
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
        return data_arr

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
        try:
            infile = open(fn, 'r', encoding='ascii')
        except Exception as e:
            raise e
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
                    cal_timerange = [float(line[1]), float(line[2])]()
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
        infile.close()
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
    # Doing the calibration
    if (data_arr.dtype.kind != 'f'):
        data_arr = float(data_arr)
    for i in range(len(index_start)):
        # Collecting the calibration factors andd errors
        calfac_act = np.empty(len(signal_proc),dtype=float)
        calfac_act_err = np.empty(len(signal_proc),dtype=float)
        for i_ch in range(len(signal_proc)):
            try:
                sig = signal_proc[i_ch].encode("ascii")
                if (len(sig) < 7):
                    sig += b' '
                i_cal = cal_channels[i].index(sig)
            except IndexError:
                if (signal_proc[i_ch][0:5] == 'ABES-'):
                    raise ValueError("No calibration data for signal "+signal_proc[i_ch])
            if (data_arr.ndim == 2):
                data_arr[index_start[i]:index_stop[i], i_ch] /= calfac[i][i_cal]
            else:
                data_arr[index_start[i]:index_stop[i]] /= calfac[i][i_cal]
    return data_arr

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
        'Start delay' : Start delay relative to the phase start in microsec
        'End delay' : Delay at the end of the phase in microsec. Positive means later times.
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
    chop_clk_in_sample = round(config['APDCAM_f_sample'] / config['Chopper clock'])
    # This is temporary, has to be corrected with the flight time and APDCAM trigger delay
    instrument_delay = 2/Decimal(1000000)
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

    try:
        start_delay = float(options['Start delay'])
        start_delay_sample = int(round(start_delay *1e-6 / float(config['APDCAM_sampletime'])))
    except KeyError:
        start_delay_sample = 0
    except ValueError:
        raise ValueError("Invalid chopper start delay.")

    try:
        stop_delay = float(options['End delay'])
        stop_delay_sample = int(round(stop_delay * 1e-6 / float(config['APDCAM_sampletime'])))
    except KeyError:
        stop_delay_sample = 0
    except ValueError:
        raise ValueError("Invalid chopper end delay.")

    if (read_samplerange is not None):
        start_sample = start_clk * chop_clk_in_sample
        period_time_sample = period_time_clk * chop_clk_in_sample
        if (start_sample < read_samplerange[0]):
            start_ind = int((read_samplerange[0] - start_sample)
                            / period_time_sample + 1)
            start_clk += period_time_clk * start_ind
            stop_clk += period_time_clk * start_ind
        stop_sample = stop_clk * chop_clk_in_sample
        n_period =  int((read_samplerange[1] - stop_sample) / period_time_sample)
        if (n_period <= 0):
            raise ValueError("No chopper intervals in time (sample) range.")
    else:
        n_period = int(round(config['APDCAM_samplenumber'] /  (period_time_clk * chop_clk_in_sample)))

    mode = flap.CoordinateMode(equidistant=True, range_symmetric=False)

    c_time = copy.deepcopy(flap.Coordinate(name='Time',
                                           unit='Second',
                                           mode=mode,
                                           start=start_clk * clock_unit + config['APDCAM_starttime']
                                                 + instrument_delay
                                                 + start_delay_sample * config['APDCAM_sampletime'],
                                           step=period_time_clk*clock_unit,
                                           value_ranges=[0, float((stop_clk-start_clk) * clock_unit)
                                                           + (stop_delay_sample - start_delay_sample)
                                                             * float(config['APDCAM_sampletime'])],
                                           dimension_list=[0]
                                           ))
    c_sample = copy.deepcopy(flap.Coordinate(name='Sample',
                                             unit='n.a.',
                                             mode=mode,
                                             start=round(start_clk * chop_clk_in_sample
                                                    + instrument_delay/config['APDCAM_sampletime'])
                                                    + start_delay_sample,
                                             step=round(period_time_clk*chop_clk_in_sample) ,
                                             value_ranges=[0,round((stop_clk-start_clk) * chop_clk_in_sample)
                                                           + stop_delay_sample - start_delay_sample],
                                             dimension_list=[0]
                                             ))

    d = copy.deepcopy(flap.DataObject(data_shape=[n_period],
                        data_unit=flap.Unit(name='Interval', unit='n.a.'),
                        coordinates=[c_sample, c_time],
                        data_source = 'W7X_ABES',
                        exp_id = config['ShotID']
                        ))

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
        'Scaling':  'Digit'
                    'Volt'
        'Offset timerange': Time range for offset subtraction
        'Datapath': Data path (string)
        'Calibration': True/False do/don't do amplitude calibration of the data
        'Calib. path': Calibration directory name
        'Calib. file': Calibration cld file name.
        For further options see Chopper_times see shopper_timing_data()

    """
    if (exp_id is None):
        raise ValueError('exp_id should be set for W7X ABES.')

    default_options = {'Datapath':'data',
                       'Scaling':'Digit',
                       'Offset timerange': None,
                       'Amplitude calibration': False,
                       'Amplitude calib. path': 'cal',
                       'Amplitude calib. file': None,
                       'Phase' : None,
                       'State' : None,
                       'Start delay': 0,
                       'End delay': 0,
                       'Spatial calibration': False
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
    try:
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
            if (coord.unit.name is 'Time'):
                if (coord.mode.equidistant):
                    read_range = [float(coord.c_range[0]),float(coord.c_range[1])]
                    if (read_range[1] <= read_range[0]):
                        raise ValueError("Invalid read timerange.")
                else:
                    raise NotImplementedError("Non-equidistant Time axis is not implemented yet.")
                break
            if coord.unit.name is 'Sample':
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
            try:
                f = open(fn,"rb")
            except OSError:
                raise OSError("Error opening file: " + fn)
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


    ndata = int(read_samplerange[1] - read_samplerange[0] + 1)

    if (no_data is False):
        if (len(ADC_proc) is not 1):
            data_arr = np.empty((ndata, len(ADC_proc)), dtype=dtype)
        for i in range(len(ADC_proc)):
            fn = os.path.join(datapath, "Channel_{:03d}.dat".format(ADC_proc[i] - 1))
            try:
                f = open(fn,"rb")
            except OSError:
                raise OSError("Error opening file: " + fn)

            try:
                f.seek(int(read_samplerange[0]) * 2, os.SEEK_SET)
            except Exception:
                raise IOError("Error reading from file: " + fn)

            if (len(ADC_proc) is 1):
                try:
                    data_arr = np.fromfile(f, dtype=np.int16, count=ndata)
                except Exception:
                    raise IOError("Error reading from file: " + fn)
                if (scale_to_volts):
                    data_arr = ((2 ** config['APDCAM_bits'] - 1) - data_arr) \
                                / (2. ** config['APDCAM_bits'] - 1) * 2
                else:
                    data_arr = (2 ** config['APDCAM_bits'] - 1) - data_arr
                if (offset_timerange is not None):
                        data_arr -= offset_data[i]
            else:
                try:
                    d = np.fromfile(f, dtype=np.int16, count=ndata)
                except Exception:
                    raise IOError("Error reading from file: " + fn)
                if (scale_to_volts):
                    d = ((2 ** config['APDCAM_bits'] - 1) - d) \
                                / (2. ** config['APDCAM_bits'] - 1) * 2
                else:
                    d = (2 ** config['APDCAM_bits'] - 1) - d
                if (offset_timerange is not None):
                        d -= offset_data[i]
                data_arr[:,i] = d
        f.close

        try:
            data_arr = calibrate(data_arr, signal_proc, read_range, exp_id=exp_id, options=_options)
        except Exception as e:
            raise e
        data_dim = data_arr.ndim    
    else:
        if (len(ADC_proc) is not 1):
            data_dim = 2
        else:
            data_dim = 1

    coord = [None]*data_dim*2
    c_mode = flap.CoordinateMode(equidistant=True)
    coord[0] = copy.deepcopy(flap.Coordinate(name='Time',
                                             unit='Second',
                                             mode=c_mode,
                                             start=read_range[0],
                                             step=config['APDCAM_sampletime'],
                                             dimension_list=[0])
                             )
    coord[1] = copy.deepcopy(flap.Coordinate(name='Sample',
                                             unit='n.a.',
                                             mode=c_mode,
                                             start=read_samplerange[0],
                                             step=1,
                                             dimension_list=[0])
                             )
    if (data_dim > 1):
        ch_proc = []
        for ch in signal_proc:
            if (ch[0:4] != 'ABES'):
                ch_proc = []
                break
            ch_proc.append(int(ch[5:]))
        if (ch_proc != []):
            c_mode = flap.CoordinateMode(equidistant=False)
            coord[2] = copy.deepcopy(flap.Coordinate(name='Channel',
                                                     unit='n.a.',
                                                     mode=c_mode,
                                                     shape=len(ch_proc),
                                                     values=ch_proc,
                                                     dimension_list=[1])
                                 )
        coord[3] = copy.deepcopy(flap.Coordinate(name='Signal name',
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
                        data_unit=data_unit,
                        coordinates=coord,
                        exp_id=exp_id,
                        data_title=data_title,
                        info={'Options':_options},
                        data_source="W7X_ABES")
    if options['Spatial calibration'] is True:
        # Getting the spatial calibration
        d = add_coordinate(d, ['Device R', 'Beam axis'],
                           options={"spatcal_dir": flap.config.get("Module W7X_ABES","Spatial calibration directory")})
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
    default_options = {'spatcal_dir':'', "Channels":''}
    _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES')

    if exp_id is None:
        exp_spatcal = spatcal.ShotSpatCal(data_object.exp_id, options=_options)
    else:
        exp_spatcal = spatcal.ShotSpatCal(exp_id, options=_options)

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
            raise ValueError("Rel Time in int(Sample) is not constant.")
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
        ct.unit.name='Sample'
        d.del_coordinate('Rel. Sample in int(Sample)')
    try:
        d.del_coordinate('Interval(Sample)')
        d.del_coordinate('Interval(Sample) sample index')
    except ValueError:
        pass

def proc_chopsignals_single(exp_id=None,timerange=None,signals='ABES-1', on_options=None,
                             off_options=None, test=None, dataobject=None, options={}):
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
    options_default = {'Average Chopping Period': True}
    options = {**options_default, **options}


    # Obtaining the chopper data
    if dataobject is not None:
        exp_id = dataobject.exp_id
        
    o = copy.deepcopy(on_options)
    if 'W7X_ABES' not in flap.list_data_sources():
        register()
    o.update({'State':{'Chop': 0, 'Defl': 0}})
    if timerange is None:
        timerange = [np.min(dataobject.coordinate('Time')), np.max(dataobject.coordinate('Time'))]
    d_beam_on=flap.get_data('W7X_ABES',
                            exp_id=exp_id,
                            name='Chopper_time',
                            coordinates={'Time':timerange},
                            options=o,\
                            object_name='Beam_on',
                            )
    o = copy.deepcopy(off_options)
    o.update({'State':{'Chop': 1, 'Defl': 0}})
    d_beam_off=flap.get_data('W7X_ABES',
                            exp_id=exp_id,
                            name='Chopper_time',
                            coordinates={'Time':timerange},
                            options=o,\
                            object_name='Beam_off',
                            )

<<<<<<< HEAD
    if (test):
        from matplotlib import pyplot as plt
        plt.close('all')
        if dataobject is None:
            flap.plot('ABES', axes='Time', plot_options={'marker': 'o'})
        else:
            temp = dataobject.slice_data(slicing={'Signal name':'ABES-20'})
            temp.plot(axes='Time', plot_options={'marker': 'o'})
#        flap.plot('ABES',axes='Time',plot_type='scatter')
        d_beam_on.plot(plot_type='scatter', axes=['Time', 2], options={'Force': True,'All': True})
        d_beam_off.plot(plot_type='scatter', axes=['Time', 0.1], options={'Force': True,'All': True})

=======
>>>>>>> a8799f7d6f673dab3d23a4dbb23ddadde63ccaba
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
            add_absolute_time(d)

        flap.add_data_object(d,'ABES_on')
        
        d = flap.slice_data('ABES',slicing={'Sample':d_beam_off})
        d = d.slice_data(summing={'Rel. Sample in int(Sample)':'Mean'})
        regenerate_time_sample(d)    
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
        flap.add_data_object(d,'ABES')
        
#        # error approximation
#        d_beam_off = flap.get_data_object('ABES_off')
##        beam_off_data = d_beam_off.slice_data(summing={'Rel. Sample in int(Sample)':'Mean'})
##        regenerate_time_sample(beam_off_data)
#        beam_off_data = beam_off_data.slice_data(slicing={'Time':d_beam_off},options={'Inter':'Linear'})
#        background_error = np.average((d_beam_off.data-beam_off_data.data.reshape(np.shape(d_beam_off.data)))**2, axis=0)\
#                           *len(beam_off_data.data)/(len(beam_off_data.data)-1)

        flap.delete_data_object(['ABES_on','ABES_off','Beam_on','Beam_off'],exp_id=exp_id)
        if (test):
            plt.figure()
            flap.plot('ABES',axes='Time')
            
        return d
    else:
        # in this case the passed dataobject is used and the only the copper data is obtained from file
        dataobject_beam_on = dataobject.slice_data(slicing={'Sample': d_beam_on})
        # For dataobject_beam_on.data the 0 dimension is along a constant 'Start Time in int(Time)' and 
        # "Rel. Time in int(Time)" varies. For dimension 1 it is 'Start Time in int(Time)' that varies
        dataobject_beam_on = process_chopped_dataobject(dataobject_beam_on, options=options)
    
        dataobject_beam_off = dataobject.slice_data(slicing={'Sample': d_beam_off}, options={'Partial intervals':False})
        dataobject_beam_off = process_chopped_dataobject(dataobject_beam_off, options={'Average Chopping Period': True})
    
        dataobject_background = dataobject_beam_off.slice_data(slicing={'Time': dataobject_beam_on},
                                                               options={'Inter': 'Linear'})
        if (test):
            from matplotlib import pyplot as plt
            plt.figure()
            plt.scatter(dataobject_beam_on.get_coordinate_object("Time").values, dataobject_beam_on.data)
            plt.scatter(dataobject_beam_off.get_coordinate_object("Time").values, dataobject_beam_off.data)
            plt.plot(dataobject.coordinate("Time")[0], dataobject.data)            
        
        dataobject_beam_on.data -= dataobject_background.data.reshape(np.shape(dataobject_beam_on.data))
    
        # calculating the error for the beam off part
        dataobject_beam_off_error = copy.deepcopy(dataobject_beam_off)
        dataobject_beam_off_error.data = np.sqrt(dataobject_beam_off.error)
        background_error = dataobject_beam_off_error.slice_data(slicing={'Time': dataobject_beam_on}, options={'Inter': 'Linear'})
        background_error = background_error.data.reshape(np.shape(dataobject_beam_on.data))
        dataobject_beam_on.error = np.asarray(np.sqrt(dataobject_beam_on.error**2 + background_error**2))
    
        if (test):
            plt.plot(dataobject_beam_on.get_coordinate_object("Time").values,
                     dataobject_beam_on.data)
    
        return dataobject_beam_on

def proc_chopsignals(exp_id=None,timerange=None,signals='ABES-[1-40]', on_options=None,
                     off_options=None, test=None, dataobject=None, options={}):
    """ Calculate signals in beam on and beam/off phases of the measurement and
        correct the beam-on phases with the beam-off. Further information in the comments of 
        proc_chopsignals_single
    """
    
    # checking, whether dataobject has the data for multiple channels
    naming_conventions = ["Channel", "Signal name", "Device R", "Beam axis"]
    channel_naming = []
    for name in naming_conventions:
        if name in dataobject.coordinate_names():
            channel_naming.append(name)
    
    if len(channel_naming)==0:
       processed_data = proc_chopsignals_single(dataobject=dataobject, timerange=timerange,
                                          test=test, on_options=on_options,
                                          off_options=off_options,  options=options)
    else:
        channels = dataobject.get_coordinate_object(channel_naming[0]).values
        for channel in channels:
            channel_data = dataobject.slice_data(slicing={channel_naming[0]: channel})
            if not ("processed_data" in locals()):
                processed_data = proc_chopsignals_single(dataobject=channel_data,
                                                         timerange=timerange,
                                                         test=test, on_options=on_options,
                                                         off_options=off_options,  options=options)
            else:
                partial_processed_data = proc_chopsignals_single(dataobject=channel_data, timerange=timerange,
                                          test=test, on_options=on_options,
                                          off_options=off_options,  options=options)
                if len(processed_data.data.shape) == 1:
                    processed_data.data=np.stack((processed_data.data, partial_processed_data.data), axis=1)
                    processed_data.error=np.stack((processed_data.error, partial_processed_data.error), axis=1)
                else:
                    processed_data.data=np.hstack((processed_data.data, np.expand_dims(partial_processed_data.data, axis=1)))
                    processed_data.error=np.hstack((processed_data.error, np.expand_dims(partial_processed_data.error, axis=1)))

#                del partial_processed_data
#                gc.collect()
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
        beam_on_error = np.nansum((dataobject.data-average)**2, axis = reltime_coord.dimension_list[0])/(reltime_size-1)
        dataobject.get_coordinate_object("Rel. Time in int(Sample)").dimension_list=[]
        dataobject.get_coordinate_object("Rel. Sample in int(Sample)").dimension_list=[]

        #the averaged density profile - #flap is currently not capable of doing this correctly
        # need to set up a time coordinate correctly
        regenerate_time_sample(dataobject)
        times = copy.deepcopy(dataobject.coordinate("Time")[0])[0]
        dataobject.coordinates = [dataobject.get_coordinate_object("Time")]

        #taking the average over the Rel.Time in int(Time) coordinate
        dataobject.data = average[0]
        dataobject.shape = dataobject.data.shape
        dataobject.error = np.sqrt(beam_on_error)
    else:
        add_absolute_time(dataobject)
        times = copy.deepcopy(dataobject.coordinate("Time")[0]).flatten()
        dataobject.data = dataobject.data.flatten()
        times = times[np.logical_not(np.isnan(dataobject.data))] # removing nans due to the padding of the slicer
        dataobject.data = dataobject.data[np.logical_not(np.isnan(dataobject.data))]
        dataobject.shape = dataobject.data.shape
        dataobject.error = np.zeros(dataobject.data.shape)
        dataobject.coordinates = [dataobject.get_coordinate_object("Time")]

    dataobject.get_coordinate_object("Time").values = times
    dataobject.get_coordinate_object("Time").shape = times.shape
    dataobject.get_coordinate_object("Time").dimension_list = [0]
    dataobject.get_coordinate_object("Time").mode.equidistant = False
    dataobject.get_coordinate_object("Time").start=None
    dataobject.get_coordinate_object("Time").step=None
    
    #ordering the data along the time coordinate
    time_data = np.asarray(sorted(zip(dataobject.get_coordinate_object("Time").values,
                                      dataobject.data)))
    dataobject.get_coordinate_object("Time").values = time_data[:,0]
    dataobject.data = time_data[:,1]
    
    return dataobject

def add_absolute_time(dataobject):
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
    dataobject.add_coordinate_object(time_coord)

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
    
def register(data_source=None):
    flap.register_data_source('W7X_ABES', get_data_func=w7x_abes_get_data, add_coord_func=add_coordinate)

