# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 18:36:43 2021

@author: mvecsei
"""
import sys
import os
from functools import partial
import xml.etree.ElementTree as ET
import time
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import signal

import flap
#import flap_apdcam_op1 as flap_apdcam
import flap_apdcam
import flap_w7x_abes
flap_apdcam.register()
flap_w7x_abes.register()


def read_fibre_conf(datapath, shotID, dataobject):
    # reading the fibre data
    xml_file = os.path.join(datapath, shotID+"_config.xml")
    tree = ET.parse(xml_file)
    root = tree.getroot()

    fibre_ID = copy.deepcopy(dataobject.get_coordinate_object('Signal name'))
    fibre_ID.values = np.array(fibre_ID.values, dtype=object) #so that long fibre names are not truncated
    fibre_ID.unit.name = 'Fibre ID'
    detector_type = copy.deepcopy(dataobject.get_coordinate_object('Signal name'))
    detector_type.values = np.array(detector_type.values, dtype=object)
    detector_type.unit.name = 'Detector type'
    signal_name_list = dataobject.get_coordinate_object('Signal name').values
    for child in root.find('Optics'):
        signum = np.where(signal_name_list == child.tag)[0][0]
        fibre_ID_curr = '-'.join(child.attrib['Value'].split('-')[:2])
        fibre_ID.values[signum] = fibre_ID_curr
        detector_type_curr = child.attrib['Value'].split('-')[-1]
        detector_type.values[signum] = detector_type_curr

    return detector_type, fibre_ID

def downsample(data_array, maxpoints=100):
    stepping = len(data_array)//maxpoints
    data_array_new = data_array[::stepping]
    data_array_new=data_array[:len(data_array)//maxpoints*maxpoints]
    data_array_new=data_array_new.reshape(maxpoints,len(data_array_new)//maxpoints)
    return np.mean(data_array_new, axis=1), np.min(data_array_new, axis=1), np.max(data_array_new, axis=1)

def plot_shot(shotID, resample=1e4, smoothen=0, channel_name=None, canvas=None,
              plot_rawdata=True, plot_powerspect=False):
	    
    datapath = os.path.join('/data', shotID)
    
    print ("Reading data")
    # dataobject = flap.get_data('APDCAM',
    #               name=channel,
    #               options={'Datapath':datapath,'Scaling':'Volt', "Camera type": "APDCAM-10G_4x16",
    #                         "Camera version": 1},
    #               object_name=channel
    #               )
    try:
        resample = int(resample)
    except TypeError:
        resample = None
    if channel_name is None:
        channel = "ADC*"
        if resample == 0 or resample is None:
            dataobject = flap.get_data('APDCAM',
                          name=channel,
                          options={'Datapath':datapath,'Scaling':'Volt',
                                   "Camera type": "APDCAM-10G_FC"},
                          coordinates={"Time":[3.0,5]},
                          object_name=channel
                          )
        else:
            dataobject = flap.get_data('APDCAM',
                          name=channel,
                          options={'Datapath':datapath,'Scaling':'Volt',
                                   "Camera type": "APDCAM-10G_FC",
                                   "Resample": resample},
                          object_name=channel
                          )
            # dataobject = dataobject.slice_data(slicing={"Time":flap.Intervals(1.0,1.0001)})
        # getting time shift
        xml_file = os.path.join(datapath,f"{shotID}_config.xml")
        xml = flap.FlapXml()
        xml.read_file(xml_file)
        delay = flap_w7x_abes.abes_get_config(xml)["TriggerTime"]
        dataobject.get_coordinate_object("Time").start += float(delay)

        # reading the fibre data from the xml
        detector_type, fibre_ID = read_fibre_conf(datapath, shotID, dataobject)
        dataobject.add_coordinate_object(detector_type)
        dataobject.add_coordinate_object(fibre_ID)
        
        # fibre_array = [fid.split('-')[-1] for fid in dataobject.get_coordinate_object('Fibre ID').values]
        # optical_channel = [fid.split('-')[0] for fid in dataobject.get_coordinate_object('Fibre ID').values]
        # adc_array = dataobject.get_coordinate_object('Signal name').values
        # res = np.asarray([adc_array, fibre_array, optical_channel])
        # res = res[:,np.argsort(res[0,:])]
        # print(res.transpose())
        
        
        print("Plotting data")
        
        coord_ax = dataobject.get_coordinate_object('ADC Channel').values
        #coord_ax = dataobject.get_coordinate_object('Channel').values

        if plot_rawdata is True:
            
            #get the figure order
            fibres = [name.split("-")[0].zfill(2) for name in dataobject.get_coordinate_object("Fibre ID").values]
            fibreindex = sorted(range(len(fibres)), key=fibres.__getitem__)
            plotindex_to_channel = np.zeros(len(fibres))

            maxval = np.max(np.max(dataobject.data))
            
            if canvas == None:
                fig = plt.figure(figsize=(8,6))
                fig.subplots_adjust(wspace=0,hspace=0)               
                fig.show()                
            else:
                fig = canvas.fig
                fig.subplots_adjust(wspace=0,hspace=0)

	    
            channels = np.arange(64)+1
            all_axes = []
            for channel in channels:
                if channel == 1:
                    channel_data = dataobject.slice_data(slicing={'ADC Channel': channel})
                    index = np.where(coord_ax == channel)[0][0]
                    plotindex = np.where(fibreindex == index)[0][0]+1
                    plotindex_to_channel[plotindex] = channel
                    #channel_data = dataobject.slice_data(slicing={'Channel': channel})
    
                else:
                    try:
                        index = np.where(coord_ax == channel)[0][0]
                        channel_data.data = dataobject.data[:, index]
                        channel_data.get_coordinate_object("Detector type").values=\
                            [dataobject.get_coordinate_object("Detector type").values[index]]
                        channel_data.get_coordinate_object("Fibre ID").values=\
                                [dataobject.get_coordinate_object("Fibre ID").values[index]]
                        plotindex = np.where(fibreindex == index)[0][0]+1
                        plotindex_to_channel[plotindex] = channel
                    except Exception as e:
                        channel_data.data =  channel_data.data*np.nan
                        plotindex_to_channel[63] = channel
                if canvas == None:
                    ax = plt.subplot(8,8,plotindex)
                else:
                    ax = canvas.fig.add_subplot(8,8,plotindex)
                color = 'black'
                if 'MPPC' in channel_data.coordinate('Detector type')[0]:
                    color = 'tab:green'
                if 'SM' in channel_data.coordinate('Fibre ID')[0][0]:
                    color = 'tab:blue'
                if ('Null' in channel_data.coordinate('Fibre ID')[0][0]) or ('Unknown' in channel_data.coordinate('Fibre ID')[0][0]):
                    color = 'tab:red'
                if ("L" in channel_data.coordinate('Fibre ID')[0][0]):
                    color = 'tab:purple'
                if ("R" in channel_data.coordinate('Fibre ID')[0][0]):
                    color = 'tab:orange'
                [downs_data, downs_data_min, downs_data_max] = downsample(channel_data.data)
                [downs_time, dows_time_min, downs_time_max] = downsample(channel_data.coordinate('Time')[0])
                if canvas == None:
                    plt.errorbar(downs_time, downs_data, np.array([downs_data-downs_data_min,downs_data_max-downs_data]), color=color,
                            alpha = 0.75, label =  channel_data.coordinate('Detector type')[0][0]+str(channel)+"-"+channel_data.coordinate('Fibre ID')[0][0].split('-')[-1] )
                    
                    #channel_data.plot(axes='Time',
                    #                  plot_options={'color':color, 'alpha':0.75, 'label': channel_data.coordinate('Detector type')[0][0]+str(channel)+"-"+channel_data.coordinate('Fibre ID')[0][0].split('-')[-1] },
                    #                  options={'Maxpoints':100})
                    # channel_data.plot(axes='Time',
                    #                   plot_options={'color':color, 'alpha':0.75, 'label': "ADC"+str(channel)+"-"+channel_data.coordinate('Fibre ID')[0][0].split('-')[-1] },
                    #                   options={'Maxpoints':100})
                    #plt.ylim([0,maxval*1.2])
                    plt.ylabel('')
                    plt.yticks([])
                    plt.tick_params(axis='y', direction='in', pad=-22)
                    plt.xlabel('')
                    plt.xticks([])
                    plt.legend()
                    all_axes += [ax]
                    #plt.show()
                    #plt.pause(0.1)
                    fig.canvas.draw_idle()
                    fig.canvas.flush_events()
                    time.sleep(0.01)
                else:

                    ax.errorbar(downs_time, downs_data, np.array([downs_data-downs_data_min,downs_data_max-downs_data]), color=color,
                            alpha = 0.75, label =  channel_data.coordinate('Detector type')[0][0]+str(channel)+"-"+channel_data.coordinate('Fibre ID')[0][0].split('-')[-1] )
                    
                    ax.set_ylabel('')
                    ax.set_yticks([])
                    ax.tick_params(axis='y', direction='in', pad=-22)
                    ax.set_xlabel('')
                    ax.set_xticks([])
                    ax.legend()
                    all_axes += [ax]
                    fig.canvas.draw_idle()
                    fig.canvas.flush_events()
                    time.sleep(0.01)

            onclick = partial(onclick_all, dataobject, all_axes, plotindex_to_channel)
            fig.canvas.mpl_connect("button_press_event", onclick)
            fig.text(0.5, 0.04, 'Time (s)', ha='center')
            fig.text(0.04, 0.5, 'Signal (V)', va='center', rotation='vertical')
            #plt.get_current_fig_manager().full_screen_toggle()
            fig.suptitle(shotID)
            #fig.tight_layout()
            #fig.subplots_adjust(wspace=0,hspace=0)
            
            # fig.show()
            # plt.show(block=True)        
    else:
        resample = int(resample)
        if resample == 0:
            dataobject = flap.get_data('APDCAM',
                          name=channel_name,
                          options={'Datapath':datapath,'Scaling':'Volt',
                                   "Camera type": "APDCAM-10G_FC"},
                          object_name=channel_name
                          )
        else:
            dataobject = flap.get_data('APDCAM',
                          name=channel_name,
                          options={'Datapath':datapath,'Scaling':'Volt',
                                   "Camera type": "APDCAM-10G_FC",
                                   "Resample": resample},
                          object_name=channel_name
                          )
        if smoothen == 1: 
            dataobject = dataobject.filter_data(options={"Type":"Int", "Tau":0.002})

        # dataobject.plot(axes='Time',
        #                   plot_options={'alpha':0.75})
        plt.plot(dataobject.coordinate("Time")[0], dataobject.data)
        # plt.title(f"ADC{channel_name} - Fibre {dataobject.coordinate('Fibre ID')[0][0].split('-')[-1]} - {dataobject.coordinate('Detector type')[0][0]}")
        plt.ylabel('Signal (V)')
        plt.xlabel('Time (s)')
        plt.show(block=True)


def onclick_all(dataobject, axes, plotindex_to_channel, event):
    channel = 1 
    for ax in axes:
        if event.inaxes == ax:
            fig = plt.figure()
            channel_data = dataobject.slice_data(slicing={'ADC Channel': channel})
            color = 'black'
            if 'MPPC' in channel_data.coordinate('Detector type')[0]:
                color = 'tab:green'
            if 'SM' in channel_data.coordinate('Fibre ID')[0][0]:
                color = 'tab:blue'
            if ('Null' in channel_data.coordinate('Fibre ID')[0][0]) or ('Unknown' in channel_data.coordinate('Fibre ID')[0][0]):
                color = 'tab:red'
            if ("L" in channel_data.coordinate('Fibre ID')[0][0]):
                color = 'tab:purple'
            if ("R" in channel_data.coordinate('Fibre ID')[0][0]):
                color = 'tab:orange'
            # plt.plot(channel_data.coordinate('Time')[0], channel_data.data,
            #             alpha = 0.75, color = color)
            channel_data.plot(axes='Time',
                              plot_options={'color':color, 'alpha':0.75},
                              options={"All points": True})
            plt.title(f"ADC{channel} - Optical channel {channel_data.coordinate('Fibre ID')[0][0].split('-')[0]} - Fibre {channel_data.coordinate('Fibre ID')[0][0].split('-')[-1]} - {channel_data.coordinate('Detector type')[0][0]}")
            plt.ylabel('Signal (V)')
            plt.xlabel('Time (s)')
            #plt.legend()
            plt.tight_layout()
            plt.show()
        channel += 1
  
def smooth(y,box_pts):
    box=np.ones(box_pts)/box_pts
    y_smooth=np.convolve(y,box,mode='same')
    return np.array(y_smooth)

def get_spectrogram(channel_data, time_resolution, timestep=None,
                    freq_max=None):
    if timestep is None:
        timestep = time_resolution/10
    freq_min = 0
    if freq_max is None:
        freq_max = 10000
    fft_data = copy.deepcopy(channel_data)
    time_vec = copy.deepcopy(channel_data.coordinate('Time')[0])
    new_time_step = time_vec[1]-time_vec[0]
    # print(new_time_step)
    # new_time_step = new_time_step[0]
    sampling = int(timestep/new_time_step)
    time_vec = time_vec[::sampling]
    time_vec2 = time_vec+time_resolution
    time_vec = time_vec[np.where(time_vec2<=np.max(channel_data.coordinate('Time')[0]))]
    time_vec2 = time_vec2[np.where(time_vec2<=np.max(channel_data.coordinate('Time')[0]))]
    for index in range(0,len(time_vec)):
        time_step = channel_data.slice_data(slicing={'Time': flap.Intervals(time_vec[index],time_vec2[index])})
        fourier = np.fft.fft(time_step.data-np.mean(time_step.data))
        n = time_step.data.size
        timestep = time_step.coordinate('Time')[0][1]-time_step.coordinate('Time')[0][0]
        freq = np.fft.fftfreq(n, d=timestep)
        fourier = fourier[np.where((freq>=freq_min)*(freq<freq_max))]
        freq = freq[np.where((freq>=freq_min)*(freq<freq_max))]
        if index == 0:
            fft_data.data = np.zeros([len(time_vec),len(fourier)])
        fourier = [x for _,x in sorted(zip(freq,fourier))]
        fsize = min(len(fourier), np.shape(fft_data.data)[1])
        fft_data.data[index,:fsize] = np.absolute(fourier[:fsize])
        fft_data.data[index,:fsize] = fft_data.data[index,:fsize]/np.max(fft_data.data[index,:fsize])
        if index%(len(time_vec)//20) == 0:
            print(f"{int(index/len(time_vec)*100)}%")
    fsize = min(len(freq), np.shape(fft_data.data)[1])
    freq=np.sort(freq)
    norm_data = abs(np.transpose(fft_data.data)/np.nanmax(channel_data.data, axis=0))**2
    return time_vec, freq, norm_data

def filter_dataobject(dataobject, flow, fhigh):
    time_vec = copy.deepcopy(dataobject.coordinate('Time')[0])
    time_step = time_vec[1]-time_vec[0]
    sos = signal.ellip(8, 1, 100, [flow, fhigh], btype="bandpass",\
                        fs=1/time_step, output='sos')
    new_dataobject = copy.copy(dataobject)
    new_dataobject.data[:] = signal.sosfilt(sos, new_dataobject.data[:])
    return new_dataobject

def plot_spectrogram(exp_id, channel, timeres=0.1, timestep=0.01,
                     timerange=None, flow=0, fhigh=None, plot_smoothening=0.5,
                     canvas=None):
    flap.config.read()
    
    print(f"Reading {channel} for {exp_id}")
    
    dataobject = flap.get_data('W7X_ABES', name=channel,
                               coordinates={'Time':timerange},
                               options={"Amplitude calibration": False},
                               exp_id=exp_id,
                               object_name='f{exp_id}/{channel}')
    try:
        spatcal = flap_w7x_abes.ShotSpatCal(exp_id)
        try:
            spatcal.generate_shotdata()
        except FileExistsError:
            pass
        spatcal.read()
        r = spatcal.data["Device R"][np.where(spatcal.data["Channel name"] == channel)][0]
        print("Spatial data obtained")
    except ValueError as e:
        print(f"Unable to read spatial calibration: {e}")
        r=None

    
    d_beam_on=flap.get_data('W7X_ABES',
                            exp_id=exp_id,
                            name='Chopper_time',
                            coordinates={'Time':timerange},
                            options={'State':{'Chop': 0, 'Defl': 0}},\
                            object_name='Beam_on',
                            )
    
    beam_on = dataobject.slice_data(slicing={'Sample':d_beam_on})
    beam_on = beam_on.slice_data(summing={'Rel. Sample in int(Sample)':'Mean'})
        
    if fhigh is not None and flow>0:
        print(f"Filtering for [{flow}Hz,{fhigh}Hz]")
        beam_on = filter_dataobject(beam_on, flow=flow, fhigh=fhigh)
    
    print(f"Obtaining spectrogram")
    flap_w7x_abes.regenerate_time_sample(beam_on)
    
    # dataobject.plot(axes='Time', options={'All':True}, plot_options={'marker': 'o'})
    # d_beam_on.plot(plot_type='scatter', axes=['Time', plt.ylim()[1]], options={'Force': True,'All': True})
    # beam_on.plot(axes='Time', options={'All':True}, plot_options={'marker': 'o'})

    time_vec, freq, norm_data = get_spectrogram(beam_on, timeres, timestep=timestep, freq_max=None)
    
    print(f"Plotting spectrogram")
    if canvas == None:
        fig = plt.figure(figsize=(4,4))
        fig.show()
        plt.pcolormesh(time_vec, freq/1000, np.sqrt((norm_data[:freq.shape[0]-1,:time_vec.shape[0]-1])**plot_smoothening))
        if r is not None:
            plt.title(f"{exp_id}/{channel} R={r}m\nTimeres: {timeres}s Freq:[{flow}Hz,{fhigh}Hz]")
        else:
            plt.title(f"{exp_id}/{channel}\nTimeres: {timeres}s Freq:[{flow}Hz,{fhigh}Hz]")
        plt.xlabel("Time [s]")
        plt.ylabel("Frequency [kHz]")
        plt.tight_layout()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()
    else:
        fig=canvas.fig
        ax = canvas.fig.add_subplot(1,1,1)
        ax.pcolormesh(time_vec, freq/1000, np.sqrt((norm_data[:freq.shape[0]-1,:time_vec.shape[0]-1])**plot_smoothening))
        if r is not None:
            canvas.fig.suptitle(f"{exp_id}/{channel} R={r}m\nTimeres: {timeres}s Freq:[{flow}Hz,{fhigh}Hz]")
        else:
            canvas.fig.suptitle(f"{exp_id}/{channel}\nTimeres: {timeres}s Freq:[{flow}Hz,{fhigh}Hz]")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Frequency [kHz]")
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

if __name__ == '__main__':
    # plot_shot("20241127.070", plot_rawdata=True, plot_powerspect=False, resample=0)

    exp_id = "20250227.026"
    timerange = None
    channel = "ABES-24"
    timeres = 0.1
    timestep = 0.01
    plot_smoothening = 0.5
    flow=10
    fhigh=10000
    plot_spectrogram(exp_id, channel, timeres=timeres, timestep=timestep, timerange=timerange,
                     flow=10, fhigh=10000, plot_smoothening=plot_smoothening)

    # globals()[sys.argv[1]](*sys.argv[2:])
