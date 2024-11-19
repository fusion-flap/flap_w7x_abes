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
        
        if plot_powerspect is True:
            print("in")
            fig = plt.figure()
            
            channels = np.arange(64)+1
            for channel in channels:
                channel_data = dataobject.slice_data(slicing={'ADC Channel': channel})
                plt.subplot(8,8,channel)
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
                ps = np.abs(np.fft.fft(channel_data.data))**2
        
                time_step = channel_data.coordinate('Time')[0][1]- channel_data.coordinate('Time')[0][0]
                freqs = np.fft.fftfreq(channel_data.data.size, time_step)
                idx = np.argsort(freqs)
                newfreq = freqs[idx]
                newps = ps[idx]
                newps = newps[np.where(newfreq>0)]
                newfreq = newfreq[np.where(newfreq>0)]
        
                plt.plot(newfreq, newps,
                         color = color, alpha=0.75, label=channel_data.coordinate('Detector type')[0][0]+str(channel)+"-"+channel_data.coordinate('Fibre ID')[0][0].split('-')[-1] )
                plt.ylim([np.min(newps)*0.8, np.max(newps)*1.2])
                try:
                    plt.yscale('log')
                    plt.xscale('log')
                except:
                    pass
                plt.ylabel('')
                plt.tick_params(axis='y', direction='in', pad=-22)
                plt.xlabel('')
                if channel < 56:
                    plt.xticks([])
                plt.legend()
                plt.show(block=False)
                plt.pause(0.001)
            #plt.get_current_fig_manager().full_screen_toggle()
            fig.text(0.5, 0.04, 'Frequency (Hz)', ha='center')
            fig.text(0.04, 0.5, 'Power spectrum', va='center', rotation='vertical')
            # fig.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
            fig.subplots_adjust(wspace=0,hspace=0)
            fig.show()
            plt.show(block=True)
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
        channel += 1
  
def smooth(y,box_pts):
    box=np.ones(box_pts)/box_pts
    y_smooth=np.convolve(y,box,mode='same')
    return np.array(y_smooth)

if __name__ == '__main__':
    # plot_shot("T20240806.018")
    plot_shot("20241105.053", plot_rawdata=True, plot_powerspect=False, resample=0)

    # globals()[sys.argv[1]](*sys.argv[2:])
