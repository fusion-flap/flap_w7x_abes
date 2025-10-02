#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 15:06:14 2025

@author: apdcam
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import copy
from scipy import interpolate

import flap
import flap_w7x_abes


voltage_difference = []
heating_current = []
beam_current_drop = []
timeval = []
extraction_time = []
extracted_charge = []
extractor_overcurrent = []
emitter_overcurrent = []

# dates = ["20240919",
#           "20240924", "20240925","20240926",
#           "20241001"] #III tested20240206
# dates = ["20241008", "20241009", "20241010",
#           "20241015", "20241016", "20241018"] #II, tested20240202
# dates = ["20241022", "20241023", "20241024",
#           "20241029",] #I
# dates = ["20241105", "20241106", "20241107",
#           "20241112", "20241113", "20241114",
#           "20241119", "20241120", "20241121",
#           "20241126", "20241127", "20241128",
#           "20241203", "20241204", "20241205"] #II, tested20240202 sanded down before, removed and  and on 2nd December 2024
# dates = ["20241209", "20241210", "20241211", "20241212"] # minimal sanding I
# dates = ["20250218", "20250219", "20250220",
#           "20250225", "20250226", "20250227",
#           "20250304", "20250305", "20250306",
#           "20250311", "20250312", "20250313",
#           "20250318", "20250319", "20250320"] #V
# dates = ["20250325", "20250326", "20250327",
#           "20250401", "20250402", "20250403",
#           "20250408", "20250409",
#           "20250415", "20250416", "20250417",
#           "20250423", "20250424",
#           "20250429"] #VI
dates = ["20250506", "20250507", "20250508",
          "20250513", "20250514", "20250515",
          "20250520", "20250521"] #IV

plt.figure()
dates_to_pop = []
for index, date in enumerate(dates):
    
    voltage_difference_date = []
    heating_current_date = []
    beam_current_drop_date = []
    emitter_overcurrent_date = []
    extracted_charge_date = []
    extractor_overcurrent_date = []

    
    try:
        beam_log = flap_w7x_abes.BORIMonitor(date=date, material="Na")
        
        neutralizer_shutter=None
        
        
        if ("Emitter overcurrent" not in list(beam_log.data.keys())) or ("Extractor overcurrent" not in list(beam_log.data.keys())):
            beam_log.get_overcurrent()
        
        em_voltage_der = np.abs(beam_log.data['HV Em Meas Voltage'].data[1:]-beam_log.data['HV Em Meas Voltage'].data[:-1])
        ex_voltage_der = np.abs(beam_log.data['HV Ex Meas Voltage'].data[1:]-beam_log.data['HV Ex Meas Voltage'].data[:-1])
        
        if neutralizer_shutter is None:
            rel_time = np.where((em_voltage_der<0.1)*\
                                (ex_voltage_der<0.1)*\
                                (beam_log.data['HV Em Meas Voltage'].data[1:]>0.1)*\
                                (beam_log.data['HV Ex Meas Voltage'].data[1:]<beam_log.data['HV Em Meas Voltage'].data[1:])*\
                                (beam_log.data['HV Em Meas Voltage'].data[1:]-beam_log.data['HV Ex Meas Voltage'].data[1:]>3.5))
        else:
            rel_time = np.where((em_voltage_der<0.1)*\
                                (ex_voltage_der<0.1)*\
                                (abs(beam_log.data['Neut Shut Closed'].data[1:] - neutralizer_shutter)<0.5)*\
                                (beam_log.data['HV Em Meas Voltage'].data[1:]>0.1)*\
                                (beam_log.data['HV Ex Meas Voltage'].data[1:]<beam_log.data['HV Em Meas Voltage'].data[1:])*\
                                (beam_log.data['HV Em Meas Voltage'].data[1:]-beam_log.data['HV Ex Meas Voltage'].data[1:]>4.))
        if len(rel_time[0]) == 0:
            dates_to_pop += [date]
        else:
            plt.scatter(beam_log.data["Emit Current A"].data[rel_time]**2*0.08, beam_log.data["Emitter overcurrent"].data[rel_time], label=date, color=(index/len(dates), np.abs(2*index/len(dates)-1), 1-index/len(dates)/2), alpha=0.15)
            # plt.scatter(beam_log.data["Emit Current A"].data[rel_time], beam_log.data["Emitter overcurrent"].data[rel_time], label=date, color=(index/len(dates), np.abs(2*index/len(dates)-1), 1-index/len(dates)/2), alpha=0.15)
            plt.show()
            plt.pause(1)
    
            #grouping neighbouring time points
            group_times = []
            for rel_time_ind in rel_time[0]:
                try:
                    if rel_time_ind - group_times[-1][-1] > 5:
                        group_times += [[rel_time_ind]]
                    else:
                        group_times[-1] += [rel_time_ind]
                except IndexError:
                    group_times += [[rel_time_ind]]
    
            for shot in group_times:
                # #only look at where the beam already stabilized
                # emitter_overcurrent = beam_log.data["Emitter overcurrent"].data[1:][shot]
                # startpoint = np.min(np.where(np.abs(emitter_overcurrent-emitter_overcurrent[-1])/emitter_overcurrent[-1]<0.1)[0])
    
                # vd = np.mean(beam_log.data['HV Em Meas Voltage'].data[1:][shot] - beam_log.data['HV Ex Meas Voltage'].data[1:][shot[startpoint:]])
                # hc = np.mean(beam_log.data['Emit Current A'].data[1:][shot[startpoint:]])
                # bc = np.mean(beam_log.data["Emitter overcurrent"].data[1:][shot[startpoint:]])
                
                vd = beam_log.data['HV Em Meas Voltage'].data[1:][shot] - beam_log.data['HV Ex Meas Voltage'].data[1:][shot][-1]
                hc = beam_log.data['Emit Current A'].data[1:][shot][-1]
                bcd = (beam_log.data["Emitter overcurrent"].data[1:][shot][0]-beam_log.data["Emitter overcurrent"].data[1:][shot][-1])/beam_log.data["Emitter overcurrent"].data[1:][shot][-1]         
                bcd = (beam_log.data["Emitter overcurrent"].data[1:][shot][0]-beam_log.data["Emitter overcurrent"].data[1:][shot][-1])/beam_log.data["Emitter overcurrent"].data[1:][shot][-1]
                bc = np.mean(beam_log.data["Emitter overcurrent"].data[1:][shot])
                
                voltage_difference_date += [vd]
                heating_current_date += [hc]
                beam_current_drop_date += [bcd]
                emitter_overcurrent_date += [bc]
                extractor_overcurrent_date  += [np.mean(beam_log.data["Extractor overcurrent"].data[shot])]
                timeval += [beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][-1]]
                try:
                    extraction_time += [extraction_time[-1]+beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][-1]-beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][0]]
                except IndexError:
                    extraction_time += [beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][-1]-beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][0]]
                try:
                    extracted_charge_date += [extracted_charge_date[-1]+np.mean(beam_log.data["Emitter overcurrent"].data[shot])*(beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][-1]-beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][0])]
                except IndexError:           
                    extracted_charge_date += [np.mean(beam_log.data["Emitter overcurrent"].data[shot])*(beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][-1]-beam_log.data["Emitter overcurrent"].coordinate("Time")[0][shot][0])]
            voltage_difference += [voltage_difference_date]
            heating_current += [heating_current_date]
            beam_current_drop += [beam_current_drop_date]
            extracted_charge += [extracted_charge_date]
            extractor_overcurrent += [extractor_overcurrent_date]
            emitter_overcurrent += [emitter_overcurrent_date]
    except ValueError:
        dates_to_pop += [date]

for date in dates_to_pop:
    dates.remove(date)

plt.legend(ncol=2)
# plt.xlabel("Heating current [A]")
plt.xlabel("Emitter heating power [W]")

plt.ylabel("Extracted current [mA]")
plt.tight_layout()
    
    
    # plt.scatter(-1/heating_current_date, beam_current_drop_date, alpha=0.2, label=date)

plt.figure()
startindex = 0
lastextraction = 0
for i, day in enumerate(dates):
    try:
        plt.scatter(extraction_time[startindex:startindex+len(extracted_charge[i])], (np.asarray(extracted_charge[i])+lastextraction)/3600,
                    alpha=0.8, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2), label=dates[i])
        lastextraction += extracted_charge[i][-1]
        startindex += len(extracted_charge[i])
    except IndexError:
        pass
plt.xlabel("Extraction time [s]")
plt.ylabel("Extracted charge [mAh]")
plt.legend(ncol=2)


import matplotlib as mpl
mpl.rcParams['text.usetex'] = False
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
fig = plt.figure(figsize=[8,4])
ax = plt.subplot(2,1,1)
startindex = 0
lastextraction = 0
for i, day in enumerate(dates):
    try:
        datatoplot = beam_current_drop[i]
        plt.scatter((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    np.asarray(beam_current_drop[i]),
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2), label=dates[i])
        plt.plot((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    np.asarray(beam_current_drop[i]),
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2))
        lastextraction += extracted_charge[i][-1]
        startindex += len(beam_current_drop[i])
    except IndexError:
        pass
# plt.xlabel("Extracted charge [mAh]")
plt.ylim(-0.1,0.6)
handles, labels=plt.gca().get_legend_handles_labels()
plt.ylabel("Relative beam current\n drop during shot [1]")
box=ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

# plt.figure()
ax = plt.subplot(2,1,2)
startindex = 0
lastextraction = 0
for i, day in enumerate(dates):
    try:
        # current_per_shot = np.concatenate([[extracted_charge[i][0]], np.asarray(extracted_charge[i])[1:]-np.asarray(extracted_charge[i])[:-1]])
        datatoplot = beam_current_drop[i]
        plt.scatter((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    -np.asarray(extractor_overcurrent[i])/np.asarray(emitter_overcurrent[i]),
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2), label=dates[i])
        plt.plot((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    -np.asarray(extractor_overcurrent[i])/np.asarray(emitter_overcurrent[i]),
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2))
        lastextraction += extracted_charge[i][-1]
        startindex += len(extractor_overcurrent[i])
    except IndexError:
        pass
plt.xlabel("Extracted charge [mAh]")
plt.ylabel(r"$\frac{\mathrm{Extractor\, overcurrent}}{\mathrm{Extracted\, beam\, current}}$", fontsize=14)
plt.yscale("log")
# plt.legend(ncol=2)
plt.suptitle("Ion source characteristics")
# lines_labels = [ax.get_legend_handles_labels for ax in fig.axes]
# lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
box=ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
plt.legend(handles, labels, bbox_to_anchor=(1.3,1.2), loc="center right")
# plt.tight_layout()

plt.figure()
startindex = 0
lastextraction = 0
for i, day in enumerate(dates):
    try:
        datatoplot = beam_current_drop[i]
        plt.scatter((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    np.asarray(extractor_overcurrent[i]),
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2), label=dates[i])
        plt.plot((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    np.asarray(extractor_overcurrent[i]),
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2))
        lastextraction += extracted_charge[i][-1]
        startindex += len(extractor_overcurrent[i])
    except IndexError:
        pass
plt.xlabel("Extracted charge [mAh]")
plt.ylabel("Extractor overcurrent [mA]")
# plt.yscale("log")
plt.legend(ncol=2)
plt.tight_layout()


plt.figure()
startindex = 0
lastextraction = 0
for i, day in enumerate(dates):
    try:
        # current_per_shot = np.concatenate([[extracted_charge[i][0]], np.asarray(extracted_charge[i])[1:]-np.asarray(extracted_charge[i])[:-1]])
        datatoplot = beam_current_drop[i]
        # heating_power = [np.asarray(heating_current_day)**2 for heating_current_day in heating_current]
        
        plt.scatter((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    np.asarray(heating_current[i])**2,
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2), label=dates[i],
                    marker="x")
        plt.plot((np.asarray(extracted_charge[i])+lastextraction)/3600,
                    np.asarray(heating_current[i])**2,
                    alpha=0.5, color=(i/len(dates), np.abs(2*i/len(dates)-1), 1-i/len(dates)/2))
        lastextraction += extracted_charge[i][-1]
        startindex += len(extractor_overcurrent[i])
    except IndexError:
        pass
plt.xlabel("Extracted charge [mAh]")
plt.ylabel("Heating current **2")
plt.yscale("log")
plt.legend(ncol=2)
beam_log.data["Emit Current A"].data[rel_time]**2*0.08