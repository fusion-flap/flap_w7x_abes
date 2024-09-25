#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:57:38 2024

@author: apdcam
"""
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import copy

import flap
import flap_w7x_abes

def overview_plot(date, reference_days=[], last_minutes=None):
     
    beam_log = flap_w7x_abes.bori_log.BORIMonitor(date=date)
    
    figure = plt.figure(figsize=[18,10], constrained_layout=True)
    gs = gridspec.GridSpec(8, 6, figure=figure, hspace=2,wspace=0.5)
    logtext = ""
    
    em_resistance, ex_resistance = beam_log.get_load_resistance()
    if not (abs(em_resistance-82)<3 and abs(ex_resistance-72)<3):
        logtext += f'CALCULATED RESISTANCE R_emit={int(beam_log.em_resistance)}MOhm  R_ext={int(beam_log.ex_resistance)}MOhm\n'
        beam_log.em_resistance = 82
        beam_log.ex_resistance = 72

    if last_minutes is not None:
        beam_log_lastminutes = beam_log.slice_time(last_minutes)
        beam_log.get_child_langmuir
    else:
        beam_log_lastminutes = beam_log

    beam_log_lastminutes.get_child_langmuir()
    beam_log.get_extracted_charge()
    #convert time axis to minutes
    for index, data in enumerate(beam_log.data.keys()):
        if index == 0:
            starttime = beam_log.data[data].get_coordinate_object("Time").values[0]
        time_vect = beam_log.data[data].get_coordinate_object("Time")
        time_vect.values -= np.max(time_vect.values)
        time_vect.values = time_vect.values/60
        time_vect.unit.unit = "min"
        beam_log.data[data].del_coordinate("Time")
        beam_log.data[data].add_coordinate_object(time_vect)
    

    if last_minutes is not None:
        beam_log_lastminutes.get_extracted_charge()
    
        #convert time axis to minutes
        for data in beam_log_lastminutes.data.keys():
            time_vect = beam_log_lastminutes.data[data].get_coordinate_object("Time")
            time_vect.values -= np.max(time_vect.values)
            time_vect.values = time_vect.values/60
            time_vect.unit.unit = "min"
            beam_log_lastminutes.data[data].del_coordinate("Time")
            beam_log_lastminutes.data[data].add_coordinate_object(time_vect)

    plt.suptitle(f'Date: {date}, R_emit={beam_log.em_resistance}MOhm  R_ext={beam_log.ex_resistance}MOhm',fontsize=14)

    ax = plt.subplot(gs[0:2,0:2])
    beam_log_lastminutes.data['Emit Current A'].plot()
    plt.title('Emitter heating current')
    
    plt.subplot(gs[0:2,2:4],sharex=ax)
    beam_log_lastminutes.data['HV Em Meas Voltage'].plot(plot_options={"label":"Emitter"})
    beam_log_lastminutes.data['HV Ex Meas Voltage'].plot(plot_options={"label":"Extractor"})
    plt.title('Accelerator voltages')
    plt.legend()
    
    plt.subplot(gs[2:4,0:2],sharex=ax)
    beam_log_lastminutes.data['TC Oven Bottom'].plot(plot_options={"label":"Bottom", "color":"tab:red"})
    beam_log_lastminutes.data['TC Oven Top'].plot(plot_options={"label":"Top", "color":"tab:purple"})
    beam_log_lastminutes.data['TC Torus Side Cone'].plot(plot_options={"label":"Torus side", "color":"tab:pink"})
    beam_log_lastminutes.data['TC Emit Side Cone'].plot(plot_options={"label":"Emitter side", "color":"tab:brown"})
    plt.title('Neutralizer temperatures')
    plt.legend()
    
    plt.subplot(gs[2:4,2:4],sharex=ax)
    beam_log_lastminutes.data['HV Em Meas Current'].plot(plot_options={"label":"Emitter"})
    beam_log_lastminutes.data['HV Ex Meas Current'].plot(plot_options={"label":"Extractor"})
    plt.title('HV PS currents')
    plt.legend()
    
    plt.subplot(gs[4:6,0:2],sharex=ax)
    beam_log_lastminutes.data['VG HighVac1'].plot(plot_options={"label":"HighVac1", "color":"tab:gray"})
    beam_log_lastminutes.data['VG HighVac2'].plot(plot_options={"label":"HighVac2", "color":"tab:olive"})
    plt.yscale('log')
    plt.title('Vacuum pressure')
    plt.legend()
    
    plt.subplot(gs[4:6,2:4],sharex=ax)
    beam_log_lastminutes.data['Emitter overcurrent'].plot(plot_options={"label":"Emitter"})
    beam_log_lastminutes.data['Extractor overcurrent'].plot(plot_options={"label":"Extractor"})
    plt.legend()
    plt.title('Extra currents')

    plt.subplot(gs[6:8,0:2],sharex=ax)
    fc2_rel_index = np.where((beam_log_lastminutes.data["FC2 in"].data>0)*\
                             (beam_log_lastminutes.data["Neut Shut Closed"].data>0)*\
                             (beam_log_lastminutes.data["FC1 in"].data == 0))
    fc2_data = beam_log_lastminutes.data['FC2 Resistor Current mA'].data[fc2_rel_index]
    fc2_reltime = beam_log_lastminutes.data['FC2 Resistor Current mA'].get_coordinate_object("Time").values[fc2_rel_index]
    plt.plot(fc2_reltime, fc2_data, label="FC2", color="tab:gray")
    fc1_rel_index = np.where((beam_log_lastminutes.data["FC1 in"].data > 0))
    fc1_data = beam_log_lastminutes.data['FC1 Resistor Current mA'].data[fc1_rel_index]
    fc1_reltime = beam_log_lastminutes.data['FC1 Resistor Current mA'].get_coordinate_object("Time").values[fc1_rel_index]
    plt.plot(fc1_reltime, fc1_data, label="FC1", color="tab:olive")
    plt.xlabel(f"Time [{beam_log_lastminutes.data['FC1 Resistor Current mA'].get_coordinate_object('Time').unit.unit}]")
    plt.ylabel(f"Current [{beam_log_lastminutes.data['FC1 Resistor Current mA'].data_unit.unit}]")
    plt.title("Faraday cup currents")
    plt.legend()

    plt.subplot(gs[6:8,2:4],sharex=ax)
    data1 = copy.deepcopy(beam_log_lastminutes.data['Emitter overcurrent'])
    data2 = beam_log_lastminutes.data['Extractor overcurrent']
    data1.data += data2.data
    gate_valve_open = np.where(beam_log_lastminutes.data['PB Feedback-Open'].data>0.5)
    datatemp = copy.deepcopy(data1)
    datatemp.data[gate_valve_open] = np.nan
    datatemp.plot(plot_options={"label":"Gate Valve Closed", "color":"black"})
    gate_valve_closed = np.where(beam_log_lastminutes.data['PB Feedback-Open'].data<0.5)
    data1.data[gate_valve_closed] = np.nan
    data1.plot(plot_options={"label":"Gate Valve Open", "color":"tab:green"})
    plt.legend()
    plt.title('Beam current')
    
    plt.subplot(gs[0:4,4:6])
    if last_minutes is not None:
        beam_log.plot_child_langmuir(plotcolor="tab:cyan", newfigure=False)

    beam_log_lastminutes.plot_child_langmuir(plotcolor="tab:blue", neutralizer_shutter=0, newfigure=False,
                                             label=f"Last {last_minutes}min - neut. shut. open", alpha=1)
    beam_log_lastminutes.plot_child_langmuir(plotcolor="tab:red", neutralizer_shutter=1, newfigure=False,
                                             label=f"Last {last_minutes}min - neut. shut. closed", alpha=1)


    for reference_day in reference_days:
        beam_log_ref = flap_w7x_abes.bori_log.BORIMonitor(date=reference_day)
        beam_log_ref.plot_child_langmuir(plotcolor="tab:gray", newfigure=False)
    plt.title('Child-Langmuir')
    
    plt.subplot(gs[4:6,4:6], sharex=ax)
    beam_log.data["Extracted charge"].plot(plot_options={"label":f"Cumulated {date}", "color":"tab:cyan"})
    plt.plot(beam_log_lastminutes.data["Extracted charge"].get_coordinate_object("Time").values,
             beam_log_lastminutes.data["Extracted charge"].data,
             color="tab:blue", label = f"Last {last_minutes}min")
    plt.legend()
    plt.title("Extracted charge from emitter")
    plt.xlim([np.min(beam_log_lastminutes.data["Extracted charge"].get_coordinate_object("Time").values),
              np.max(beam_log_lastminutes.data["Extracted charge"].get_coordinate_object("Time").values)])
    
    
    plt.subplot(gs[6:8,4:6])
    beam_log.plot_neutralizer(plotcolor="tab:cyan", newfigure=False)
    for reference_day in reference_days:
        beam_log_ref.plot_neutralizer(plotcolor="tab:gray", newfigure=False)
    plt.title('Neutralization')
    
    plt.text(0, -0.15, logtext, color="tab:red")

if __name__ == "__main__":
    # overview_plot("20240924", last_minutes=20, reference_days=["20240923"])
    # overview_plot("20240924", reference_days=["20240923"])
    #This can be run from command line by running python test_bori_log.py 20240924
    datetoget=sys.argv[1]
    overview_plot(datetoget,last_minutes=20, reference_days=["20240923"])
    plt.savefig("/home/apdcam/Measurement/borilog_test.png", dpi=150)