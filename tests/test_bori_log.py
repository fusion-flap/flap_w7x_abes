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

import flap
import flap_w7x_abes

def overview_plot(date, reference_days=[], last_minutes=20):
     
    beam_log = flap_w7x_abes.bori_log.BORIMonitor(date=date)
    
    figure = plt.figure(figsize=[18,10], constrained_layout=True)
    gs = gridspec.GridSpec(4, 6, figure=figure, hspace=0.5,wspace=0.5)
    logtext = ""
    
    em_resistance, ex_resistance = beam_log.get_load_resistance()
    if not (abs(em_resistance-82)<3 and abs(ex_resistance-72)<3):
        logtext += f'Calculated resistances R_emit={beam_log.em_resistance}MOhm  R_ext={beam_log.ex_resistance}MOhm\n'
        beam_log.em_resistance = 82
        beam_log.ex_resistance = 72
     
    beam_log_lastminutes = beam_log.slice_time(last_minutes)

    plt.suptitle(f'Date: {beam_log.date}, R_emit={beam_log.em_resistance}MOhm  R_ext={beam_log.ex_resistance}MOhm',fontsize=14)

    beam_log_lastminutes.get_child_langmuir()

    ax = plt.subplot(gs[0,0:2])
    beam_log_lastminutes.data['Emit Current A'].plot()
    plt.title('Emitter heating current')
    plt.xlabel("")
    
    plt.subplot(gs[0,2:4],sharex=ax)
    beam_log_lastminutes.data['HV Em Meas Voltage'].plot(plot_options={"label":"Emitter"})
    beam_log_lastminutes.data['HV Ex Meas Voltage'].plot(plot_options={"label":"Extractor"})
    plt.title('Accelerator voltages')
    plt.legend()
    plt.xlabel("")
    
    plt.subplot(gs[1,0:2],sharex=ax)
    legend = []
    beam_log_lastminutes.data['TC Oven Bottom'].plot(plot_options={"label":"Bottom"})
    beam_log_lastminutes.data['TC Oven Top'].plot(plot_options={"label":"Top"})
    beam_log_lastminutes.data['TC Torus Side Cone'].plot(plot_options={"label":"Torus side"})
    beam_log_lastminutes.data['TC Emit Side Cone'].plot(plot_options={"label":"Emitter side"})
    plt.title('Neutralizer temperatures')
    plt.legend()
    plt.xlabel("")
    
    plt.subplot(gs[1,2:4],sharex=ax)
    beam_log_lastminutes.data['HV Em Meas Current'].plot(plot_options={"label":"Emitter"})
    beam_log_lastminutes.data['HV Ex Meas Current'].plot(plot_options={"label":"Extractor"})
    plt.title('HV PS currents')
    plt.legend()
    plt.xlabel("")
    
    plt.subplot(gs[2,0:2],sharex=ax)
    beam_log_lastminutes.data['VG HighVac1'].plot(plot_options={"label":"HighVac1"})
    beam_log_lastminutes.data['VG HighVac2'].plot(plot_options={"label":"HighVac2"})
    plt.xlabel("")
    plt.yscale('log')
    plt.title('Vacuum pressure')
    plt.legend()
    
    plt.subplot(gs[2,2:4],sharex=ax)
    data1 = beam_log_lastminutes.data['Emitter overcurrent']
    data2 = beam_log_lastminutes.data['Extractor overcurrent']
    data1.data += data2.data
    data1.plot()
    plt.xlabel("")
    plt.title('Beam current')


    plt.subplot(gs[3,0:2],sharex=ax)
    beam_log_lastminutes.data['Emitter overcurrent'].plot(plot_options={"label":"Emitter"})
    beam_log_lastminutes.data['Extractor overcurrent'].plot(plot_options={"label":"Extractor"})
    plt.legend()
    plt.title('Extra currents')
    
    plt.subplot(gs[3,2:4],sharex=ax)
    fc2_rel_index = np.where((beam_log_lastminutes.data["FC2 in"].data>0)*\
                             (beam_log_lastminutes.data["Neut Shut Closed"].data>0)*\
                             (beam_log_lastminutes.data["FC1 in"].data == 0))
    fc2_data = beam_log_lastminutes.data['FC2 Resistor Current mA'].data[fc2_rel_index]
    fc2_reltime = beam_log_lastminutes.data['FC2 Resistor Current mA'].get_coordinate_object("Time").values[fc2_rel_index]
    plt.plot(fc2_reltime, fc2_data, label="FC2")
    fc1_rel_index = np.where((beam_log_lastminutes.data["FC1 in"].data > 0))
    fc1_data = beam_log_lastminutes.data['FC1 Resistor Current mA'].data[fc1_rel_index]
    fc1_reltime = beam_log_lastminutes.data['FC1 Resistor Current mA'].get_coordinate_object("Time").values[fc1_rel_index]
    plt.plot(fc1_reltime, fc1_data, label="FC1")
    plt.xlabel(f"Time {beam_log_lastminutes.data['FC1 Resistor Current mA'].get_coordinate_object('Time').unit.unit}")
    plt.ylabel(f"Time {beam_log_lastminutes.data['FC1 Resistor Current mA'].data_unit.unit}")
    plt.title("Faraday cup currents")
    plt.legend()
    
    plt.subplot(gs[0:2,4:6])
    beam_log.plot_child_langmuir(plotcolor="tab:blue", neutralizer_shutter=0, newfigure=False)
    beam_log.plot_child_langmuir(plotcolor="tab:orange", neutralizer_shutter=1, newfigure=False)
    beam_log_lastminutes.plot_child_langmuir(plotcolor="tab:green", neutralizer_shutter=0, newfigure=False, label="last 20min")
    beam_log_lastminutes.plot_child_langmuir(plotcolor="tab:red", neutralizer_shutter=0, newfigure=False, label="last 20min")

    for reference_day in reference_days:
        beam_log_ref = flap_w7x_abes.bori_log.BORIMonitor(date=reference_day)
        beam_log_ref.plot_child_langmuir(plotcolor="tab:pink", newfigure=False)
    plt.title('Child-Langmuir')
    
    plt.subplot(gs[2:4,4:6])
    beam_log.plot_neutralizer(plotcolor="tab:blue", newfigure=False)
    for reference_day in reference_days:
        beam_log_ref.plot_neutralizer(plotcolor="tab:pink", newfigure=False)
    
    plt.text(0, -0.15, logtext, color="tab:red")
    gs.tight_layout(figure)

if __name__ == "__main__":
    # overview_plot("20240924", last_minutes=20, reference_days=["20240923"])
    #This can be run from command line by running python test_bori_log.py 20240924
    datetoget=sys.argv[1]
    overview_plot(datetoget,last_minutes=20, reference_days=["20240923"])
    plt.savefig("/home/apdcam/Measurement/beam_20min.png")