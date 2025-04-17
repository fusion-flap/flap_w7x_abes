# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:23:49 2018

@author: Miklos Vecsei, Fusion Plasma Physics Department, Centre for Energy Research

Spatial calibration of the for W7-X alkali BES diagnostic module for flap
"""

import os
import sys
import warnings
import datetime
from scipy import fft
import h5py
import numpy as np
import flap

def read_fibre_config(exp_id=None, year=None):
    #finding the configuration file
    if exp_id != None:
        currdir = os.path.dirname(os.path.abspath(__file__))
        if int(exp_id[:4]) <= 2023:
            year = 2021 
        elif int(exp_id[:4]) <= 2024:
             year = 2021
        else:
            year = 2025
    config_file = os.path.join(currdir, 'spatcal', str(year), 'cxrs_fiber_patchconfig.dat')
    
    # reading the data
    with open(config_file, "r") as fibre_config:
        lines = fibre_config.readlines()
    # reading the patching from Domoknos's patch box
    patch_fibers_oc = {}
    # patch_fibers_oc_na = {}
    for line in lines[1:12]:
        for connection_point in line.split("\t")[1:-1]:
            patch_fibers_oc[connection_point.split("/")[1]] = connection_point.split("/")[0]
    dead_index = 0
    for connection_point in lines[13].split("\t")[1:]:
        dead_index += 1
        curr_id = connection_point.split('/')[1].split('\n')[0]
        patch_fibers_oc[f"BR1.{dead_index}"] = connection_point.split("/")[0]
    
    # temp = {}
    # for key in sorted(patch_fibers_oc.keys()):
    #     temp[patch_fibers_oc[key].zfill(2)] = key
    # relkeys = [key for key in sorted(temp.keys()) if key != "0S"]
    # index = 1
    # allkeys = np.flip(np.asarray(sorted(relkeys)))
    # for key in range(len(allkeys)//2):
    #     print(f"{index} {allkeys[2*key]} {temp[allkeys[2*key]]} \t {allkeys[2*key+1]} {temp[allkeys[2*key+1]]}")
    #     index += 1

    # reading the second configuration of the second patch panel
    patch_fibers_spectf = {}
    for line in lines[16:-1]:
        for connection_point in line.split("\t")[1:-1]:
            patch_fibers_spectf[connection_point.split("/")[1]] = connection_point.split("/")[0]
        connection_point = line.split("\t")[-1]
        patch_fibers_spectf[connection_point.split("/")[1].split("\n")[0]] = connection_point.split("/")[0]
    NA = [dead.split("\n")[0] for dead in lines[-1].split("\t")[1:]]
    BR2_index = 0
    for connection_point in NA:
        patch_fibers_spectf[connection_point] = f"BR2.{BR2_index}"
        BR2_index += 1

    # merging the two
    patch_oc_spectf = {}
    for key in patch_fibers_oc.keys():
        patch_oc_spectf[patch_fibers_oc[key]] = 'N.A.'     
    BR2_index = 0
    for fiber in patch_fibers_spectf.keys():
        try:
            if fiber in NA:
                BR2_index += 1
                patch_oc_spectf[patch_fibers_oc[fiber]] = f"BR2.{BR2_index}"
            else:
                patch_oc_spectf[patch_fibers_oc[fiber]] = patch_fibers_spectf[fiber]
        except KeyError as e:
            pass
    br1 = [patch_fibers_oc[key] for key in patch_fibers_oc.keys() if "BR" in key]
    BR1_index = 0
    for broken in br1:
        patch_oc_spectf[broken] = f"BR1.{BR1_index}"
        BR1_index += 1
    
    return patch_oc_spectf, patch_fibers_spectf, patch_fibers_oc

def plot_fibre_config(patch_oc_spectf):
    xstep = 3
    ystep = 1
    
    locations = {}
    for oc in np.arange(90)+1:
        locations[str(oc)] = [(oc-1)%2-0.5, (oc-1)//2]

    
    naming = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    for side_oc_id in range(len(naming)):
        starter_loc_x = [-3, 2]
        starter_loc_y = [9, 16.5, 23.5,30.5, 38]
        starter_loc = np.array([starter_loc_x[side_oc_id%2], starter_loc_y[side_oc_id//2]])
        for fibers in np.arange(4)+1:
            if (fibers-1)//2 == 0:
                locations[naming[side_oc_id]+str(fibers)] = [starter_loc[0]+(fibers-1)%2, starter_loc[1]+(fibers-1)//2]
            else:
                locations[naming[side_oc_id]+str(fibers)] = [starter_loc[0]+(fibers)%2, starter_loc[1]+(fibers-1)//2]

    from matplotlib import pyplot as plt
    plt.figure(figsize=[12,9])
    color = {'A':'tab:blue', 'H':'tab:green', 'HF': 'red', 'Z':'chocolate', 'NE':'tab:pink', 'BES':'gold'}
    pointid=0
    x = [xstep*locations[channel][0] for channel in locations.keys()]
    y = [ystep*locations[channel][1] for channel in locations.keys()]
    plt.clf()
    sc = plt.scatter(x, y, s=0, facecolors="none", edgecolors="black")

    for channel in locations.keys():
        try:
            plt.text(xstep*locations[channel][0]-0.5, ystep*locations[channel][1], f"{channel}/{patch_oc_spectf[str(channel)].split('.')[1]}",
                     c=color[patch_oc_spectf[str(channel)].split('.')[0]])
        except Exception as e:
            try:
                if "BR" in patch_oc_spectf[str(channel)].split('.')[0]:
                    t = plt.text(xstep*locations[channel][0]-0.5, ystep*locations[channel][1], f"{channel}/{patch_oc_spectf[str(channel)].split('.')[0]}", c="white", backgroundcolor="black", weight="bold")
                    t.set_bbox(dict(facecolor='black', alpha=0.25))
                else:
                    plt.text(xstep*locations[channel][0]-0.5, ystep*locations[channel][1], f"{channel}/N.A.", c="black")
            except Exception as e:
                    t = plt.text(xstep*locations[channel][0]-0.5, ystep*locations[channel][1], f"{channel}/BR1", c="white", backgroundcolor="black", weight="bold")
                    t.set_bbox(dict(facecolor='black', alpha=0.25))

        pointid += 1
    plt.text(0,-2*ystep,"INBOARD", weight="bold")
    plt.text(0,46*ystep,"OUTBOARD", weight="bold")
    plt.axis('equal')
    plt.xlim([-5,5])
    plt.axis("off")
    plt.text(5, 0, 'Optical channel / Alkali fiber number', c=color['A'], weight="bold")
    plt.text(5, 1, 'Optical channel / Helium fiber number', c=color['H'], weight="bold")
    plt.text(5, 2, 'Optical channel / Helium Filterscope fiber number', c=color['HF'], weight="bold")
    plt.text(5, 3, 'Optical channel / Zeff fiber number', c=color['Z'], weight="bold")
    plt.text(5, 4, 'Optical channel / Neon fiber number', c=color['NE'], weight="bold")
    plt.text(5, 5, 'Optical channel / BES fiber number', c=color['BES'], weight="bold")
    t = plt.text(5, 6, 'BR1 - detached or broken at port / BR2 - broken form port to patchbox', c="white", weight="bold")
    t.set_bbox(dict(facecolor='black', alpha=0.25))
    plt.tight_layout()
    # plt.gca().invert_yaxis()
    plt.ylim([46*ystep, -2*ystep])

def plot_patchpanel_spectrometer_config(spect_config):
    from matplotlib import pyplot as plt
    color = {'A':'tab:blue', 'H':'tab:green', 'HF': 'red', 'Z':'chocolate', 'NE':'tab:pink', 'BES':'gold', 'BR2':'white'}
    plt.figure(figsize=[15,6])
    for key in spect_config.keys():
        stringkey = str(key)
        if "N.A" not in stringkey:
            panel = int(stringkey.split(".")[0])
            location = int(stringkey.split(".")[1])
            plt.scatter(panel*11+(location-1)%8*1.2-10, -location//8, alpha=0)
            t = plt.text(panel*11+(location-1)%8*1.2-10, -location//8, str(location)+"/"+spect_config[key],
                         color=color[spect_config[key].split('.')[0]], weight="bold")
            if "BR" in spect_config[key]:
                    t.set_bbox(dict(facecolor='black', alpha=0.25))

    plt.title("Patchpanel locations/Spectrometer channels")
    plt.axis('equal')
    plt.tight_layout()
    plt.axis('off')
    
def plot_patchpanel_optical_channel_config(spect_config, patchp_config):
    from matplotlib import pyplot as plt
    color = {'A':'tab:blue', 'H':'tab:green', 'HF': 'red', 'Z':'chocolate', 'NE':'tab:pink', 'BES':'gold', 'BR2':'white'}
    plt.figure(figsize=[12,6])
    for key in patchp_config.keys():
        stringkey = str(key)
        if ("S" not in stringkey) and ("BR1" not in stringkey):
            panel = int(str(key).split(".")[0])
            location = int(str(key).split(".")[1])
            plt.scatter(panel*10+(location-1)%8, -location//8, alpha=0)
            # plt.text(panel*10+(location-1)%8, -location//8, str(location)+"/"+spect_config[key],
            #          color=color[spect_config[key].split('.')[0]])
            if key in spect_config.keys():
                if spect_config[key].split('.')[0] != "BR2":
                    plt.text(panel*10+(location-1)%8, -location//8, str(location)+"/"+patchp_config[key],
                             color=color[spect_config[key].split('.')[0]], weight="bold")
                else:
                    t = plt.text(panel*10+(location-1)%8, -location//8, str(location)+"/BR2",
                             color=color[spect_config[key].split('.')[0]], weight="bold")
                    t.set_bbox(dict(facecolor='black', alpha=0.25))

            else:
                plt.text(panel*10+(location-1)%8, -location//8, str(location)+"/"+patchp_config[key],
                          color="black", weight="bold")
    plt.title("Patchpanel locations/Optical channels")
    plt.axis('equal')
    plt.tight_layout()
    plt.axis('off')

def spect_calib_cxrs_op21():
    dataset = dict()
    dataset["grating_check"] = [[datetime.datetime(2023, 4, 27, 10, 43, 00),
                                 datetime.datetime(2023, 4, 27, 10, 56, 00)],
                                ["P04-grat1-540","P04-grat1-541","P04-grat0-540","P04-grat0-541","P04-grat2-540","P04-grat2-541"]]
    dataset["540_2400"] = [[datetime.datetime(2023, 4, 26, 13, 57, 00),
                            datetime.datetime(2023, 4, 26, 14, 13, 00)],
                           ["P01-1200-540","P02-1200-540","P03-1200-540","P04-1200-540","P05-1200-540","P06-1200-540"]]
    dataset["540_1800"] = [[datetime.datetime(2023, 4, 26, 14, 30, 00),
                            datetime.datetime(2023, 4, 26, 14, 54, 00)],
                           ["P01-2400-540","P02-2400-540","P03-2400-540","P04-2400-540","P05-2400-540","P06-2400-540"]]
    dataset["540_1200"] = [[datetime.datetime(2023, 4, 26, 15, 00, 00),
                            datetime.datetime(2023, 4, 26, 15, 18, 00)],
                           ["P01-1800-540","P02-1800-540","P03-1800-540","P04-1800-540","P06-1800-540","P05-1800-540"]]

    #checking the gratings for grat0/grat1/grat2
    grat_dict = define_grating(dataset["grating_check"])
        
def define_grating_op21(dataset_curr):
    #checking the gratings for grat1/grat2/grat3
    qsi_data = spect_calib_cxrs_get_data(dataset_curr)

    #For each grating, the shift of the spectra is determined as the central wavelength is changed 
    gratings = np.unique(qsi_data.get_coordinate_object("Grating").values)
    meas_resolution = dict()
    for grating in gratings:
        grat_data = qsi_data.slice_data(slicing={"Grating":grating})
        grat_data = grat_data.slice_data(slicing={"ROI":np.unique(grat_data.coordinate("Relevant ROI")[0])})
        wavelengths = np.sort(np.unique(grat_data.coordinate("Central Wavelength")[0]))
        ref = grat_data.slice_data(slicing={"Central Wavelength": wavelengths[0]})
        shift = np.zeros(len(wavelengths))
        wave_index = 1
        for wavelength in wavelengths[1:]:
            curr = grat_data.slice_data(slicing={"Central Wavelength": wavelength})
            cross_corr = np.abs(fft.ifft(fft.fft(np.concatenate([curr.data,
                        np.zeros(len(curr.data))]))*np.conj(fft.fft(np.concatenate([np.zeros(len(ref.data)),
                        ref.data])))))
            shift[wave_index] = np.argmax(cross_corr)-len(ref.data)
            # from matplotlib import pyplot as plt

            # plt.plot(cross_corr, label=f"{curr.coordinate('Grating')[0][0]} Shift: {shift[wave_index]}")
            # plt.legend()

            # plt.plot(ref.coordinate("Coord 2")[0], ref.data, label="540nm")
            # plt.plot(ref.coordinate("Coord 2")[0], curr.data, label="541nm")
            # plt.ylabel("Intensity [n.a]")
            # plt.xlabel("Pixel No.")
            # plt.legend()
            wave_index += 1

        meas_resolution[grating] = np.abs(1000/np.polyfit(wavelengths,shift,1)[0])
    #FInding gratings with closest nominal resolution
    grat_dict = dict()
    nominal_resolution = {29.3: 1200,
                          17.7: 1800,
                          11.3: 2400 }
    res_keys = np.asarray(list(nominal_resolution.keys()))
    for grating in gratings:
        error = np.min(np.abs(res_keys-meas_resolution[grating]))
        grat_dict[grating] = [nominal_resolution[res_keys[np.argmin(np.abs(res_keys-meas_resolution[grating]))]], f"Resolution error {error}"]
    
    # Testing if the results make any sense
    outputs = np.asarray([grat_dict[grating][0] for grating in grat_dict.keys()])
    if len(outputs) != len(np.unique(outputs)):
        raise ValueError("Same grating obtained for different spectrometer settings, unable to define grating automatically")
    return grat_dict
    
def spect_calib_cxrs_get_data_op21(dataset_curr):
    cxrs_data = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              options={'Scale Time': True,
                                       'Cache Data': False,
                                       'Time Start': dataset_curr[0][0],
                                       'Time Stop': dataset_curr[0][1]},
                              object_name='CXRS_grating_check')
    
    roi_coord_flap = flap.Coordinate(name='ROI', unit=1, mode=flap.CoordinateMode(equidistant=False), shape=(6),
                                  values=["P01","P02","P03","P04","P05","P06"], dimension_list=[1])
    cxrs_data.add_coordinate_object(roi_coord_flap)
    cxrs_data.del_coordinate("Coord 1")

    rel_roi_coord = np.array([desc.split("-")[0] for desc in dataset_curr[1]])
    rel_roi_coord_flap = flap.Coordinate(name='Relevant ROI', unit=1, mode=flap.CoordinateMode(equidistant=False), shape=rel_roi_coord.shape,
                                  values=rel_roi_coord, dimension_list=[0])
    cxrs_data.add_coordinate_object(rel_roi_coord_flap)
    wavelength_coord = np.array([int(desc.split("-")[2]) for desc in dataset_curr[1]])
    wavelength_coord_flap = flap.Coordinate(name='Central Wavelength', unit="nm", mode=flap.CoordinateMode(equidistant=False), shape=wavelength_coord.shape,
                                  values=wavelength_coord, dimension_list=[0])
    cxrs_data.add_coordinate_object(wavelength_coord_flap)
    try:
        grat_coord = np.array([int(desc.split("-")[1]) for desc in dataset_curr[1]])
        grat_coord_flap = flap.Coordinate(name='Grating', unit="gr/mm", mode=flap.CoordinateMode(equidistant=False), shape=wavelength_coord.shape,
                                      values=wavelength_coord, dimension_list=[0])
    except ValueError:
        grat_coord = np.array([desc.split("-")[1] for desc in dataset_curr[1]])
        grat_coord_flap = flap.Coordinate(name='Grating', unit="1", mode=flap.CoordinateMode(equidistant=False), shape=grat_coord.shape,
                                      values=grat_coord, dimension_list=[0])
    cxrs_data.add_coordinate_object(grat_coord_flap)
    
   
    return cxrs_data

# def autocorr(dataobject):
#     for ROI in range(dataobject.get_coordinate_object("Channel").shape):
        
#     pass
    
