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
try:
    import flap_w7x_webapi as webapi
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
except ModuleNotFoundError:
    warnings.warn("flap_w7x_webapi not found, unable to download cxrs data")
    
def w7x_abes_cxrs_get_data(exp_id=None, data_name=None, no_data=False, options=None, coordinates=None, data_source=None):
    """ Data read function for the W7-X Alkali BES CXRS diagnostic
    data_name: arbitrary string

    exp_id: Experiment ID, YYYYMMDD.xxx
            or a list of datetime.datetime type of variables, e.g. [datetime.datetime(2023, 2, 21, 13, 42, 00),
                                                                    datetime.datetime(2023, 2, 21, 14, 00, 00)]
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
    pass

def spect_calib_cxrs():
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
    print(grat_dict)
        
def define_grating(dataset_curr):
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
    
def spect_calib_cxrs_get_data(dataset_curr):
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

    
spect_calib_cxrs()