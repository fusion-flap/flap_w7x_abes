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
import copy

from . import spatcal
from . import cxrs_util
from . import cxrs

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
        "Cache data": Whether to cache the data to webapi
        "Spatial calibration": Whether to obtain the spatial location of the channels
        For further options see Chopper_times see shopper_timing_data()

    """
    if (exp_id is None):
        raise ValueError('exp_id should be set for W7X ABES.')
    default_options = {'Cache data': True,
                       "Spatial calibration": True
                       }
    _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES_CXRS')
    d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              exp_id=exp_id,
                              options={'Scale Time': True,
                                       'Cache Data': _options["Cache data"]},
                              object_name=data_name)
    time_coord = copy.deepcopy(d.get_coordinate_object("Time"))
    time_coord.start = time_coord.values[0]
    time_coord.step = np.median(time_coord.values[1:]-time_coord.values[:-1])
    time_coord.mode.equidistant = True
    d.del_coordinate("Time")
    d.add_coordinate_object(time_coord)
    d.data_source = data_source
    channel_coord = copy.deepcopy(d.get_coordinate_object("Coord 1"))
    channel_coord.unit.name = "Channel"
    d.del_coordinate("Coord 1")
    d.add_coordinate_object(channel_coord)
    
    
    data_source = [line for line in d.info.split("\n") if "Source:" in line][0].split("Source:")[1][1:]
    
    # getting the grating, central wavelength setting, slit width
    central_wavelength = np.mean(webapi.get_data(exp_id,
                                                 data_name = f"{data_source.split('DATASTREAM')[0]}PARLOG/parms/IsoPlaneControl/config/centralWavelength/").data)
    grating = int(np.mean(webapi.get_data(exp_id,
                                      data_name = f"{data_source.split('DATASTREAM')[0]}PARLOG/parms/IsoPlaneControl/config/gratingNumber/").data))
    grating_dict={2: "1200g_per_mm", 
                  1: "1800g_per_mm",
                  0: "2400g_per_mm"}
    slit_width = np.mean(webapi.get_data(exp_id,
                                         data_name = f"{data_source.split('DATASTREAM')[0]}PARLOG/parms/IsoPlaneControl/config/slitWidth/").data)
    
    d.config = {"Central wavelength": central_wavelength,
                "Grating": grating_dict[grating],
                "Slit width": slit_width}

    
    if _options['Spatial calibration'] is True:
         # Getting the spatial calibration
         try:
             spatcal_dir = flap.config.get("Module W7X_ABES","Spatial calibration directory")
             d = cxrs_add_coordinate(d, ['Device R', 'Device x', "Device y"],
                                options={"Shot spatcal dir": flap.config.get("Module W7X_ABES","Spatial calibration directory"),
                                         "Cache data": _options["Cache data"]})
         except:
             d = cxrs_add_coordinate(d, ['Device R', 'Device x', "Device y"],
                                options={"Cache data": _options["Cache data"]})
         
    return d

def cxrs_add_coordinate(data_object,
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
                       "Cache data": True
                       }
    _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES')

    if exp_id is None:
        exp_spatcal = spatcal.ShotSpatCalCXRS(data_object.exp_id, options=_options)
    else:
        exp_spatcal = spatcal.ShotSpatCalCXRS(exp_id, options=_options)
    try:
        exp_spatcal.read(options=_options)
    except (OSError, IOError):
        # if the file does not exist
        exp_spatcal.generate_shotdata(options={'Plot': False, 'Overwrite': False, "Shot spatcal dir": _options["Shot spatcal dir"]})
        exp_spatcal.read(options=_options)

    # getting the dimension of the channel coordinate, this should be the same as the spatial coordinate
    data_coord_list = np.array([coord.unit.name for coord in data_object.coordinates])

    #Adding the signal names
    if exp_id is None:
        exp_id = data_object.exp_id

    #Spatial coordinates
    data_source = [line for line in data_object.info.split("\n") if "Source:" in line][0].split("Source:")[1][1:]
    d = flap.get_data('W7X_WEBAPI', name=f"{data_source.split('DATASTREAM')[0]}PARLOG/parms/PICam/config/Rois_configString/",
                      exp_id=exp_id,
                      options={'Scale Time': True,
                               'Cache Data': _options["Cache data"]},
                      object_name="ROIs")
    patch_config, temp1, temp2 = cxrs_util.read_fibre_config(exp_id)
    roi_names = [ROI.split(",")[0] for ROI in d.data[0].split("{")[2:] if "Enabled" in ROI]

    if int(exp_id[:4])<2024:
        roi_names = ["A."+str(roi[-1]).zfill(2) for roi in roi_names]
    if "channel" in roi_names[0]: #At some point Barnabas changed the naming convention of the ROIs
        roi_names = [ROI.split(",")[0] for ROI in d.data[0].split("{")[1:] if "Enabled" in ROI]
        roi_names = ["A."+str(roi.split("channel ")[1]).zfill(2) for roi in roi_names]
    signal_names = copy.deepcopy(data_object.get_coordinate_object('Channel'))
    signal_names.unit.name = "Fiber name"
    signal_names.values = np.asarray(["N.A." for index in range(len(roi_names))])
    signal_names.mode.equidistant = False
    optical_channel = copy.deepcopy(signal_names.values)
    for key in patch_config.keys():
        try:
            signal_names.values[np.where(patch_config[key] == np.asarray(roi_names))[0]] = patch_config[key]
            optical_channel[np.where(patch_config[key] == np.asarray(roi_names))[0]] = key
        except IndexError:
            # the referred channel was not measured for
            pass
    optical_channel_coord = copy.deepcopy(signal_names)
    optical_channel_coord.unit.name = "Optical channel"
    optical_channel_coord.values = optical_channel
    
    data_object.add_coordinate_object(signal_names)
    data_object.add_coordinate_object(optical_channel_coord)

    channel_coordinate_dim = data_object.get_coordinate_object("Fiber name").dimension_list[0]
    channel_names = data_object.get_coordinate_object("Fiber name").values

    for coord_name in coordinates:
        coord_object = exp_spatcal.create_coordinate_object(channel_coordinate_dim, coord_name,
                                                         channel_names=channel_names)
        data_object.add_coordinate_object(coord_object)

    #Wavelength
    data_object = cxrs_add_wavelength(data_object)

    return data_object

def cxrs_add_wavelength(data_object):
       
    grid = data_object.config["Grating"]
    w_s = data_object.config["Central wavelength"]
    try:
        datapath_base = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cxrs_util_data')
    except NameError:
        datapath_base = "cxrs_util_data"
    if(int(data_object.exp_id[:8]) < 20231231):
        wl_values = np.zeros((6, 1024))
        for i in range(1, 7):
            wl_values[i-1, :] = cxrs.wavelength_grid_generator_op21(
                grid, w_s, i, datapath_base)

        wl_coord_flap = flap.Coordinate(name='Wavelength', unit="nm",
                        mode=flap.CoordinateMode(equidistant=False),
                        shape=(6,1024),values=wl_values,dimension_list=[1,2])
        data_object.add_coordinate_object(wl_coord_flap)
        data_object.del_coordinate("Coord 2")
        
    else:
        wl_values = np.zeros((data_object.coordinate("Device R")[0].shape[1], 1024))
        for i in range(1, data_object.coordinate("Device R")[0].shape[1]+1):
            wl_values[i-1, :] = cxrs.wavelength_grid_generator_op22(
                grid, w_s, datapath_base)

        wl_coord_flap = flap.Coordinate(name='Wavelength', unit="nm",
                        mode=flap.CoordinateMode(equidistant=False),
                        shape=(data_object.coordinate("Device R")[0].shape[1],1024),
                        values=wl_values,dimension_list=[1,2])
        data_object.add_coordinate_object(wl_coord_flap)
        data_object.del_coordinate("Coord 2")
    
    return data_object