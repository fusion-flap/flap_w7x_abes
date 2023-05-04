# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:23:49 2018

@author: Miklos Vecsei, Fusion Plasma Physics Department, Centre for Energy Research

Spatial calibration of the for W7-X alkali BES diagnostic module for flap
"""

import os
import sys
import warnings
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
                                ["P04-grat1-540","P04-grat1-541","P04-grat2-540","P04-grat2-541","P04-grat3-540","P04-grat3-541"]]
    dataset["540_2400"] = [[datetime.datetime(2023, 4, 26, 13, 57, 00),
                            datetime.datetime(2023, 4, 26, 14, 13, 00)],
                           ["P01","P02","P03","P04","P05","P06"]]
    dataset["540_1800"] = [[datetime.datetime(2023, 4, 26, 14, 30, 00),
                            datetime.datetime(2023, 4, 26, 14, 54, 00)],
                           ["P01","P02","P03","P04","P05","P06"]]
    
    #checking the gratings for grat1/grat2/grat3
    grating_check = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                              options={'Scale Time': True,
                                       'Cache Data': False,
                                       'Time Start': dataset["grating_check"][0][0],
                                       'Time Stop': dataset["grating_check"][0][1]},
                              object_name='CXRS_grating_check')