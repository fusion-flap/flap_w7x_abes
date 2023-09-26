#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 15:16:49 2023

@author: mvecsei
"""

import os
import warnings
import time
import matplotlib.pyplot as plt
import numpy as np

import flap
from .w7x_abes import abes_get_config
from . import spatcal

#functions to add to the cmos dataobject
class CMOS_data(flap.DataObject):
    def get_chopstate(self,chop=0, defl=0):
        o={'State':{'Chop': chop, 'Defl': defl}}
        d_beam_st=flap.get_data('W7X_ABES',
                                exp_id=self.exp_id,
                                name='Chopper_time',
                                coordinates={'Time':[np.min(self.coordinate('Time')[0]), np.max(self.coordinate('Time')[0])]},
                                options=o,\
                                object_name=f'CMOS chop {chop} defl {defl}',
                                )
        times = d_beam_st.coordinate('Time')[1]/2+d_beam_st.coordinate('Time')[2]/2
        d = self.slice_data(slicing={'Time':times})
        return d

def w7x_abes_cmos_get_data(exp_id=None, data_name="CMOS", no_data=False, options=None, coordinates=None, data_source=None):
    """ Data read function for CMOS camera of the W7-X Alkali BES diagnostic
    data_name: Not required
      
    exp_id: Experiment ID, YYYYMMDD.xxx
    Unix style regular expressions are allowed with * and []
                       Can also be a list of data names, eg. ['ABES-1','ABES-4']
    coordinates: Not yet implemented: List of flap.Coordinate() or a single flap.Coordinate
                 Defines read ranges. The following coordinates are interpreted:
                     'Time': The read times
                     Only a single equidistant range is interpreted in c_range.
    options:
        'Datapath': Data path (string)
        For further options see Chopper_times see shopper_timing_data()
      
    """
    if (exp_id is None):
        raise ValueError('exp_id should be set for W7X ABES.')
    default_options = {'Datapath': 'data',
                       'Spatial calibration': False,
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
    if os.path.exists(xmlfile) is False:
        raise ValueError('XML file '+xmlfile+' does not exist')
    try:
        while os.access(xmlfile, os.R_OK) is False:
            time.sleep(0.1)
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
    
    images = get_images(datapath)
    image_array = read_images(datapath, images)
    time_coord = get_time(xml, image_array.shape[0])
    
    unit = flap.Unit(name='Intensity', unit='a.u.')
    d = CMOS_data(data_array=image_array, data_unit=unit,
                  coordinates=[time_coord], exp_id=exp_id,
                  data_title='CMOS', data_shape=image_array.shape)
    
    if _options['Spatial calibration'] is True:
        # Getting the spatial calibration
        d = w7x_abes_cmos_add_coordinate(d, ['Device x', 'Device y', 'Device R', 'Beam axis'],
                           options={"Shot spatcal dir": flap.config.get("Module W7X_ABES","Spatial calibration directory")})
    
    return d

def get_images(datapath):
    'reads all bmp or png images in the folder and orders them according to the id of the taken image'
    files = os.listdir(datapath)
    bmp_images = [filename for filename in files if ".bmp" in filename]
    png_images = [filename for filename in files if ".png" in filename]
    
    if len(bmp_images) > 0 and len(png_images)>0:
       warnings.warn(f"Both png and bmp images are present in {datapath}. Using only the bmp images")
       images = bmp_images
    elif len(bmp_images)>0:
        images = bmp_images
    elif len(png_images)>0:
        images = png_images
    else:
        raise OSError(f"No bmp or png files in {datapath}")
    
    image_id = [int(filename.split("_")[0]) for filename in images]
    
    images = [image for _,image in sorted(zip(image_id, images))]
    
    return images

def read_images(datapath, images):
    data_array = np.asarray([np.sum(np.asarray(plt.imread(os.path.join(datapath, filename))), axis=2)
                  for filename in images], dtype="int")
    return data_array


def get_time(xml, length):
    frametime = float(xml.get_element(section="CMOS", element="Frametime")["Value"])/1000
    time_coord = flap.Coordinate(name='Time', unit='Second', start=frametime/2, step=[frametime],
                                 shape=[length],
                                 mode=flap.CoordinateMode(equidistant=True),
                                 dimension_list=[0])

    return time_coord

def register(data_source=None):
    flap.register_data_source('W7X_ABES_CMOS', get_data_func=w7x_abes_cmos_get_data, add_coord_func=w7x_abes_cmos_add_coordinate)


def w7x_abes_cmos_add_coordinate(data_object,
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
        exp_spatcal = spatcal.ShotSpatCalCMOS(data_object.exp_id, options=_options)
    else:
        exp_spatcal = spatcal.ShotSpatCalCMOS(exp_id, options=_options)
    try:
        exp_spatcal.read(options=_options)
    except OSError:
        # if the file does not exist
        exp_spatcal.generate_shotdata(options={'Plot': False, 'Overwrite': False, "Shot spatcal dir": _options["Shot spatcal dir"]})
        exp_spatcal.read(options=_options)

    #Adding the signal names
    if exp_id is None:
        exp_id = data_object.exp_id
    
    #Spatial coordinates
    for coord_name in coordinates:
        coord_object = exp_spatcal.create_coordinate_object([1, 2], coord_name)
        data_object.add_coordinate_object(coord_object)


    return data_object
