# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:23:49 2018

@author: Vecsei

Spatial calibration of the for W7-X alkali BES diagnostic module for flap
"""

import os.path
import numpy as np
import h5py
import flap

class ShotSpatCal(flap.DataObject):
    def __init__(self, shotID, options=None):
        # read the spatcal hdf5 file
        if options is not None:
            options_list = options.keys()
        else:
            options_list = dict()
        if 'spatcal_dir' in options_list:
            spatcal_dir = options['spatcal_dir']
        else:
            spatcal_dir ='./'
        filename = shotID + '_spat.cal'
        fn = os.path.join(spatcal_dir,filename)
        with h5py.File(fn, "r", libver='earliest') as h5File: 
            channel_names = np.array(h5File['/Channels'].value)
            fibres = h5File["/Fibres"].value
            self.calibration = h5File["/Calibration"].value[0].decode("utf-8") 
            self.beam_start_im = h5File["/Beam_start_im"].value
            self.beam_end_im = h5File["/Beam_end_im"].value
            fibre_coords_im = h5File["/Fibre_coords_im"].value
            self.beam_start_beam = h5File["/Beam_start_beam"].value
            self.beam_end_beam = h5File["/Beam_end_beam"].value
            fibre_coords_beam = h5File["/Fibre_coords_beam"].value
            self.beam_start_xyz = h5File["/Beam_start_xyz"].value
            self.beam_end_xyz = h5File["/Beam_end_xyz"].value
            fibre_coords_xyz = h5File["/Fibre_coords_xyz"].value
            self.beam_plane_vector = h5File["/Beam_plane_vector"].value
        
        if 'Channels' in options.keys():
            channels = [bytes(channel, "utf-8") for channel in options["Channels"]]
        else:
            channels=[bytes('ABES-'+str(index), 'utf-8') for index in range(1,41)]

        relative_loc = np.zeros(len(channel_names))
        index=0
        for ch in channels:
            relative_loc[np.where(channel_names==ch)[0][0]] = index+1
            index=index+1
        relative_loc = [x-1 if x > 0 else len(channel_names) for x in relative_loc]
        relative_loc = np.argsort(relative_loc)

        channel_names = channel_names[relative_loc]
        self.fibres = fibres[relative_loc]
        self.fibre_coords_im = fibre_coords_im[:, relative_loc]
        fibre_coords_beam = fibre_coords_beam[:,relative_loc]
        fibre_coords_xyz = fibre_coords_xyz[:,relative_loc]
        
        self.data={
                "Channel name": np.array(channel_names),
                "Device x": fibre_coords_xyz[0,:],
                "Device y": fibre_coords_xyz[1,:],
                "Device z": fibre_coords_xyz[2,:],
                "Device R": np.sqrt(fibre_coords_xyz[0,:]**2+fibre_coords_xyz[1,:]**2+fibre_coords_xyz[2,:]**2),
                "Device Z": fibre_coords_xyz[2,:],
                "Beam axis": fibre_coords_beam[0,:]
                }
        return

    def create_coordinate_object(self, dimension_list, coord_name, channel_names=None):
        if not (channel_names is None):
            data = np.zeros(len(channel_names))
            index = 0
            for channel in channel_names:
                data[index] = (self.data[coord_name])[np.where(self.data["Channel name"] == bytes(channel, 'utf-8'))[0][0]]
                index = index+1
        else:
            data = self.data[coord_name]

        coord_object = flap.Coordinate(name=coord_name, unit='m', values=data, shape=np.shape(data),
                          mode=flap.CoordinateMode(equidistant=False), dimension_list=dimension_list)

        return coord_object

if __name__ == '__main__':
    a = ShotSpatCal('20180912.040', options={"spatcal_dir": "./tests/"})
    channel_names = ['ABES-1', 'ABES-2']
    b=a.create_coordinate_object([0], 'Beam axis', channel_names=channel_names)