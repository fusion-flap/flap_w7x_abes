#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:10:03 2020

@author: mvecsei
"""

import flap
import flap_w7x_abes
import os
import numpy as np

if __name__ == '__main__':
#    a = ShotSpatCal('20180912.040', options={"spatcal_dir": "./tests/"})
#    channel_names = ['ABES-1', 'ABES-2']
#    b=a.create_coordinate_object([0], 'Beam axis', channel_names=channel_names)
    a = flap_w7x_abes.ShotSpatCal('20181016.008')
#    a = flap_w7x_abes.ShotSpatCal('20171207.024')

    
#    a.lab_calib(options)
    options = {'Get CMOS to machine': False,'Get APDCAM to CMOS': True, 'Spherical symmetry': True}
#    a.full_calib(options=options)
#    a.generate_shotdata()
    import h5py
    old = h5py.File('/media/mvecsei/DATA/data/W7-X/APDCAM/spatcal/20181016.008_spat.cal', 'r')
    new = flap.load('/DATA/repos/flap/modules/flap_w7x_abes/spatcal/20181016.008_spat.cal')

#    h5File = h5py.File('/media/mvecsei/DATA/data/W7-X/APDCAM/spatcal/20171207.024_spat.cal', 'r')
    from matplotlib import pyplot as plt
#    data = h5File['Fibre_coords_xyz'].value
#    names = h5File['Channels'].value
#    index=0
#    data_rel_x=[]
#    data_rel_y=[]
#    for chan in names:
#        if (chan[:4] == b'ABES'):
#            data_rel_x+=[data[0][index]]
#            data_rel_y+=[data[1][index]]
#        index=index+1
#    data_rel_x= np.asarray(data_rel_x)
#    data_rel_y= np.asarray(data_rel_y)
#
#    points = np.sort(np.sqrt(data_rel_x**2+data_rel_y**2))
#    plt.plot(points, label='old')
#    plt.scatter(data[0], data[1])
#    plt.legend()

    plt.plot([1312-old['Beam_start_im'].value[0],1312-old['Beam_end_im'].value[0]],1082-np.asarray([old['Beam_start_im'].value[1],old['Beam_end_im'].value[1]]), color='white')
    plt.scatter(1312-old['Beam_start_im'].value[0],1082-old['Beam_start_im'].value[1], color='white')
    plt.plot([new['Beam_start_im'][0], new['Beam_end_im'][0]],np.asarray([new['Beam_start_im'][1], new['Beam_end_im'][1]]), color='red')
    plt.scatter(new['Beam_start_im'][0], new['Beam_start_im'][1], color='red')
    index=0
    im=(np.asarray(plt.imread('./20181016.008.png')))
    plt.imshow(im)
    for chan in h5File['Channels'].value:
        plt.scatter(1312-old['Fibre_coords_im'].value[0],1082-old['Fibre_coords_im'].value[1], color='white', marker = 'x')
        plt.scatter(new['Fibre_coords_im'][0], new['Fibre_coords_im'][1], color='red', marker = 'o')
    plt.show()

    
#    plt.scatter([h5File['Beam_start_xyz'].value[0]], [h5File['Beam_start_xyz'].value[1]], color='blue', marker='x')
#    plt.plot([h5File['Beam_start_xyz'].value[0],h5File['Beam_end_xyz'].value[0]], [h5File['Beam_start_xyz'].value[1],h5File['Beam_end_xyz'].value[1]], color='blue')
#    plt.legend()
#    index=0
#    for chan in h5File['Channels'].value:
#        if chan == b'ABES-1':
#            plt.scatter(h5File['Fibre_coords_xyz'].value[0],h5File['Fibre_coords_xyz'].value[1], color='blue', marker = 'x')
#        index+=1