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
from scipy.ndimage import median_filter


if __name__ == '__main__':

    shotID = '20230307.047'
    a = flap_w7x_abes.ShotSpatCal(shotID)
    images = a.calc_chan_range()
    adf
    # full_calib = True
    # if full_calib is True:
    #     a = flap_w7x_abes.ShotSpatCal(shotID)
    #     options = {'Get CMOS to machine': True,'Get APDCAM to CMOS': False, 'Get CXRS to CMOS': False, 'Circular symmetry':True, 'Elliptical symmetry':False,'Noise limit': 200}
    #     a.full_calib(options=options)


    # filename = os.path.join('/media/mvecsei/DATA/repos/flap/modules/flap_w7x_abes/spatcal/2021/Geometry', 'apdcam_to_cmos.hdf5')
    # apdcam_to_cmos = flap.load(filename)
    # print(apdcam_to_cmos['APDCAM channel centers'])
    # filename = os.path.join('/media/mvecsei/DATA/repos/flap/modules/flap_w7x_abes/spatcal/2021/Geometry', 'cxrs_to_cmos.hdf5')
    # apdcam_to_cmos = flap.load(filename)
    # print(apdcam_to_cmos['CXRS channel centers'])
    # asdf

    shot_calib = True
    if shot_calib is True:
        # a = flap_w7x_abes.ShotSpatCalCXRS(shotID)
        a = flap_w7x_abes.ShotSpatCalCMOS(shotID)
        a.generate_shotdata(options={'Plot': True, 'Overwrite': True})
        a.read()
        a = flap_w7x_abes.ShotSpatCal(shotID)

    raise ValueError('stop')


#    a.lab_calib(options)
    options = {'Get CMOS to machine': True ,'Get APDCAM to CMOS': False, 'Circular symmetry':True,
               'Flip horizontally': True, 'Noise limit': 200}
#    a.full_calib(options=options)
#    a.generate_shotdata(options={'Overwrite': True})
    import h5py
    old = h5py.File('/media/mvecsei/DATA/data/W7-X/APDCAM/spatcal/20181016.008_spat.cal', 'r')
    old_unflipped = h5py.File('/media/mvecsei/DATA/repos/flap/modules/flap_w7x_abes/tests/20181016.008_spat.cal', 'r')
    new = flap.load('/media/mvecsei/DATA/repos/flap/modules/flap_w7x_abes/spatcal/20240222.012_cxrs_spat.cal')


    whatnot = dict()
    index = 0
    for channel in new['Channels']:
        whatnot[channel] = [new['Fibre_coords_xyz'][0][index], new['Fibre_coords_xyz'][1][index]]
        index += 1
    keys = sorted(whatnot.keys())
    for key in keys:
        print(f"{key} {whatnot[key]}")
        plt.scatter(whatnot[key][0], whatnot[key][1])
        plt.text(whatnot[key][0], whatnot[key][1], key)
#    h5File = h5py.File('/media/mvecsei/DATA/data/W7-X/APDCAM/spatcal/20171207.024_spat.cal', 'r')
    from matplotlib import pyplot as plt
#
#    plt.style.use('dark_background')    
#    plt.plot([1312-old['Beam_start_im'].value[0],1312-old['Beam_end_im'].value[0]],1082-np.asarray([old['Beam_start_im'].value[1],old['Beam_end_im'].value[1]]), color='white')
#    plt.scatter(1312-old['Beam_start_im'].value[0],1082-old['Beam_start_im'].value[1], color='white')
#    plt.plot([old_unflipped['Beam_start_im'].value[0],old_unflipped['Beam_end_im'].value[0]],1082-np.asarray([old_unflipped['Beam_start_im'].value[1],old_unflipped['Beam_end_im'].value[1]]), color='magenta')
#    plt.scatter(old_unflipped['Beam_start_im'].value[0],1082-old_unflipped['Beam_start_im'].value[1], color='magenta')
#    plt.plot([new['Beam_start_im'][0], new['Beam_end_im'][0]],np.asarray([new['Beam_start_im'][1], new['Beam_end_im'][1]]), color='red')
#    plt.scatter(new['Beam_start_im'][0], new['Beam_start_im'][1], color='red')
#    im=(np.asarray(plt.imread('./20181016.008.png')))
#    plt.imshow(im)
#    plt.scatter(old_unflipped['Fibre_coords_im'].value[0],1082-old_unflipped['Fibre_coords_im'].value[1], color='magenta', marker = 'x')
#    plt.scatter(1312-old['Fibre_coords_im'].value[0],1082-old['Fibre_coords_im'].value[1], color='white', marker = 'x')
#    plt.scatter(new['Fibre_coords_im'][0], new['Fibre_coords_im'][1], color='red', marker = 'o', alpha = 0.6)
#    plt.show()
#
#    plt.style.use('dark_background')    
#    plt.scatter(new['Beam_start_xyz'][0], [new['Beam_start_xyz'][1]], color='red', marker='o')
#    plt.plot([new['Beam_start_xyz'][0], new['Beam_end_xyz'][0]], [new['Beam_start_xyz'][1], new['Beam_end_xyz'][1]], color='red', marker='o')
#    plt.scatter([old['Beam_start_xyz'][0], old['Beam_end_xyz'][0]], [old['Beam_start_xyz'][1], old['Beam_end_xyz'][1]], color='white', marker='x')
#    plt.scatter(old['Fibre_coords_xyz'].value[0], old['Fibre_coords_xyz'].value[1], color='white', marker = 'x')
#    plt.scatter(new['Fibre_coords_xyz'][0],new['Fibre_coords_xyz'][1], color='red', marker = 'o')
#    plt.scatter([old_unflipped['Beam_start_xyz'][0], old['Beam_end_xyz'][0]], [old_unflipped['Beam_start_xyz'][1], old_unflipped['Beam_end_xyz'][1]], color='magenta', marker='x')
#    plt.scatter(old_unflipped['Fibre_coords_xyz'].value[0], old_unflipped['Fibre_coords_xyz'].value[1], color='magenta', marker = 'x')
#    plt.show()
#    
#    plt.style.use('dark_background')    
#    index=0
#    chp=0
#    for chan in new['Channels']:
#        if chan[:4] == 'ABES':
#            plt.scatter(int(chan[5:]),np.sqrt(new['Fibre_coords_xyz'][0][index]**2+new['Fibre_coords_xyz'][1][index]**2), color='red', marker = 'o')
#            chp +=1
#        index+=1
#    index=0
#    chp=0
#    for chan in old['Channels'].value:
#        if chan[:4].decode('utf-8') == 'ABES':
#            plt.scatter(int(chan[5:].decode('utf-8')),np.sqrt(old['Fibre_coords_xyz'].value[0][index]**2+ old['Fibre_coords_xyz'].value[1][index]**2), color='white', marker = 'x')
#            chp +=1
#        index+=1
#    index=0
#    chp=0
#    for chan in old_unflipped['Channels'].value:
#        if chan[:4].decode('utf-8') == 'ABES':
#            plt.scatter(int(chan[5:].decode('utf-8')),np.sqrt(old_unflipped['Fibre_coords_xyz'].value[0][index]**2+ old_unflipped['Fibre_coords_xyz'].value[1][index]**2), color='magenta', marker = 'x')
#            chp +=1
#        index+=1
#    plt.show()
#
#    with open('../spatcal/2018/fibre_calib_images.dat') as file:
#         from scipy.ndimage import median_filter
#         lines = file.readlines()
#         for line_index in range(1, 61-3):
#             line_data = [data for data in lines[line_index].split("\t") if data is not ""]
#             filename = '../spatcal/2018/'+line_data[-1].rstrip()[:-3]+'png'
#             image=(np.asarray(plt.imread(filename)))
#             image = image/np.var(image)
#             image = image-np.median(image)
#             image[image<100]=0
#             image = median_filter(image, size=3)
#             if 'all_data' not in locals():
#                 all_data = image
#             else:
#                 all_data=all_data+image
#    plt.imshow(all_data)
#    plt.scatter(new['Fibre_coords_im'][0], new['Fibre_coords_im'][1], color='red', marker = 'o', alpha = 0.6)
#    plt.show()

    # plotting the noise to check whether the images were flipped
#    im=(np.asarray(plt.imread('./20181016.008.png')))
#    im = im/np.var(im)
#    image = median_filter(im, size=5)
#    im = np.abs(im-image)
#    limit=180
#    im[im<limit]=0
#    im[im>1.5*limit]=1
#    im[im>limit]=0.5
#    im[1:-1,1:-1] = im[1:-1,1:-1]+\
#                    im[2:,1:-1]+im[:-2,1:-1]+im[1:-1,2:]+im[1:-1,:-2]+\
#                    im[2:,2:]+im[2:,:-2]+im[:-2,2:]+im[:-2,:-2]
#    im[im>limit]=1
#    dx=200
#    dy=200
#    color=0.5
#    im[:,:-1:dy] = color
#    im[:,1::dy] = color
#    im[:-1:dx,:] = color
#    im[1::dx,:] = color
#    plt.subplot(1,2,1)
#    plt.title('stellarator image')
#    plt.imshow(im, cmap='gray', interpolation='none', vmin=0, vmax=np.max(im))
#    im=(np.asarray(plt.imread('../spatcal/2018/00000013_00000002B43D3F36_238941562.png')))
#    im = im/np.var(im)
#    image = median_filter(im, size=5)
#    im = np.abs(im-image)
#    limit=900
#    im[im<limit]=0
#    im[im>1.5*limit]=1
#    im[im>limit]=0.5
#    im[1:-1,1:-1] = im[1:-1,1:-1]+\
#                    im[2:,1:-1]+im[:-2,1:-1]+im[1:-1,2:]+im[1:-1,:-2]+\
#                    im[2:,2:]+im[2:,:-2]+im[:-2,2:]+im[:-2,:-2]
#    im[im>limit]=1
#    dx=200
#    dy=200
#    color=0.5
#    im = np.fliplr(im)
#    im[:,:-1:dy] = color
#    im[:,1::dy] = color
#    im[:-1:dx,:] = color
#    im[1::dx,:] = color
#    im = np.fliplr(im)
#    plt.subplot(1,2,2)
#    plt.title('laboratory image')
#    plt.imshow(im, cmap='gray', interpolation='none', vmin=0, vmax=np.max(im))
#    plt.show()
    
    # Plotting the IDL fudicial points on the IDL calibration image using the flap spatial calibration
    fiducial_curve = []
    fiducial_curve += [[[6.28495, -0.844860, 73.6651],
                       [6.27878, -0.849542, 73.4806],
                       [6.26751, -0.853739, 73.3160],
                       [6.25187, -0.857165, 73.1824],
                       [6.23290, -0.859586, 73.0889],
                       [6.21187, -0.860837, 73.0421]]]
    fiducial_curve += [[[6.21684, -0.824006, 74.5203],
                       [6.23751, -0.825265, 74.4618],
                       [6.25582, -0.827693, 74.3578],
                       [6.27053, -0.831125, 74.2157],
                       [6.28068, -0.835326, 74.0452],
                       [6.28559, -0.840011, 73.8576]]]
    fiducial_curve += [[[6.19519, -0.824001, 74.5291],
                       [6.17404, -0.825252, 74.4870],
                       [6.15485, -0.827673, 74.3963],
                       [6.13893, -0.831099, 74.2630],
                       [6.12741, -0.835296, 74.0962],
                       [6.12110, -0.839978, 73.9073]]]
    fiducial_curve += [[[6.12044, -0.844827, 73.7097],
                       [6.12547, -0.849512, 73.5174],
                       [6.13585, -0.853713, 73.3440],
                       [6.15083, -0.857145, 73.2016],
                       [6.16938, -0.859573, 73.1000],
                       [6.19021, -0.860832, 73.0457]]]
    fiducial_curve += [[[6.1117644227326, -0.7874695, 74.2139156147943],
                       [6.1053344084946, -0.7824952, 74.2315295925234],
                       [6.0989049724530, -0.7775208, 74.2491807090633],
                       [6.0924761164383, -0.7725465, 74.266869078637],
                       [6.0860478422887, -0.7675722, 74.2845948159218],
                       [6.0796201518498, -0.7625978, 74.3023580360513],
                       [6.0731930469751, -0.7576235, 74.3201588546178]]]
##    
    fiducial_curve = []
    fiducial_curve += [[[6.27608, -0.845267, 73.6868],
                       [6.27050, -0.849499, 73.5198],
                       [6.26030, -0.853291, 73.3709],
                       [6.24615, -0.856387, 73.2502],
                       [6.22900, -0.858574, 73.1658],
                       [6.20999, -0.859704, 73.1238]]]
    fiducial_curve += [[[6.21451, -0.826423, 74.4600],
                       [6.23319, -0.827561, 74.4072],
                       [6.24974, -0.829755, 74.3133],
                       [6.26305, -0.832856, 74.1849],
                       [6.27223, -0.836653, 74.0306],
                       [6.27667, -0.840886, 73.8609]]]
    fiducial_curve += [[[6.19495, -0.826420, 74.4678],
                       [6.17584, -0.827550, 74.4295],
                       [6.15850, -0.829737, 74.3475],
                       [6.14413, -0.832833, 74.2270],
                       [6.13373, -0.836625, 74.0763],
                       [6.12803, -0.840856, 73.9058]]]
    fiducial_curve += [[[6.12743, -0.845238, 73.7275],
                       [6.13197, -0.849471, 73.5539],
                       [6.14133, -0.853268, 73.3973],
                       [6.15486, -0.856369, 73.2686],
                       [6.17161, -0.858563, 73.1766],
                       [6.19042, -0.859700, 73.1273]]]
    fiducial_curve += [[[6.1117644227326, -0.7874695, 74.2139156147943],
                       [6.1053344084946, -0.7824952, 74.2315295925234],
                       [6.0989049724530, -0.7775208, 74.2491807090633],
                       [6.0924761164383, -0.7725465, 74.266869078637],
                       [6.0860478422887, -0.7675722, 74.2845948159218],
                       [6.0796201518498, -0.7625978, 74.3023580360513],
                       [6.0731930469751, -0.7576235, 74.3201588546178]]]
#    
    options = {
    'Spatial calib source dir': '/DATA/repos/flap/modules/flap_w7x_abes/spatcal/'}
    filename = os.path.join(options['Spatial calib source dir'], '2018/'
                                'Geometry', 'cmos_to_real.hdf5')
    projections = flap.load(filename)
    real_to_cmos = projections['Machine to CMOS projection']
#    real_to_cmos[0] = real_to_cmos[0]+30
#    real_to_cmos[4] = real_to_cmos[4]+45
    calibconf = flap_w7x_abes.spatcal.MachineCalibConfig(calibration_id='2018', options=options)
    calibconf.get_optical_axis_midplane_crosspoint()
    im=(np.asarray(plt.imread('./idl_calib_image_2017_2.png')))
#    im=(np.asarray(plt.imread('./20181016.008.png')))
    plt.imshow(im, cmap='gray')
    for curve in fiducial_curve:
        points_rzt = np.asarray(curve)
        point_XY = np.zeros([points_rzt.shape[0], 2])
        for index in range(points_rzt.shape[0]):
            point_xyz = np.asarray([points_rzt[index, 0]*np.cos(points_rzt[index, 2]*np.pi/180),
                                    points_rzt[index, 0]*np.sin(points_rzt[index, 2]*np.pi/180),
                                    points_rzt[index, 1]])
            # Calculating the XY coordinates - the coordinates in the plane perpendicular to the optical axis
            point_XY[index, :]=np.asarray(calibconf.get_image_XY_coord(point_xyz))
        cmos_check = flap_w7x_abes.spatcal.get_points_projection(point_XY, real_to_cmos)
        plt.scatter(cmos_check[:,0], cmos_check[:,1], marker='x')
    plt.scatter(real_to_cmos[0], real_to_cmos[4])
    plt.show()