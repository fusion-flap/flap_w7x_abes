#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 14:28:09 2023

@author: apdcam
"""

from  flap_w7x_abes import cmos_main
import flap_w7x_abes
import flap
from matplotlib import pyplot as plt
import copy

flap_w7x_abes.register()

if __name__ == "__main__":
    shotID = "20230316.089"
    flap.config.read()

    # testing removing the background
    # cmos = cmos_main.w7x_abes_cmos_get_data(exp_id=shotID)
    # cmos_on = cmos.get_chopstate()
    # cmos_off = cmos.get_chopstate(chop=1)
    # for timeindex in range(cmos_on.data.shape[0]):
    #     plt.clf()
    #     plt.imshow(cmos_on.data[timeindex,:,:]-cmos_off.data[timeindex,:,:])
    #     plt.show()
    #     plt.pause(0.1)

    # testing adding the spatial coordinate   
    # cmos_main.register()
    # cmos=flap.get_data('W7X_ABES_CMOS', exp_id=shotID,
    #                 name="QSI", options={'Spatial calibration': True})   
    # timeindex = 10
    # dataplot = cmos.data[timeindex,:,:]
    # plt.tripcolor(cmos.get_coordinate_object("Device x").values.flatten(), cmos.get_coordinate_object("Device y").values.flatten(), dataplot.flatten())
    # plt.axis('equal')
    
    #testing adding the channel locations
    cmos_main.register()
    cmos=flap.get_data('W7X_ABES_CMOS', exp_id=shotID,
                    name="QSI", options={'Spatial calibration': True})   
    timeindex = 10
    dataplot = cmos.data[timeindex,:,:]
    # plt.tripcolor(cmos.get_coordinate_object("Device x").values.flatten(), cmos.get_coordinate_object("Device y").values.flatten(), dataplot.flatten())
    # plt.axis('equal')
    a = flap_w7x_abes.ShotSpatCal(shotID)
    a.generate_shotdata(options={'Plot': False, 'Overwrite': True})
    a.read()
    chimages = a.calc_chan_range()
    
    import numpy as np
    firstimage = 1
    image = []
    for key in chimages.keys():
        currimage = copy.deepcopy(chimages[key])
        currimage[np.where(chimages[key]<np.median(chimages[key])/2)] = np.nan
        currimage*= 255/np.nanmax(currimage)
        if firstimage == 1:
            firstimage = 0
        image += [currimage]
    image_backup = copy.deepcopy(image)
    print("Reading data finished")

    # calibimage = plt.imread("/DATA/repos/flap/modules/flap_w7x_abes/spatcal/2021/14_52baf8bfe_1636047267756.png")
    # plotting = np.zeros([calibimage.shape[0], calibimage.shape[1], 3])
    # plotting[:,:,1] = image[:calibimage.shape[0], :calibimage.shape[1]]/np.nanmax(image)*3
    # plotting[:,:,0] = cmos.data[timeindex,:,:]/np.max(cmos.data[timeindex,:,:])
    # plotting[:,:,2] = calibimage/np.max(calibimage)*2*0

    
    plt.tripcolor(cmos.get_coordinate_object("Device x").values.flatten(), cmos.get_coordinate_object("Device y").values.flatten(), (cmos.data[timeindex,:,:]-cmos.data[timeindex-1,:,:]/2-cmos.data[timeindex+1,:,:]/2).flatten(), cmap="gray")
    plt.axis('equal')
    # plotting[:,:,1][560:600,751:800]=0
    # plt.imshow(plotting[:,:,1])
    for currimage in image:
        currimage[560:600,751:800]=np.nan
        currimage[np.where(currimage<np.nanmax(currimage)/5)] = np.nan
        currimage[np.where(np.isnan(currimage) == False)]=1
        currimage[np.where(np.isnan(currimage))] = 0
        toplot = np.sign(np.abs(np.gradient(currimage)[0])+np.abs(np.gradient(currimage)[1]))
        toplot[np.where(toplot>0)] = 1
        toplot[np.where(toplot<1)] = np.nan
        plt.scatter(cmos.get_coordinate_object("Device x").values.flatten(), cmos.get_coordinate_object("Device y").values.flatten(), toplot.flatten(), color="tab:green", alpha=0.125)
        plt.show()
        plt.pause(0.1)
    plt.xlim([1.7,2.1])
    plt.ylim([5.7, 6.05])
    plt.axis("equal")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.title(f"Spatial data for {shotID} at {cmos.coordinate('Time')[0][10,0,0]}s")
    plt.tight_layout()
    plt.savefig("res.pdf")


