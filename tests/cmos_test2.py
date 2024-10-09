#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:13:29 2024

@author: apdcam
"""

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
import numpy as np

flap_w7x_abes.register()

if __name__ == "__main__":
    shotID = "20230316.089"
    shotID1="T20240806.018"
    shotID2 = '20230307.059'
    flap.config.read()

    # testing removing the background
    cmos1 = cmos_main.w7x_abes_cmos_get_data(exp_id=shotID1)
    cmos1_on = cmos1.get_chopstate()
    cmos1_off = cmos1.get_chopstate(chop=1)
    plt.figure()
    for timeindex in range(cmos_on.data.shape[0]):
        plt.clf()
        plt.imshow(cmos1_on.data[timeindex,:,:]-cmos1_off.data[timeindex,:,:])
        plt.title(timeindex)
        plt.show()
        plt.pause(0.1)
        
    timeindex_rel = 20
    plt.imshow(cmos1_on.data[timeindex_rel,:,:]-cmos1_off.data[timeindex_rel,:,:])
    plt.title(timeindex_rel)
    plt.figure()
        
    # testing removing the background
    cmos2 = cmos_main.w7x_abes_cmos_get_data(exp_id=shotID2)
    cmos2_on = cmos2.get_chopstate()
    cmos2_off = cmos2.get_chopstate(chop=1)
    plt.figure()
    for timeindex in range(cmos_on.data.shape[0]):
        plt.clf()
        plt.imshow(cmos2_on.data[timeindex,:,:]-cmos2_off.data[timeindex,:,:])
        plt.title(timeindex)
        plt.show()
        plt.pause(0.1)
    
    # testing adding the spatial coordinate   
    cmos_main.register()
    cmos=flap.get_data('W7X_ABES_CMOS', exp_id=shotID1,
                    name="QSI", options={'Spatial calibration': True})   
    timeindex = 20
    dataplot = cmos.data[timeindex,:,:]
    plt.tripcolor(cmos.get_coordinate_object("Device x").values.flatten(), cmos.get_coordinate_object("Device y").values.flatten(), dataplot.flatten())
    plt.axis('equal')
    
    #testing adding the channel l    cmos = cmos_main.w7x_abes_cmos_get_data(exp_id=shotID)
    cmos_on = cmos.get_chopstate()
    cmos_off = cmos.get_chopstate(chop=1)
    for timeindex in range(cmos_on.data.shape[0]):
        plt.clf()
        plt.imshow(cmos_on.data[timeindex,:,:]-cmos_off.data[timeindex,:,:])
        plt.show()
        plt.pause(0.1)

    plt.tripcolor(cmos.get_coordinate_object("Device x").values.flatten(), cmos.get_coordinate_object("Device y").values.flatten(), cmos.data[timeindex,:,:].flatten(), cmap="gray")
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
        coord_x = cmos.get_coordinate_object("Device x").values[np.where(np.isnan(toplot) == False)].flatten()
        coord_y = cmos.get_coordinate_object("Device y").values[np.where(np.isnan(toplot) == False)].flatten()
        plt.scatter(coord_x,coord_y, toplot[np.where(np.isnan(toplot) == False)].flatten(), color="tab:green", alpha=0.125)
        plt.show()
        plt.pause(0.1)
    plt.xlim([1.7,2.1])
    plt.ylim([5.7, 6.05])
    plt.axis("equal")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.title(f"Spatial data for {shotID} at {cmos.coordinate('Time')[0][10,0,0]}s")
    plt.tight_layout()
    # plt.savefig("res.pdf")
    
    #obtaining the light profile at the channels
    cmos_on = cmos.get_chopstate()
    cmos_off = cmos.get_chopstate(chop=1)
    cmos_on.data -= cmos_off.data
    light_profile=cmos_on.get_channel_signal()
    light = np.zeros([light_profile.data.shape[0], 40])
    light_l = np.zeros([light_profile.data.shape[0], 8])
    light_r = np.zeros([light_profile.data.shape[0], 8])

    mainindex = 0
    for channel in light_profile.get_coordinate_object('Channel').values:
        mainindex += 1
        try:
            index = int(channel)-1
            light[:, index] = light_profile.data[:, mainindex-1]
        except:
            if 'L' in channel:
                index = int(channel[1:])-1
                light_l[:, index] = light_profile.data[:, mainindex-1]
            else:
                index = int(channel[1:])-1
                light_r[:, index] = light_profile.data[:, mainindex-1]
    # bg = light_r/2+light_l/2
    # rel_bg_channel = (np.arange(40)+1)//6
    # for channel in range(40):
    #     light[:,channel] -= bg[:,rel_bg_channel[channel]]
    # light = light-np.min(cmos.data)
    
    timeindex=10
    fig=plt.figure(figsize=[6.5,3])
    ax = fig.add_subplot(1,1,1)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax.set_ylabel('Signal [a.u.]')
    ax.yaxis.set_label_coords(-0.125,0.5)
    plt.subplot(2,1,1)
    plt.plot(np.arange(40)+1,light[timeindex,:], color="tab:green", label=f"Time: {cmos.coordinate('Time')[0][timeindex,0,0]}s")
    plt.legend()
    # plt.ylabel('Signal [a.u.]')
    plt.subplot(2,1,2)
    plt.plot(np.arange(40)+1,np.mean(light[-6:,:],axis=0), color="tab:green", label="During gas injection")
    exp_params = np.polyfit(np.arange(40)+1, np.log(np.mean(light[-6:,:],axis=0)), 1)
    plt.plot(np.arange(40)+1, np.exp(exp_params[1]+(np.arange(40)+1)*exp_params[0]), ls="--", label="Fitted exponential")
    plt.xlabel('Channel number')
    # plt.ylabel('Signal [a.u.]')
    plt.suptitle(f"Light intensity of CMOS image at the channel locations\n for {shotID}")
    plt.legend()
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

