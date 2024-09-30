#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:15:16 2024

@author: zoletnik
"""
import matplotlib.pyplot as plt

import flap
import flap_w7x_abes
import numpy as np

flap_w7x_abes.register()

def beam_position_plot_apd(exp_id,time=[0,5]):
    blocks = {}
    blocks['1'] = {'L':'L1','R':'R1','Beam channels':[1,2,3]}  
    blocks['2'] = {'L':'L2','R':'R2','Beam channels':[7,8,9]}
    blocks['3'] = {'L':'L3','R':'R3','Beam channels':[13,14,15]}
    blocks['4'] = {'L':'L4','R':'R4','Beam channels':[19,20,21]}
    blocks['5'] = {'L':'L5','R':'R5','Beam channels':[25,26,27]}
    blocks['6'] = {'L':'L6','R':'R6','Beam channels':[32,33,34]}
    blocks['7'] = {'L':'L7','R':'R7','Beam channels':[38,39,40]}
    
    plt.close('all')
    plt.figure(figsize=(4,8))
    plt.suptitle(exp_id)
    ax = None
    t = None
    options = {'Resample':1e3}
    data_found = False
    err = None
    for ib in range(len(blocks)):
        maxarr = [0]
        try:
            d_l = flap.get_data('W7X_ABES',
                                 exp_id=exp_id,
                                 name=blocks[str(ib + 1)]['L'],
                                 coordinates={'Time':time},
                                 options=options
                                 )
            if (t is None):           
                t = d_l.coordinate('Time',options={'Change':True})[0]
            maxarr.append(np.max(d_l.data))    
            data_found = True
        except ValueError as e:
            err = e
            d_l = None        
        try:
            d_r = flap.get_data('W7X_ABES',
                                 exp_id=exp_id,
                                 name=blocks[str(ib + 1)]['R'],
                                 coordinates={'Time':time},
                                 options=options
                                 )
            if (t is None):           
                t = d_r.coordinate('Time',options={'Change':True})[0]
            maxarr.append(np.max(d_r.data))    
            data_found = True
        except ValueError  as e:
            d_r = None
            err = e         
        try:
            d_c = flap.get_data('W7X_ABES',
                             exp_id=exp_id,
                             name=['ABES-'+str(ch) for ch in blocks[str(ib + 1)]['Beam channels']],
                             coordinates={'Time':time},
                             options=options
                             ).slice_data(summing={'Channel':'Mean'})
            if (t is None):           
                t = d_c.coordinate('Time',options={'Change':True})[0]
            maxarr.append(np.max(d_c.data))    
            data_found = True
        except ValueError as e:
            d_c = None
            err = e
        if (not data_found):
            raise err

        m = max(maxarr)
        
        if (d_l is not None):
            if (ax is None):
                ax = plt.subplot(len(blocks),3,ib * 3 + 1)
            else:
                plt.subplot(len(blocks),3,ib * 3 + 1,sharex=ax)
            d_l.plot()
            plt.title('Block ' + str(ib+1) + 'L')
            plt.ylim(0,m)
        
        if (d_c is not None):
            if (ax is None):
                ax = plt.subplot(len(blocks),3,ib * 3 + 2)
            else:
                plt.subplot(len(blocks),3,ib * 3 + 2,sharex=ax)                
            d_c.plot()
            plt.title('Block '+str(ib+1))
            plt.ylim(0,m)
        
        if (d_r is not None):
            if (ax is None):
                ax = plt.subplot(len(blocks),3,ib * 3 + 3)
            else:
                plt.subplot(len(blocks),3,ib * 3 + 3,sharex=ax)                
            d_r.plot()
            plt.title('Block '+str(ib+1) + 'R')
            plt.ylim(0,m)
    plt.show()
    plt.pause(0.1)      
    plt.tight_layout()
        

plt.rcParams['figure.titlesize'] = 8
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.labelsize'] = 8  
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['ytick.minor.size'] = 1
plt.rcParams['legend.fontsize'] = 8         

#beam_position_plot_apd('20240917.032',time=[2,3])       

        
        
    
    
    