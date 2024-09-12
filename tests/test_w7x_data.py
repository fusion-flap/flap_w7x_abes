# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:50:28 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt
import os
import time

import flap
import flap_w7x_abes


def test_W7X_data():
    flap.delete_data_object('*')
    print("\n------- test data read with W7-X ABES data --------")
    print("***** Reading 39 channels (1-40, omitting 32) in time interval 3-4 s.")
    d=flap.get_data('W7X_ABES',exp_id='20181018.003',name=['ABES-[1-31]','ABES-[33-40]'],\
                    object_name='ABES',coordinates={'Time':[3,4]})
    print("**** Storage contents")
    flap.list_data_objects()
    plt.close('all')
    print("**** Plotting ABES-20")
    flap.plot('ABES',slicing={'Signal name':'ABES-20'},axes='Time')
    plt.figure()
    print("**** Plotting mean signal in 3.26-3.28s time interval.")
    flap.plot('ABES',slicing={'Time':flap.Intervals(3.26,3.28)},summing={'Time':'Mean'},axes='Channel')


# Reading configuration file in the test directory
thisdir = os.path.dirname(os.path.realpath(__file__))
fn = os.path.join(thisdir,"w7x_config.cfg")
flap.config.read(file_name=fn)

test_W7X_data()

import numpy as np
mod = np.zeros(d.data.shape[0])
for time in range(d.data.shape[0]):
    mod[time] = np.var(d.data[time,3,400:450])
plt.plot(d.coordinate("Time")[0][:,0,0], mod)