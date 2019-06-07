# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:59:59 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt
import numpy as np
import os

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def test_w7x_coordinate_data():
    print("\n------- test coordinate data with W7-X ABES module --------")
    print("**** Reading two channels 10 microsec signal")
    d=flap.get_data('W7X_ABES',exp_id='20181018.003',name=['ABES-1[5,6]'],\
                    options={'Scaling':'Volt'},object_name='ABES',\
                    coordinates=flap.Coordinate(name='Time',c_range=[1,1.00001]))
    print("**** Storage contents")
    flap.list_data_objects()
    print("**** Getting Time for row 0")
    c,cl,ch = d.coordinate('Time',[...,0])
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    print(" Coordinates: "+str(c))
    if (cl is not None):
        print("  Low range shape:"+str(cl.shape)+" range: "+str(np.min(cl))+" - "+str(np.max(cl)))
    if (ch is not None):
        print("  High range shape:"+str(ch.shape)+" range: "+str(np.min(ch))+" - "+str(np.max(ch)))

    print("***** Getting Signal name for column 0")
    c,cl,ch = d.coordinate('Signal name',[0,...])
    if (max(c.shape) < 30):
        print(" Coordinates: "+str(c))

    print("***** Getting Channel for all data")
    c,cl,ch = d.coordinate('Channel',...)
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    print(" Coordinates: "+str(c))

    print("***** Getting Time for all data")
    c,cl,ch = d.coordinate('Time',...)
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    print(" Coordinates: "+str(c))

# Reading configuration file in the test directory
thisdir = os.path.dirname(os.path.realpath(__file__))
fn = os.path.join(thisdir,"w7x_config.cfg")
flap.config.read(file_name=fn)

test_w7x_coordinate_data()