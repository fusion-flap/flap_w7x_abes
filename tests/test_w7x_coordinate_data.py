# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:59:59 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

flap.config.read(file_name="w7x_config.cfg")

def test_w7x_coordinate_data():
    print("\n------- test coordinate data  --------")
#    dp = 'c:/Data/W7-X_ABES/'
    d=flap.get_data('W7X_ABES',exp_id='20181018.003',name=['ABES-1[5,6]'],\
                    options={'Scaling':'Volt'},object_name='ABES',\
                    coordinates=flap.Coordinate(name='Time',c_range=[1,1.00001]))
    print("**** Storage contents")
    flap.list_data_objects()
    print("Getting Time for row 0")
    c,cl,ch = d.coordinate('Time',[...,0])
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    if (max(c.shape) < 30):
        print(" Coordinates: "+str(c))
    if (cl is not None):
        print("  L range shape:"+str(cl.shape)+" range: "+str(np.min(cl))+" - "+str(np.max(cl)))
    if (ch is not None):
        print("  H range shape:"+str(ch.shape)+" range: "+str(np.min(ch))+" - "+str(np.max(ch)))

    print("Getting Signal name for column 0")
    c,cl,ch = d.coordinate('Signal name',[0,...])
    if (max(c.shape) < 30):
        print(" Coordinates: "+str(c))

    print("Getting Channel for all data")
    c,cl,ch = d.coordinate('Channel',...)
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    if (max(c.shape) < 30):
        print(" Coordinates: "+str(c))

    print("Getting Time for all data")
    c,cl,ch = d.coordinate('Time',...)
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    if (max(c.shape) < 30):
        print(" Coordinates: "+str(c))

    d=flap.get_data('W7X_ABES',exp_id='20181018.003',name=['ABES-12'],\
                    options={'Scaling':'Volt'},object_name='ABES12',\
                    coordinates=flap.Coordinate(name='Time',c_range=[5,6]))
    print("Getting Time for all data of ABES12")
    c,cl,ch = d.coordinate('Time',...)
    print("  shape:"+str(c.shape)+" range: "+str(np.min(c))+" - "+str(np.max(c)))
    if (max(c.shape) < 30):
        print(" Coordinates: "+str(c))

test_w7x_coordinate_data()