# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:50:28 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

flap.config.read(file_name="w7x_config.cfg")

def test_W7X_data():
    print("\n------- test data read with W7-X ABES data --------")
#    dp = 'c:/Data/W7-X_ABES/'
    d=flap.get_data('W7X_ABES',exp_id='20181018.008',name=['ABES-[1-31]','ABES-[33-40]'],options={'Scaling':'Volt'},\
                    object_name='ABES',coordinates={'Time':[3,4]})
    print("**** Storage contents")
#    flap.list_data_objects()
#    flap.slice_data('ABES',slicing={'Time':flap.Intervals(3.1,3.11)},summing={'Time':'Mean'})
#    flap.list_data_objects()
    plt.close('all')
    flap.plot('ABES',slicing={'Signal name':'ABES-20'},axes='Time')
    plt.figure()
    flap.plot('ABES',slicing={'Time':flap.Intervals(3.26,3.28)},summing={'Time':'Mean'},axes='Channel')

test_W7X_data()