# -*- coding: utf-8 -*-
"""
Created on Thu May 30 08:06:35 2019

@author: Sandor Zoletnik

Test chopper selection.

"""

import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def test_w7x_select_multislice():
    plt.close('all')
    flap.delete_data_object('*',exp_id='*')
    exp_id = '20181018.008'
    # Reading data
    d=flap.get_data('W7X_ABES',name='ABES-[8-20]',
                    exp_id=exp_id,
                    options={'Scaling':'Volt'},
                    object_name='ABES')
    # Reading beam on times
    d_beam_on=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'State':{'Chop': 0, 'Defl': 0}}, object_name='Beam_on')
    # Plotting the mean data in chopper intervals as a function of time at the beginning of the chopper intervals
    flap.plot('ABES',axes='Start Time in int(Sample)',
              slicing={'Sample':d_beam_on},
              summing={'Interval(Sample) sample index':'Mean'},
              options={'Y sep':3})
    # Reading beam off times
    d_beam_off=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'State':{'Chop': 1, 'Defl': 0}}, object_name='Beam_off')
    # Plotting the mean data in beam off intervals as a function of time at the beginning of the chopper intervals
    flap.plot('ABES',axes='Start Time in int(Sample)',
              slicing={'Sample':d_beam_off},
              summing={'Interval(Sample) sample index':'Mean'},
              options={'Y sep':3})
    flap.list_data_objects()

test_w7x_select_multislice()