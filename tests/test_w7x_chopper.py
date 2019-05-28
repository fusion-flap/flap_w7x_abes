# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:54:29 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt

import flap
import flap_w7x_abes

flap_w7x_abes.register()

flap.config.read(file_name="w7x_config.cfg")

def test_w7x_chopper(test_type):
    plt.figure()
    flap.delete_data_object('*')
    # Camera chopper
    if (test_type == 'timed'):
        # timed chopper
        exp_id = '20181017.024'
        timerange = [3,3.001]
    elif (test_type == 'camera'):
        # Camera chopper
        exp_id = '20181018.008'
        timerange = [4,6]
    else:
        raise ValueError("Invalid test type.")

    d_beam_on=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'State':{'Chop': 0, 'Defl': 0}},\
                         object_name='Beam_on',
                         coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'State':{'Chop': 1, 'Defl': 0}},\
                         object_name='Beam_off',
                         coordinates={'Time': timerange})

    if (test_type == 'camera'):
        d_defl_on=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 1}},\
                             object_name='Beam_defl',
                             coordinates={'Time': timerange})

    d_chop_period=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'Phase': ..., 'Start delay':0},\
                         object_name='Chop_period',
                         coordinates={'Time': timerange})
    d=flap.get_data('W7X_ABES',exp_id=exp_id,name=['ABES-20'],
                    options={'Scaling':'Volt'},\
                    coordinates={'Time': timerange},
                    object_name='ABES-20')
    legend = []
    flap.plot('ABES-20',axes=['Time'])
    legend.append('ABES-20')
    flap.plot('Beam_on',axes=['Time',1.25],plot_type='scatter',options={'Force':True})
    legend.append('Beam on')
    flap.plot('Beam_off',axes=['Time',1.3],plot_type='scatter',options={'Force':True})
    legend.append('Beam off')
    if (test_type == 'camera'):
        flap.plot('Beam_defl',axes=['Time',1.35],plot_type='scatter',options={'Force':True})
        legend.append('Beam deflected')
    flap.plot('Chop_period',axes=['Time',1.4],plot_type='scatter',options={'Force':True})
    legend.append('Chopper period')
    plt.legend(legend)

    plt.figure()
    flap.plot('ABES-20',axes=['Sample'])
    legend.append('ABES-20')
    flap.plot('Beam_on',axes=['Sample',1.25],plot_type='scatter',options={'Force':True})
    legend.append('Beam on')

    flap.slice_data('ABES-20',slicing={'Sample':d_beam_on},
                    summing={'Interval(Sample) sample index':'Mean'}, output_name='ABES-20_on')
    flap.slice_data('ABES-20',slicing={'Sample':d_beam_off},
                    summing={'Interval(Sample) sample index':'Mean'}, output_name='ABES-20_off')
    flap.list_data_objects()

    plt.figure()
    flap.plot('ABES-20_on',axes='Start Time in int(Sample)',plot_type='scatter')
    flap.plot('ABES-20_off',axes='Start Time in int(Sample)',plot_type='scatter')
    legend = ['ABES-20 beam on','ABES-20 beam off']
    plt.legend(legend)

    flap.slice_data('ABES-20',
              slicing={'Sample':d_chop_period},
              output_name='ABES_20_chop_slice')
    flap.list_data_objects()
    plt.figure()
    flap.plot('ABES-20',
              slicing={'Sample':d_chop_period},
              summing={'Interval(Sample)':'Mean'},axes='Rel. Time in int(Sample)')
    flap.plot('ABES-20',
              slicing={'Sample':d_chop_period},
              summing={'Interval(Sample)':'Min'},axes='Rel. Time in int(Sample)')
    flap.plot('ABES-20',
               slicing={'Sample':d_chop_period},
               summing={'Interval(Sample)':'Max'},axes='Rel. Time in int(Sample)')
    legend = ['Mean','Minimum','Maximum']
    plt.legend(legend)

plt.close('all')
test_w7x_chopper('camera')
test_w7x_chopper('timed')