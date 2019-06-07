# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:54:29 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt
import os

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def test_w7x_chopper(test_type):
    
    plt.figure()
    flap.delete_data_object('*')
    # Camera chopper
    if (test_type == 'timed'):
        print("***** Testing timed chopper")
        # timed chopper
        exp_id = '20181017.024'
        timerange = [3,3.001]
    elif (test_type == 'camera'):
        # Camera chopper
        print("***** Testing camera chopper")
        exp_id = '20181018.008'
        timerange = [4,6]
    else:
        raise ValueError("Invalid test type.")

    print("***** Reading beam-on, beam-off intervals into DataObjects.")
    d_beam_on=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'State':{'Chop': 0, 'Defl': 0}},\
                         object_name='Beam_on',
                         coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'State':{'Chop': 1, 'Defl': 0}},\
                         object_name='Beam_off',
                         coordinates={'Time': timerange})

    if (test_type == 'camera'):
        print("***** Reading beam deflection intervals into a DataObject.")
        d_defl_on=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 1}},\
                             object_name='Beam_defl',
                             coordinates={'Time': timerange})
    print("***** Reading the full chopper period timing into a DataObject.")
    d_chop_period=flap.get_data('W7X_ABES',exp_id=exp_id,name='Chopper_time',
                         options={'Phase': ..., 'Start delay':0},\
                         object_name='Chop_period',
                         coordinates={'Time': timerange})
    print("***** Reading ABES-20")
    d=flap.get_data('W7X_ABES',exp_id=exp_id,name=['ABES-20'],
                    options={'Scaling':'Volt'},\
                    coordinates={'Time': timerange},
                    object_name='ABES-20')
    legend = []
    print("***** Plotting ABES-20 as a function of time.")
    flap.plot('ABES-20',axes=['Time'])
    legend.append('ABES-20')
    print("***** Overplotting different time intervals.")
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
    print("***** Plotting ABES-20 as a function of sample number.")
    flap.plot('ABES-20',axes=['Sample'])
    legend.append('ABES-20')
    flap.plot('Beam_on',axes=['Sample',1.5],plot_type='scatter',options={'Force':True})
    legend.append('Beam on')

    print("***** Calculating and plotting mean data in beam on and beam off times.")
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

    print("***** Plotting minimum and maximum signal over chopper periods as a function of time in chopper period.")
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

# Reading configuration file in the test directory
thisdir = os.path.dirname(os.path.realpath(__file__))
fn = os.path.join(thisdir,"w7x_config.cfg")
flap.config.read(file_name=fn)


plt.close('all')
test_w7x_chopper('camera')
test_w7x_chopper('timed')