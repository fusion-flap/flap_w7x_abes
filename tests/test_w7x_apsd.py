# -*- coding: utf-8 -*-
"""
Created on Thu May 30 11:25:33 2019

@author: Sandor Zoletnik
"""
import matplotlib.pyplot as plt
import os

def test_w7x_apsd():
    plt.close('all')
    flap.delete_data_object('*')
    print("***** Reading ABES-10, beam-on and beam-off times.")
    flap.get_data('W7X_ABES',name='ABES-10',exp_id='20181018.008',object_name='ABES')
    d_beam_on=flap.get_data('W7X_ABES',exp_id='20181018.008',name='Chopper_time',
                     options={'State':{'Chop': 0, 'Defl': 0},'Start':1000,'End':-1000},\
                     object_name='Beam_on', coordinates={'Time':[4,6]})
    d_beam_off=flap.get_data('W7X_ABES',exp_id='20181018.008',name='Chopper_time',
                     options={'State':{'Chop': 1, 'Defl': 0},'Start':1000,'End':-1000},\
                     object_name='Beam_off', coordinates={'Time':[4,6]})
    print("***** Plotting signal and time intervals.")
    legend = []
    flap.plot('ABES',axes='Time',)
    legend.append('ABES')
    flap.plot('Beam_on',axes=['Time',1],plot_type='scatter',options={'Force':True})
    legend.append('Beam on')
    flap.plot('Beam_off',axes=['Time',1.1],plot_type='scatter',options={'Force':True})
    legend.append('Beam off')
    plt.legend(legend)

    
    print("***** Calculating power spectra with logarithmic freqeuency binning for beam-on, beam-off times.")
    flap.apsd('ABES',output_name='ABES_APSD',intervals=d_beam_on,options={'Log':True,'Res':100,'Rang':[1e3,1e6]})
    flap.apsd('ABES',output_name='ABES_APSD_off',intervals=d_beam_off,options={'Log':True,'Res':100,'Rang':[1e3,1e6]})
    print("***** Plotting spectra.")
    plt.figure()
    flap.plot('ABES_APSD',options={'Log x':True, 'Log y': True, 'Y sep': 100})
    flap.plot('ABES_APSD_off',options={'Log x':True, 'Log y': True, 'Y sep': 100})
    legend = ['Beam on', 'Beam off']
    plt.legend
    flap.list_data_objects()

# Reading configuration file in the test directory
thisdir = os.path.dirname(os.path.realpath(__file__))
fn = os.path.join(thisdir,"w7x_config.cfg")
flap.config.read(file_name=fn)

test_w7x_apsd()