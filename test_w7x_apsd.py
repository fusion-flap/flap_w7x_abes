# -*- coding: utf-8 -*-
"""
Created on Thu May 30 11:25:33 2019

@author: Sandor Zoletnik
"""

def test_w7x_apsd():
    plt.close('all')
    flap.get_data('W7X_ABES',name='ABES-10',exp_id='20181018.008',object_name='ABES')
    d_beam_on=flap.get_data('W7X_ABES',exp_id='20181018.008',name='Chopper_time',
                     options={'State':{'Chop': 0, 'Defl': 0},'Start':1000,'End':-1000},\
                     object_name='Beam_on', coordinates={'Time':[4,6]})
    d_beam_off=flap.get_data('W7X_ABES',exp_id='20181018.008',name='Chopper_time',
                     options={'State':{'Chop': 1, 'Defl': 0},'Start':1000,'End':-1000},\
                     object_name='Beam_off', coordinates={'Time':[4,6]})
    legend = []
    flap.plot('ABES',axes='Time',options={'Y sep': 1})
    legend.append('ABES')
    flap.plot('Beam_on',axes=['Time',1],plot_type='scatter')
    legend.append('Beam on')
    flap.plot('Beam_off',axes=['Time',1.1],plot_type='scatter')
    legend.append('Beam off')
    plt.legend(legend)

    plt.figure()
    flap.apsd('ABES',output_name='ABES_APSD',intervals=d_beam_on,options={'Log':True,'Res':100,'Rang':[1e3,1e6]})
    flap.plot('ABES_APSD',options={'Log x':True, 'Log y': True, 'Y sep': 100})
    flap.apsd('ABES',output_name='ABES_APSD_off_Han',intervals=d_beam_off,options={'Han':True,'Log':True,'Res':100,'Rang':[1e3,1e6]})
    flap.apsd('ABES',output_name='ABES_APSD_off',intervals=d_beam_off,options={'Han':False,'Log':True,'Res':100,'Rang':[1e3,1e6]})
    flap.plot('ABES_APSD_off',options={'Log x':True, 'Log y': True, 'Y sep': 100})
    flap.plot('ABES_APSD_off_Han',options={'Log x':True, 'Log y': True, 'Y sep': 100})
    flap.list_data_objects()

test_w7x_apsd()