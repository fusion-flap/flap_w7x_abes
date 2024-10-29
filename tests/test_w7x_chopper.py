# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:54:29 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import matplotlib.pyplot as plt
import os
import numpy as np

import flap
import flap_w7x_abes

flap_w7x_abes.register()
    
def test_w7x_chopper(test_type,exp_id=None,timerange=None):

    plt.figure()
    flap.delete_data_object('*')
    # Camera chopper
    if (test_type == 'timed'):
        print("***** Testing timed chopper")
        # timed chopper
        if (exp_id is None):
            exp_id = '20181017.024'
        if (timerange is None):
            timerange = [3,3.001]
    elif (test_type == 'camera'):
        # Camera chopper
        print("***** Testing camera chopper")
        if (exp_id is None):
            exp_id = '20181018.008'
        if (timerange is None):
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

def test_chopper_timing(exp_id=None, timerange=None,signal='ABES-15',resample=1e3,x_axis='Time',start_shift=None,end_shift=None,
                        beam_on_start_delay=None,beam_on_end_delay=None,beam_off_start_delay=None,beam_off_end_delay=None):
    """
    Test the chopper timing visually. Plots a signal and the chopper on and off periods.
    Also calculates the mean signal in on/off intervals and plots them.

    Parameters
    ----------
    exp_id : string
        The experiment ID.
    timerange : 2-element list, optional
        The time range. The default is None, which means all data.
    signal : string, optional
        The signal name. The default is 'ABES-15'.
    resample : float, optional
        The resampling during data read. The default is 1e3. If None no resampling is done.
    x_axis : string, optional
        'Time' or 'Sample'. The default is 'Time'.
    start_shift : float
        The shift of the chopper start points [microsec]
    end_shift : float
        The shift of the chopper end points [microsec]

    Returns
    -------
    None.

    """
    if (resample is not None):
        options = {'Resample':resample}
    else:
        options = None
    d=flap.get_data('W7X_ABES',
                    exp_id = exp_id,
                    name = signal,
                    options = options,
                    coordinates = {'Time': timerange}
                    )
    chopper_mode,beam_on_time,beam_off_time,period_time,d_beam_on,d_beam_off = flap_w7x_abes.chopper_parameters(exp_id,
                                                                                                                timerange=timerange,
                                                                                                                beam_on_start_delay=beam_on_start_delay,
                                                                                                                beam_on_end_delay=beam_on_end_delay,
                                                                                                                beam_off_start_delay=beam_off_start_delay,
                                                                                                                beam_off_end_delay=beam_off_end_delay
                                                                                                                )
    print('Chopper mode: '+ chopper_mode,flush=True)
    on1,on2,on3 = d_beam_on.coordinate('Time')
    off1,off2,off3 = d_beam_on.coordinate('Time')
    print('Beam on {:7.3f} [ms], Beam off {:7.3f} [ms]'.format((on3[0]-on2[0]) * 1000,(off3[0]-off2[0]) * 1000))
    d.plot(axes=x_axis,options={'All':True})
    d_beam_on.plot(plot_type='scatter',axes=[x_axis,np.min(d.data)],options={'Force':True})
    d_beam_off.plot(plot_type='scatter',axes=[x_axis,np.min(d.data)],options={'Force':True})

    d_on = d.slice_data(slicing={'Time':d_beam_on},
                        summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True}
                        )
    d_off = d.slice_data(slicing={'Time':d_beam_off},
                         summing={'Rel. Sample in int(Time)':'Mean'},
                         options={'Regenerate':True}
                         )
    d_on.plot(axes=x_axis,plot_options={'marker':'o'})
    d_off.plot(axes=x_axis,plot_options={'marker':'s'})
    plt.legend([signal,'Beam on','Beam off','Mean signal on','Mean signal off'])


#plt.close('all')
#test_w7x_chopper('camera',exp_id='20181018.003')
#test_w7x_chopper('timed',exp_id='20181018.003')