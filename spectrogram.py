import flap
import flap_w7x_abes
flap_w7x_abes.register()
import numpy as np
import copy
import matplotlib.pyplot as plt

def rewrite_time_coordinate(dataobject):
    time=dataobject.get_coordinate_object('Time')
    step=np.mean(np.diff(time.values))
    start=time.values[0]
    dataobject.del_coordinate('Time')
    dataobject.add_coordinate_object(flap.Coordinate(name='Time',start=start, step=step, dimension_list=[0]))
    return dataobject

def rewrite_time_dimension_list(dataobject):
    time=dataobject.get_coordinate_object('Time')
    step=time.step
    start=time.start
    dataobject.del_coordinate('Time')
    dataobject.add_coordinate_object(flap.Coordinate(name='Time',start=start, step=step, dimension_list=[0]))
    return dataobject

def background_correction(datapath, exp_id, channel_name, timerange):
    d=flap.get_data('W7X_ABES', exp_id = exp_id, name = channel_name, coordinates = {'Time': timerange},
                    options={'Datapath':datapath,'Scaling':'Volt'})

    d_beam_on=flap.get_data('W7X_ABES', exp_id=exp_id, name='Chopper_time',
                            options={'Datapath':datapath,'State':{'Chop': 0, 'Defl': 0},'Start':1,'End':1},)
    d_beam_off=flap.get_data('W7X_ABES', exp_id=exp_id, name='Chopper_time',
                            options={'Datapath':datapath,'State':{'Chop': 1, 'Defl': 0},'Start':1,'End':1}) 
    d_on = d.slice_data(slicing={'Time':d_beam_on}, summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True})
    d_off = d.slice_data(slicing={'Time':d_beam_off}, summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True})
    d_off = d_off.slice_data(slicing={'Time':d_on},
                            options={'Interpolation':'Linear'})
    d_on_data = d_on.data
    d_off_data = d_off.data
    d_corr_data = d_on_data - d_off_data
    d_corr=copy.deepcopy(d_on)
    d_corr.data=d_corr_data

    d_off=rewrite_time_coordinate(d_off)
    d_on==rewrite_time_coordinate(d_on)
    d_corr=rewrite_time_coordinate(d_corr)

    return d_off, d_on, d_corr

def plot_spectrogram(d,window_size,time_step,freq_range,z_range):
    nperseg=window_size//d.get_coordinate_object('Time').step[0]
    noverlap=(window_size-time_step)//d.get_coordinate_object('Time').step[0]
    stft=flap.time_frequency_analysis._stft(d, coordinate='Time',
                                        options={'nperseg': nperseg, 
                                                 'noverlap':noverlap, 
                                                 'nfft': nperseg*2})
    dps = (stft.data.conjugate() * stft.data).real
    stft.data=np.log(dps)
    fig=plt.figure(figsize=(13,10))
    stft.plot(plot_type='image',options={'Log x': True, 'X range':freq_range,'Z range': z_range, 'Clear': True})

def get_beam_on_data(datapath, exp_id, channel_name, timerange):
    d=flap.get_data('W7X_ABES', exp_id = exp_id, name=channel_name, coordinates={'Time':timerange}, 
                    options={'Datapath':datapath,'Scaling':'Volt'})
    d_beam_on=flap.get_data('W7X_ABES', exp_id=exp_id, name='Chopper_time', coordinates={"Time":timerange},
                            options={'Datapath':datapath,'State':{'Chop': 0, 'Defl': 0},'Start':1,'End':1})
    d_on = d.slice_data(slicing={'Time':d_beam_on}, summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True})
    d_on=rewrite_time_coordinate(d_on)
    return d_on

if __name__ == '__main__':
    exp_id = '20241105.053'  
    datapath='/data2/W7-X/APDCAM' 
    timerange=[3,8]
    channel='ABES-14' 
    channel='ABES-*'
    background_corr=1
    window_size=0.02
    time_step=0.1
    freq_range=[100,1e4]
    z_range=[-20,-15]

    if background_corr:
        d_off, d_on, d_corr=background_correction(datapath=datapath, exp_id=exp_id, channel_name=channel, timerange=timerange)
        d_plot=d_corr
    else:
        d_plot=get_beam_on_data(datapath, exp_id, channel, timerange)

    try:
        d_plot.get_coordinate_object('Signal name')
        for ch in d_plot.get_coordinate_object('Signal name').values:
            d_plot_channel=d_plot.slice_data(slicing={'Signal name':ch})
            rewrite_time_dimension_list(d_plot_channel)
            plot_spectrogram(d_plot_channel,window_size,time_step,freq_range,z_range)
    except ValueError:
        plot_spectrogram(d_plot,window_size,time_step,freq_range,z_range)


