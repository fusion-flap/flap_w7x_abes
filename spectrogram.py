import flap
flap.config.read('flap_defaults.cfg')
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

def background_correction(exp_id, channel_name, timerange):
    d=flap.get_data('W7X_ABES', exp_id = exp_id, name = channel_name, coordinates = {'Time': timerange})

    d_beam_on=flap.get_data('W7X_ABES', exp_id=exp_id, name='Chopper_time',
                            options={'State':{'Chop': 0, 'Defl': 0},'Start':1,'End':1},)
    d_beam_off=flap.get_data('W7X_ABES', exp_id=exp_id, name='Chopper_time',
                            options={'State':{'Chop': 1, 'Defl': 0},'Start':1,'End':1}) 
    d_on = d.slice_data(slicing={'Time':d_beam_on},
                        summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True})
    d_off = d.slice_data(slicing={'Time':d_beam_off},
                        summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True})
    d_off = d_off.slice_data(slicing={'Time':d_on},
                            options={'Interpolation':'Linear'})
    d_on_data = d_on.data
    d_off_data = d_off.data
    ind = np.nonzero(np.logical_and(np.isfinite(d_off_data), np.isfinite(d_on_data)))[0]
    d_on_data = d_on_data[ind]
    d_off_data = d_off_data[ind]
    d_corr_data = d_on_data - d_off_data
    d_corr=copy.deepcopy(d_on)
    d_corr.data=d_corr_data

    d_off=rewrite_time_coordinate(d_off)
    d_on==rewrite_time_coordinate(d_on)
    d_corr=rewrite_time_coordinate(d_corr)

    return d_off, d_on, d_corr

def plot_spectrogram(d, name):
    nperseg=4096
    nfft=nperseg*2
    stft=flap.time_frequency_analysis._stft(d, 'Time',options={'nperseg': nperseg, 'noverlap':nperseg//2, 'nfft': nfft})
    dps = (stft.data.conjugate() * stft.data).real
    stft.data=np.log(dps)
    fig=plt.figure(figsize=(13,10))
    stft.plot(plot_type='image',options={'Log x': True, 'X range':[50,1e4],'Z range': [-20,-17], 'Clear': True})
    fig.savefig('spectrogram_'+name+'.png',dpi=300)

def get_beam_on_data(exp_id, channel_name):
    d=flap.get_data('W7X_ABES', exp_id = exp_id, name = channel_name)
    d_beam_on=flap.get_data('W7X_ABES', exp_id=exp_id, name='Chopper_time',
                            options={'State':{'Chop': 0, 'Defl': 0},'Start':0,'End':-3})
    d_on = d.slice_data(slicing={'Time':d_beam_on},
                        summing={'Rel. Sample in int(Time)':'Mean'},
                        options={'Regenerate':True})
    d_on==rewrite_time_coordinate(d_on)
    return d_on

exp_id = '20241121.010'   
timerange=[7,10]
for channel in ['ABES-'+str(a) for a in range(10,26)]:
    print('Working on', channel)
    d_off, d_on, d_corr=background_correction(exp_id=exp_id, channel_name=channel, timerange=timerange)
    plot_spectrogram(d_corr, channel)
