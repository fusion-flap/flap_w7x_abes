# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 10:17:19 2025

@author: zoletnik
"""

import matplotlib.pyplot as plt
import numpy as np
import math

import flap
import flap_w7x_abes

flap_w7x_abes.register() 

def cxrs_summary(exp_id,linewidth=0.05,test_spectrum=False,test_corr=False,integration_time=1,c1_limit=-0.2,c2_limit=0.7,list_result=True):
    """
    Finds lines in the ABES CXRS spectrum integrated for all fibres and times. Calculates the time evolulution of 
    the line intensities in each fibre and line. Selects the lines, fibres and time intervals where the intensity
    is modulated in correlation with beam modulation. 

    Parameters
    ----------
    exp_id : str
        Experiment ID.
    linewidth : float, optional
        The smooth lengt in nm for finding lines. 
        Twice this value will be used for integrating the line intensity. The default is 0.05 nm.
    test_spectrum : bool, optional
        Plot the reusult of line selection. The default is False.
    test_corr : bool, optional
        Plot spectra around lines and the correlations vs time. The default is False.
    list_result : bool, optional
        List the result.
    integration_time : float, optional
        The integration time in the discharge for processing data. The default is 1.
    c1_limit : float, optional
        The beam modulation is accepted if the correlation at 1 frame is below this limit. The default is -0.2.
    c2_limit : float, optional
        The beam modulation is accepted if the correlation at 2 frames is above this limit. The default is 0.7.

    Returns
    -------
    active_line_list: list
        A list of dictinoaries, each describing one case when the correlation is good enough.

    """
       
    d = flap.get_data('W7X_ABES_CXRS', exp_id=exp_id,name="QSI_CXRS")
#    flap.list_data_objects(d)
#    print(d.config)
    # Summing for channels and fibres to get the lines
    d_sum = d.slice_data(summing={'Time':'Mean','Channel':'Mean'})
    wavelength = d_sum.coordinate('Wavelength')[0]
    wres = abs(wavelength[1] - wavelength[0])
    w_line_global,a_line_global,offset,noiselevel = find_lines(wavelength,d_sum.data,test=test_spectrum,linewidth=linewidth,auto_offset=True)
    w_interval_start = [w - linewidth for w in w_line_global]
    w_interval_stop = [w + linewidth for w in w_line_global]
    t = d.coordinate('Time',options={'Change only':True})[0].flatten()
    tres = t[1] - t[0]
    optical_channels = d.coordinate('Optical channel',options={'Change only':True})[0].flatten()
    channels = d.coordinate('Channel',options={'Change only':True})[0].flatten()
    R = d.coordinate('Device R',options={'Change only':True})[0].flatten()
    X = d.coordinate('Device x',options={'Change only':True})[0].flatten()
    Y = d.coordinate('Device y',options={'Change only':True})[0].flatten()
    line_data = np.zeros((len(t),len(w_line_global),len(optical_channels)))
    active_line_list = []
    kernel_length = int(round(integration_time / (t[1] - t[0]))) // 2 * 2 + 1
    kernel = np.ones(kernel_length) / kernel_length
    if (test_corr):
        plt.figure()
        spectrum_plot = plt.gcf()
        plt.figure()
        corr_plot = plt.gcf()
    for i in range(len(w_line_global)):
        for i_ch in range(len(optical_channels)):  
            if (test_corr):
                plt.figure(spectrum_plot.number)
                plt.clf()
                # The wavelength scale is not exact on this plot. This is a flap.plot problem.
#                d.plot(slicing={'Optical channel':optical_channels[i_ch]},axes=['Wavelength','Time'],plot_type='image')
                d.plot(slicing={'Optical channel':optical_channels[i_ch],'Time':5},axes=['Wavelength'])
                plt.plot([w_interval_start[i]]*2,plt.ylim())
                plt.plot([w_interval_stop[i]]*2,plt.ylim())
                plt.xlim(w_line_global[i]-1,w_line_global[i]+1)
                plt.title("Channel: {:d}".format(channels[i_ch]))
                plt.pause(0.1)
            line_data[:,i,i_ch] = d.slice_data(slicing={'Optical channel':optical_channels[i_ch],
                                                        'Wavelength':flap.Intervals(w_interval_start[i],w_interval_stop[i])},
                                               summing={'Wavelength':'Mean'}
                                               ).data
            line_data[:,i,i_ch] -= line_data[0,i,i_ch]
            mean_smooth_0 = np.convolve(line_data[:-2,i,i_ch],kernel,mode='valid')
            mean_smooth_1 = np.convolve(line_data[1:-1,i,i_ch],kernel,mode='valid') 
            mean_smooth_2 = np.convolve(line_data[2:,i,i_ch],kernel,mode='valid')
            norm_0 = np.sqrt(np.convolve((line_data[kernel_length // 2 : -(kernel_length // 2) - 2,i,i_ch] - mean_smooth_0) ** 2,kernel,mode='valid'))
            norm_1 = np.sqrt(np.convolve((line_data[kernel_length // 2 + 1 : -(kernel_length // 2) - 2 + 1,i,i_ch] - mean_smooth_1) ** 2,kernel,mode='valid'))
            norm_2 = np.sqrt(np.convolve((line_data[kernel_length // 2 + 2 : -(kernel_length // 2) - 2 + 2,i,i_ch] - mean_smooth_2) ** 2,kernel,mode='valid'))
            c_0 = np.convolve((line_data[kernel_length // 2 : -(kernel_length // 2) - 2,i,i_ch] - mean_smooth_0) \
                              * (line_data[kernel_length // 2 : -(kernel_length // 2) - 2,i,i_ch] - mean_smooth_0),kernel,mode='valid')         
            c_0 /= norm_0 * norm_0
            c_1 = np.convolve((line_data[kernel_length // 2 : -(kernel_length // 2) - 2,i,i_ch] - mean_smooth_0) \
                              * (line_data[kernel_length // 2 + 1: -(kernel_length // 2) - 2 + 1,i,i_ch] - mean_smooth_1),kernel,mode='valid')   
            c_1 /= norm_0 * norm_1
            c_2 = np.convolve((line_data[kernel_length // 2 : -(kernel_length // 2) - 2,i,i_ch] - mean_smooth_0) \
                        * (line_data[kernel_length // 2 + 2: -(kernel_length // 2) - 2 + 2,i,i_ch] - mean_smooth_2),kernel,mode='valid')   
            c_2 /= norm_0 * norm_2
            t_corr = np.convolve(np.convolve(t[:-2],kernel,mode='valid'),kernel,mode='valid')
            if (test_corr):
                plt.figure(corr_plot.number)
                plt.clf()
                plt.plot(t_corr,c_0)
                plt.plot(t_corr,c_1)
                plt.plot(t_corr,c_2)
                plt.legend(['c0','c1','c2'])
                plt.title("W: {:5.1f}, ch: {:s}".format(w_line_global[i],optical_channels[i_ch]))
                plt.pause(0.05)
                
             
            ind = np.nonzero(np.logical_and(c_2 > c2_limit,c_1 < c1_limit))[0]
            if (len(ind) != 0):  
                start_ind = ind[0]
                while start_ind < ind[-1]:
                    diff_ind = np.diff(ind[start_ind:])
                    ind_diff = np.nonzero(diff_ind > 1)[0]
                    if (len(ind_diff) == 0):
                        stop_ind = ind[-1]
                    else:
                        stop_ind = ind[ind_diff[0] - 1]
                    active_line_list.append({'w':w_line_global[i],
                                             'och':optical_channels[i_ch],
                                             'ch' : channels[i_ch],
                                             'trange':[t_corr[start_ind],t_corr[stop_ind]],
                                             'c_1':c_1,
                                             'c_2':c_2,
                                             'c_time': t_corr,
                                             'amp': line_data[:,i,i_ch],
                                             't': t,
                                             'R': R[i_ch],
                                             'x': X[i_ch],
                                             'y': Y[i_ch]
                                             }
                                            )  
                    start_ind = stop_ind + 1
    if (list_result):
        for l in active_line_list:
            print("Wavelength: {:5.1f}nm, Optical Ch: {:s}, time range: [{:4.1f},{:4.1f}]".format(l['w'],l['och'],*l['trange']))                
    return active_line_list


def test_proc(exp_id):
    d = flap.get_data('W7X_ABES_CXRS', exp_id=exp_id,name="QSI_CXRS")
    dd = d.slice_data(slicing={'Time':2,'Channel':10})
    w = dd.coordinate('Wavelength')[0]
    w,a = find_lines(w,dd.data,test=True)
    for i in range(len(w)):
        print("{:5.1f}[nm]: {:7.1f}".format(w[i],a[i]))
    
def find_lines(wavelength,data,test=True,linewidth=0.05,fit_order=None,auto_offset=True,min_noiselevel=5):
    """
    Finds spectral lines in a spectrum. 
    The expected FWHM linewidth is used as inpupt parameter, 
    will ignore too wide (>10xlinewidth) and too narrow (<2xlinewidth) lines.

    Parameters
    ----------
    wavelength : numpy array
        The wavelength scale.
    data : numpy array
        The spectrum.
    test : bool, optional
        If True will make test plot. The default is True.
    linewidth : float, optional
        Expected FWHM line width in nm. The default is 0.4.
    fit order : None or int
        If not None will fit a polinomial with fit_order to the whole spectrum and subtract.

    Returns
    -------
    linelist: list
        List of line wavelength [nm]
    line_amp_list: list
        List of line peak amplitudes.
    auto_offset: float or None
        The automatically determined offset
    noiselevel: float
        The noise level. 

    """
    
    if (test):
        plt.figure()
        plt.plot(wavelength,data)
    w_res = abs(wavelength[1] - wavelength[0])
    n_kernel = int(2 * round(linewidth / w_res)) // 2 * 2 + 1
    kernel = np.exp(-((np.arange(n_kernel) - (n_kernel - 1) / 2) * w_res) ** 2 / linewidth ** 2)
    data_smooth = np.convolve(data,kernel,mode='same') / np.sum(kernel)
    data_smooth = data_smooth[n_kernel : - n_kernel]
    wavelength_smooth = wavelength[n_kernel: -n_kernel]
    if (test):
        plt.plot(wavelength_smooth,data_smooth)
    if (fit_order is not None):
        p = np.polyfit(wavelength_smooth,data_smooth,2)
        fitdata = p[-1]
        for i in range(fit_order):
            fitdata += p[i] * wavelength_smooth ** (fit_order - i)
            if (test):
                plt.plot(wavelength_smooth,fitdata) 
        data_proc = data_smooth - fitdata
    else:
        data_proc = data_smooth
        fitdata = 0
    binsize = (np.max(data_proc) -  np.min(data_proc)) / 10    
    while True:
        bin_num = int((np.max(data_proc) -  np.min(data_proc)) / binsize)
        hist, bin_edges = np.histogram(abs(data_proc),bins=bin_num)
        if (np.max(hist) < len(data_proc) / 50 ):
            binsize *= 2
            continue
        if (np.max(hist) > len(data_proc) / 10  ):
            binsize /= 2
            continue
        ind_max = np.argmax(hist)
        ind = np.nonzero(hist[ind_max:] < hist[ind_max] / 10)[0]
        noiselevel = bin_edges[ind_max + ind[0]]
        if (auto_offset):
            fitdata = np.mean(bin_edges[ind_max:ind_max+2])
            data_proc -= fitdata
            noiselevel -= fitdata
            auto_offset_value  = fitdata
        else:
            auto_offset_value = None
        break
    if (noiselevel < min_noiselevel):
        noiselevel = min_noiselevel
    if (test):
        print("Noise level:{:f}, auto offset:{:f}".format(noiselevel,auto_offset_value))
        if (fit_order is not None):
            plt.plot(wavelength_smooth,fitdata + noiselevel,linestyle='dashed')
        else:
            plt.plot(plt.xlim(),[noiselevel + fitdata]*2,linestyle='dashed')
    linelist = []
    line_amp_list = []
    act_ind = 0
    while True:
        if (act_ind > len(data_proc) - 3):
            break
        ind = np.nonzero(data_proc[act_ind:] > noiselevel)[0]
        if (len(ind) == 0):
            break
        ind_diff = np.diff(ind)
        ind_lineend = np.nonzero(ind_diff > 1)[0]
        if (len(ind_lineend) == 0):
            ind_lineend = len(ind)
        else:
            ind_lineend = ind_lineend[0] - 1
        if ((ind_lineend < linewidth / w_res * 10) and (ind_lineend > linewidth / w_res / 2)):
            line_data = data_proc[act_ind + ind[0] : act_ind + ind[0] + ind_lineend]
            line_w = wavelength_smooth[act_ind + ind[0] : act_ind + ind[0] + ind_lineend]
            linelist.append(np.sum(line_w * line_data) / np.sum(line_data))
            line_amp_list.append(np.max(line_data))
        act_ind += ind[0] + ind_lineend + 2
    if (test):
        if (len(linelist) != 0):
            for w in linelist:
                plt.plot([w] * 2, plt.ylim(),color='black')
    if (wavelength[-1] < wavelength[0]):
       linelist.reverse()
       line_amp_list.reverse()
    return linelist, line_amp_list,auto_offset_value,noiselevel
            
plt.close('all')   
# cxrs_summary('20230316.072',min_line_amp=0.1,integration_time=6)
#cxrs_summary('20250401.026',min_line_amp=20,integration_time=2)
#cxrs_summary('20250403.018',integration_time=1,linewidth=0.1)
#cxrs_summary('20240926.028',integration_time=1,linewidth=0.1,test_spectrum=True)
#cxrs_summary('20250403.018',integration_time=1,linewidth=0.05,test_spectrum=True)