# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 10:17:19 2025

ToDO:
    Check beam chopping time against spectoroscopy camera timestep. Avoid fast chopping measurements.

@author: zoletnik
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import numpy as np
import math

import flap
import flap_w7x_abes

flap_w7x_abes.register() 

def cxrs_summary(exp_id,linewidth=0.05,test=True,integration_time=1,
                 list_result=True,min_spectral_noiselevel=5,min_modulation=5,min_mod_SNR=3):
    """
    Finds lines in the ABES CXRS spectrum integrated for all fibres and times. Calculates the time evolulution of 
    the line intensities in each fibre and line. Selects the lines, fibres and time intervals where the intensity
    is modulated in correlation with beam modulation. 
    
    Works only with slow camera sync chopping.
    
    The algorithm is the following.
    - We add all spectra in all optical channels and find spectral lines in it.
    - In each fibre we determine the intensity of each line as a function of time and subtract an offset before the plasma.
    
    This is the old algorithm, not the actual:

    - We calculate the autocovariance function up to <tau_max> time lag of the line intensity in <integration_time> long time interval sliding 
      through the whole measeurement time. With this we assume the alkali beam modulation and the spectroscopy camera exposure are synchronized
      and every second spectrum contains CX signal from the beam.
    - In an ideal case the autocovariance function of a fully modulated line is P,0, P, 0, ..., where P is the modulation amplitude squared.
      A constant passive line intensity would add an offset to this. A modulated passive line would add some correlation function falling as a 
      function of time lag. 
    - We take the beam-times of the signal and interpolate for all timepoints. We use this to approximate the autovariance function of the 
      passive line intensity.
    - We subtract the passive autocovariance function from the full autocovariance function so as to determine the autocovariance function 
      of the active line intensity.
    - In an ideal case P/2, -P/2, P/2,-P/2 .. remains, where P is the power of the active line intensity modulation
    - We multiply the remaining autocovariance function with 1,-1,1,-1,... and take the mean times 2 to get P.
    - The relative line modulation amplitude is sqrt(P)/mean(signal)
    - We assess the error of P by calculating the variance of P * [1,-1,1,-1,..]. From the error of P we asses the error of the relative modulation.

    Parameters
    ----------
    exp_id : str
        Experiment ID.
    linewidth : float, optional
        The smooth lengt in nm for finding lines. 
        Twice this value will be used for integrating the line intensity. The default is 0.05 nm.
    min_spectral_noiselevel : float, optional
        The minimum noiselevel for the line finding algorithm. [digit]
    test : bool, optional
        Make test plot if True. The default is False.
    list_result : bool, optional
        List the result.
    integration_time : float, optional
        The integration time in the discharge for processing data. The default is 1 [sec]. 
        It will be in any case a multiple of 4 of the spectrometer time resolution.
    min_modulation : float, optional
        The minimum relative modulation amplitude [%] to list the line, time interval and channel as modulated.
        The default is 5.
    min_mid_SNR : float
        The minimum Signal-to-Noise ratio to list the line, time interval and channel as modulated.
        
        
    Returns
    -------
    active_line_list: list
        A list of dictiooaries, each describing one case when the correlation is good enough.

    """
       
    d = flap.get_data('W7X_ABES_CXRS', exp_id=exp_id,name="QSI_CXRS")
    # Summing for channels and fibres to get the lines
    d_sum = d.slice_data(summing={'Time':'Mean','Channel':'Mean'})
    wavelength = d_sum.coordinate('Wavelength')[0]
    wres = abs(wavelength[1] - wavelength[0])
    if (test):
        fig = plt.figure(figsize=(20,20))      
        gs = GridSpec(2, 3, figure=fig)
        spectrum_ax = fig.add_subplot(gs[0,0])
        line_ax = fig.add_subplot(gs[1,0])
        signal_ax = fig.add_subplot(gs[0,1:])
        mod_ax = fig.add_subplot(gs[1,1:])
        plt.sca(spectrum_ax)

    w_line_global,a_line_global,offset,noiselevel = find_lines(wavelength,d_sum.data,test=test,new_figure=False,linewidth=linewidth,auto_offset=True,min_noiselevel=min_spectral_noiselevel)
    w_interval_start = [w - linewidth for w in w_line_global]
    w_interval_stop = [w + linewidth for w in w_line_global]
    t = d.coordinate('Time',options={'Change only':True})[0].flatten()
    tres = t[1] - t[0]
    if (tres < 0.007):
        raise('Too short chopper time, spectrometer could not integrate.')
    optical_channels = d.coordinate('Optical channel',options={'Change only':True})[0].flatten()
    channels = d.coordinate('Channel',options={'Change only':True})[0].flatten()
    R = d.coordinate('Device R',options={'Change only':True})[0].flatten()
    X = d.coordinate('Device x',options={'Change only':True})[0].flatten()
    Y = d.coordinate('Device y',options={'Change only':True})[0].flatten()
    line_data = np.zeros((len(t),len(w_line_global),len(optical_channels)))
    active_line_list = []
    kernel_length = int(round(integration_time / (t[1] - t[0]))) // 4 * 4
    print('Kernel length:{:d}'.format(kernel_length))
    kernel = np.ones(kernel_length) / kernel_length

    for line_i in range(len(w_line_global)):
        for i_ch in range(len(optical_channels)): 
            print('Wavelength: {:5.1f}nm, channel:{:s}'.format(w_line_global[line_i],optical_channels[i_ch]),flush=True)
            if (test):
                line_ax.cla()
                plt.sca(line_ax)
                # The wavelength scale is not exact on this plot. This is a flap.plot problem.
#                d.plot(slicing={'Optical channel':optical_channels[i_ch]},axes=['Wavelength','Time'],plot_type='image')
                d.plot(slicing={'Optical channel':optical_channels[i_ch]},summing={'Time':'Mean'},axes=['Wavelength'])
                plt.plot([w_interval_start[line_i]]*2,plt.ylim())
                plt.plot([w_interval_stop[line_i]]*2,plt.ylim())
                plt.xlim(w_line_global[line_i]-1,w_line_global[line_i]+1)
                plt.title("Channel: {:d}".format(channels[i_ch]))
            line_data[:,line_i,i_ch] = d.slice_data(slicing={'Optical channel':optical_channels[i_ch],
                                                             'Wavelength':flap.Intervals(w_interval_start[line_i],w_interval_stop[line_i])
                                                             },
                                                    summing={'Wavelength':'Mean'}
                                                    ).data
            # Offset correction
            line_data[:,line_i,i_ch] -= line_data[0,line_i,i_ch]

            line_data_smooth = np.convolve(line_data[:,line_i,i_ch],kernel,mode='valid')
            mask = ((np.arange(len(line_data[:,line_i,i_ch])) + 1) % 2) * 2 - 1
            line_data_mod = np.convolve((line_data[:,line_i,i_ch] * mask),kernel,mode='valid') * 2
            mask1 = np.sign((np.arange(len(line_data[:,line_i,i_ch])) % 4) - 1.8)
            # This is the error estimate
            line_data_mod1 = np.abs(np.convolve((line_data[:,line_i,i_ch] * mask1),kernel,mode='valid') * 2)
            t_smooth = np.convolve(t,kernel,mode='valid')
            if (test):
                signal_ax.cla()
                plt.sca(signal_ax)
                plt.plot(t,line_data[:,line_i,i_ch])
                plt.plot(t_smooth,line_data_smooth)
                plt.legend(['Data','Smooth data'])
                plt.title('Wavelength: {:5.1f}nm, channel:{:s}'.format(w_line_global[line_i],optical_channels[i_ch]))
                mod_ax.cla()
                plt.sca(mod_ax)
                plt.plot(t_smooth,np.clip(line_data_mod / line_data_smooth * 100,-100,100))
                plt.plot(t_smooth,np.clip(line_data_mod1 / line_data_smooth * 100,-100,100))
                plt.plot(plt.xlim(),[0,0],color='black',linestyle='dotted')
                plt.ylim(-105,105)
                plt.legend(['Modulation','Modulation error'])
                plt.show()
                plt.pause(0.5)
                
         
#             mean_smooth_0 = np.convolve(line_data[:-tau_max,line_i,i_ch],kernel,mode='valid')
# #            norm_0 = np.sqrt(np.convolve((line_data[kernel_length // 2 : -(kernel_length // 2) - tau_max,line_i,i_ch] - mean_smooth_0) ** 2,kernel,mode='valid'))
#             t_corr = np.convolve(np.convolve(t[:- tau_max],kernel,mode='valid'),kernel,mode='valid')
#             # This will contain the correlation as a function of time
#             corr = np.zeros((tau_max + 1, len(t_corr)),dtype='float')
#             for i_tau in range(0,tau_max + 1):
#                 mean_smooth = np.convolve(line_data[i_tau : line_data.shape[0] - (tau_max - i_tau),line_i,i_ch],kernel,mode='valid')
# #                norm = np.sqrt(np.convolve((line_data[kernel_length // 2 + i_tau : -(kernel_length // 2) - tau_max + i_tau,line_i,i_ch] - mean_smooth) ** 2,kernel,mode='valid'))
#                 c = np.convolve((line_data[kernel_length // 2 : -(kernel_length // 2) - tau_max,line_i,i_ch] - mean_smooth_0) \
#                                   * (line_data[kernel_length // 2 + i_tau: -(kernel_length // 2) - tau_max + i_tau,line_i,i_ch] - mean_smooth),kernel,mode='valid')
#                 corr[i_tau,:] = c #/ (norm_0 * norm)
#             if (False):
#                 # Fitting a curve to the autocovariance function to determine background plasma correlations
#                 tau = np.arange(tau_max + 1)
#                 order = 1
#                 if (tau_max > 5):
#                     order = 2
#                 corr_fit_poly = np.polynomial.polynomial.polyfit(tau, corr, deg=order, full = False)
#                 fitcurves = np.zeros(corr.shape,dtype='float')
#                 x = np.broadcast_to(np.arange(tau_max + 1),corr.shape[::-1]).transpose()
#                 for i in range(order + 1):
#                     fitcurves += corr_fit_poly[i,:] * x ** i 
#                 CXRS_corr = corr - fitcurves
#             if (True):
#                 # Interpolating beam-off time to whole time and calculating autocovariance of passive light
#                 beam_off_data = line_data[1::2,line_i,i_ch]
#                 beam_off_data_interpol = np.interp(t, t[1::2], beam_off_data)
#                 mean_smooth_0_passive = np.convolve(beam_off_data_interpol[:-tau_max],kernel,mode='valid')
#                 corr_passive = np.zeros((tau_max + 1, len(t_corr)),dtype='float')
#                 for i_tau in range(0,tau_max + 1):
#                     mean_smooth_passive = np.convolve(beam_off_data_interpol[i_tau : beam_off_data_interpol.shape[0] - (tau_max - i_tau)],kernel,mode='valid')
#                     c = np.convolve((beam_off_data_interpol[kernel_length // 2 : -(kernel_length // 2) - tau_max] - mean_smooth_0_passive) \
#                                       * (beam_off_data_interpol[kernel_length // 2 + i_tau: -(kernel_length // 2) - tau_max + i_tau] - mean_smooth_passive),kernel,mode='valid')
#                     corr_passive[i_tau,:] = c
#                 CXRS_corr = corr - corr_passive 
#                 fitcurves = corr_passive
            
#             mask = np.broadcast_to(1 - (np.arange(tau_max + 1,dtype=int) % 2 ) * 2,corr.shape[::-1]).transpose()
#             P = np.clip(np.mean(CXRS_corr * mask,axis=0),0,None)
#             P_error = np.sqrt(np.mean(((CXRS_corr * mask) - P) ** 2,axis=0))
#             CXRS_rel_mod = np.sqrt(P) / mean_smooth_0[kernel_length // 2 : -(kernel_length // 2)]
#             CXRS_rel_mod_error = np.clip(np.sqrt(P) - np.sqrt(P - P_error) / mean_smooth_0[kernel_length // 2 : -(kernel_length // 2)],0,None)
# #            CXRS_rel_mod_error = np.clip(np.sqrt(P_error) / mean_smooth_0[kernel_length // 2 : -(kernel_length // 2)],0,None)
# #            if (estimate_error):
# #                error_kernel_length = int(round(integration_time / (t[1] - t[0]))) // 4 * 2 + 1
#             #     error_kernel_length = kernel_length
#             #     if (error_kernel_length < 3):
#             #         raise ValueError('Integration time is too short for CXRS modulation error estimation.')
#             #     error_kernel = np.ones(error_kernel_length) / error_kernel_length
#             #     CXRS_rel_smooth = np.convolve(CXRS_rel_mod,error_kernel,mode='valid')
#             #     extend_len = (len(CXRS_rel_mod) - len(CXRS_rel_smooth))
#             #     CXRS_rel_smooth = np.concatenate((np.full(extend_len // 2,CXRS_rel_smooth[0]),CXRS_rel_smooth,np.full(extend_len - extend_len // 2,CXRS_rel_smooth[-1])))
#             #     CXRS_rel_error = np.sqrt(np.convolve((CXRS_rel_mod - CXRS_rel_smooth) ** 2,error_kernel,mode='valid') / error_kernel_length)
#             #     CXRS_rel_error = np.concatenate((np.full(extend_len // 2,CXRS_rel_error[0]),CXRS_rel_error,np.full(extend_len - extend_len // 2,CXRS_rel_error[-1])))
#             # else:
#             #      CXRS_rel_error = None
            
#             if (test_corr):
#                 plt.figure(correlation_plot.number)    
#                 plt.clf()
#                 plt.subplot(3,2,1)
#                 absmax = max([abs(np.amax(corr)),abs(np.amin(corr))])
#                 plt.imshow(corr,aspect='auto',extent=[min(t_corr),max(t_corr),0,tau_max],origin='lower',cmap='bwr',vmin=-absmax,vmax=absmax)
#                 plt.title('Correlations')
#                 plt.xlabel('Time [s]')
#                 plt.ylabel('tau')
#                 plt.subplot(3,2,2)
#                 plt.imshow(fitcurves,aspect='auto',extent=[min(t_corr),max(t_corr),0,tau_max],origin='lower',cmap='bwr',vmin=-absmax,vmax=absmax)
#                 plt.title('Fitted correlations')
#                 plt.xlabel('Time [s]')
#                 plt.ylabel('tau')
#                 plt.subplot(3,2,3)
#                 plt.imshow(CXRS_corr,aspect='auto',extent=[min(t_corr),max(t_corr),0,tau_max],origin='lower',cmap='bwr',vmin=-absmax,vmax=absmax)
#                 plt.title('Correlation with beam modulation')
#                 plt.xlabel('Time [s]')
#                 plt.ylabel('tau')
#                 plt.subplot(3,2,4)
#                 if (CXRS_rel_mod_error is not None):
#                     plt.errorbar(t_corr,CXRS_rel_mod * 100, yerr=CXRS_rel_mod_error * 100)
# #                    plt.plot(t_corr,CXRS_rel_smooth * 100)
#                 else:
#                     plt.plot(t_corr,CXRS_rel_mod * 100)
#                 plt.title('Relative CXRS modulation')
#                 plt.xlabel('Time [s')
#                 plt.ylabel('[%]')
#                 plt.ylim(0,105)
#                 plt.suptitle("W: {:5.1f}, ch: {:s}".format(w_line_global[line_i],optical_channels[i_ch]))
#                 plt.subplot(3,2,5)
#                 plt.plot(t,line_data[:,line_i,i_ch])
#                 plt.xlabel('Time [s')
#                 plt.title('Line intensity')
#                 plt.tight_layout()
#                 plt.show()
#                 plt.pause(1)
#                 pass
           
            
            
            # if (len(ind) != 0):  
            #     start_ind = ind[0]
            #     while start_ind < ind[-1]:
            #         diff_ind = np.diff(ind[start_ind:])
            #         ind_diff = np.nonzero(diff_ind > 1)[0]
            #         if (len(ind_diff) == 0):
            #             stop_ind = ind[-1]
            #         else:
            #             stop_ind = ind[ind_diff[0] - 1]
            #         active_line_list.append({'w':w_line_global[i],
            #                                  'och':optical_channels[i_ch],
            #                                  'ch' : channels[i_ch],
            #                                  'trange':[t_corr[start_ind],t_corr[stop_ind]],
            #                                  'c_1':c_1,
            #                                  'c_2':c_2,
            #                                  'c_time': t_corr,
            #                                  'amp': line_data[:,i,i_ch],
            #                                  't': t,
            #                                  'R': R[i_ch],
            #                                  'x': X[i_ch],
            #                                  'y': Y[i_ch]
            #                                  }
            #                                 )  
            #         start_ind = stop_ind + 1
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
    
def find_lines(wavelength,data,test=True,new_figure=True,linewidth=0.05,fit_order=None,auto_offset=True,min_noiselevel=5):
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
    new_figure : bool, optional
        Prepare new figure for the test plot if True.
    linewidth : float, optional
        Expected FWHM line width in nm. The default is 0.4.
    fit order : None or int
        If not None will fit a polinomial with fit_order to the whole spectrum and subtract.
    auto_offset : bool
        If True determine the offset automatically from the signal statistics.
    min_noiselevel : float
        The minimum noiselevel [digit] to use. If the automatically detected noiselevel is below this
        this value will be used.

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
        if (new_figure):
            plt.figure()
        else:
            plt.cla()
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
#cxrs_summary('20250401.026',integration_time=1,linewidth=0.1,test_spectrum=True,test_corr=True,tau_max=4)   # 529 nm
#cxrs_summary('20250403.018',integration_time=1,linewidth=0.1)  # 530 nm
#cxrs_summary('20240926.028',integration_time=1,linewidth=0.1,test_spectrum=True,test_corr=True,tau_max=4)  # 529 nm
cxrs_summary('20250402.028',integration_time=2,linewidth=0.02,test=True,min_spectral_noiselevel=150)   # 584 nm  Na line
#cxrs_summary('20250402.064',integration_time=1,linewidth=0.05,test_spectrum=True,test_corr=True,tau_max=4)   # 585 nm
