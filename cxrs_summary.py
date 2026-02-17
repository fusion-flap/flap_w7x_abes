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
from matplotlib.patches import Rectangle
from matplotlib import cm
import pickle
import numpy as np
import math

import flap
import flap_w7x_abes

flap_w7x_abes.register() 

def cxrs_summary(exp_id,linewidth=0.05,test=True,integration_time=1,
                 list_result=True,min_spectral_noiselevel=5,min_modulation=5,min_mod_SNR=3,savefile='ABES_CXRS_save.dat',
                 show=True):
    """
    Finds lines in the ABES CXRS spectrum integrated for all fibres and times. Calculates the time evolulution of 
    the line intensities in each fibre and line. Selects the lines, fibres and time intervals where the intensity
    is modulated in correlation with beam modulation. 
    
    Works only with slow camera sync chopping when consequitive spectra are beam-on, beam-off.
    
    The algorithm is the following.
    - We add all spectra in all optical channels and find spectral lines in it.
    - In each fibre we determine the intensity of each line as a function of time and subtract an offset before the plasma.
    - We calculate the sum of intensity after multiplying with 1, -1, 1, -1,... This is the amplitude of the change at the modulation frequency, 
      it gives the modulation amplitude.
    - We calculate the sum of intensity after multiplying with 1, 1, -1, -1, 1, 1, -1, -1... This is the amplitude of the change at 2 times the 
      modulation frequency. This is considered as the error of modulation.
    

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
    min_mod_SNR : float
        The minimum Signal-to-Noise ratio to list the line, time interval and channel as modulated.
    savefile : str
        The file where data will be saved.
    show : Plot the result using show_modulation()
        
        
    Returns
    -------
    active_line_list: list
        A list of dictionaries, each describing one case when the correlation is good enough.

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

    res = {}
    for line_i in range(len(w_line_global)):
        res_line = {}
        for i_ch in range(len(optical_channels)): 
            if (optical_channels[i_ch] == 'N.A.'):
                continue
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
            line_data_mod_error = np.abs(np.convolve((line_data[:,line_i,i_ch] * mask1),kernel,mode='valid') * 2)
            t_smooth = np.convolve(t,kernel,mode='valid')
            smooth_tres = tres * len(kernel)
            if (test):
                signal_ax.cla()
                plt.sca(signal_ax)
                plt.plot(t,line_data[:,line_i,i_ch])
                plt.plot(t_smooth,line_data_smooth)
                plt.legend(['Data','Smooth data'])
                plt.title('Wavelength: {:5.1f}nm, channel:{:s}'.format(w_line_global[line_i],optical_channels[i_ch]))
                mod_ax.cla()
                plt.sca(mod_ax)
                plt.plot(t_smooth,np.clip(line_data_mod / line_data_smooth  * 100,-100,100))
                plt.plot(t_smooth,np.clip(line_data_mod_error / line_data_smooth * 100,-100,100))
                plt.plot(plt.xlim(),[0,0],color='black',linestyle='dotted')
                plt.ylim(-105,105)
                plt.legend(['Modulation','Modulation error'])
                plt.show()
                plt.pause(0.5)
                            
            res_line[optical_channels[i_ch]] = {'w':w_line_global[line_i],
                                                'och':optical_channels[i_ch],
                                                'ch' : channels[i_ch],
                                                'time':t_smooth,
                                                'tres':smooth_tres,
                                                'amp': line_data_smooth,
                                                'line_data' : line_data_smooth,
                                                'line_mod' : line_data_mod,
                                                'line_mod_error' : line_data_mod_error,                                             
                                                't': t_smooth,
                                                'R': R[i_ch],
                                                'x': X[i_ch],
                                                'y': Y[i_ch]
                                                }
                                                
            
            ind = np.nonzero(np.logical_and(line_data_mod > line_data_mod_error * min_mod_SNR,
                                            line_data_mod > min_modulation
                                            )
                             )[0]
            if (len(ind) != 0):  
                ind_ind = 0
                while ind_ind < len(ind):
                    start_ind = ind[ind_ind]
                    diff_ind = np.diff(ind[ind_ind:])
                    ind_diff = np.nonzero(diff_ind > 1)[0]
                    if (len(ind_diff) == 0):
                        stop_ind = ind[-1]
                        ind_ind = len(ind)
                    else:
                        stop_ind = ind[ind_ind:][ind_diff[0]]
                        ind_ind += ind_diff[0] + 1
                    active_line_list.append({'w':w_line_global[line_i],
                                             'och':optical_channels[i_ch],
                                             'ch' : channels[i_ch],
                                             'trange':[t_smooth[start_ind],t_smooth[stop_ind]],
                                             'amp': line_data_smooth,
                                             't': t_smooth,
                                             'R': R[i_ch],
                                             'x': X[i_ch],
                                             'y': Y[i_ch]
                                             }
                                            )  
            res[w_line_global[line_i]] = res_line
            
    
    if (savefile is not None):
        with open(savefile,"wb") as f:
            pickle.dump([exp_id,res],f)
    
    if (show):
        show_modulation(res=res,title=exp_id)
            
    if (list_result):
        for l in active_line_list:
            print("Wavelength: {:5.1f}nm, Optical Ch: {:s}, channel: {:d}, R:{:5.3f}, time range: [{:4.1f},{:4.1f}]".format(l['w'],l['och'],l['ch'],l['R'],*l['trange']))                
    return active_line_list

def show_modulation(res=None,file=None,min_mod_SNR=3,min_modulation=5,R_range=[6.1,6.3],
                    trange=None,mod_range=[0,100],title=None):
    """
    Plot the modulation of the lines as a function of time and R. Will sum up the signals 
    in all fibres located at the same R.
    Only plots data where the modulation of the line is at least <min_modulation> and
    the Signal to Noise Ratio of the modulation is at least <min_mod_SNR>.

    Parameters
    ----------
    res : dict, optional
        The result returned by cxrs_summary() . The default is None.
    file :str, optional
        File name with the data written by cxrs_summary(). The default is None.
        If res is None and file is not None data will be read from the file.
    min_mod_SNR : float, optional
        The minimum Signal to Noise Ratio of the modulation to plot. The default is 3.
    min_modulation : float, optional
        The minimum modulation for plotting in %. The default is 5.
    R_range : list of two floats, optional
        The major radius range of the plot. The default is [6.1,6.3].
    trange : list of two floats, optional
        The time range of the plot in second. The default is None, it will plot all data.
    mod_range : lst of two floats, optional
        The modulation range for the plot in %. The default is [0,100].
    title :str, optional
        The title of the plot. The default is None, in this case the exp_id and the plot parameters will be listed.modulation

    Returns
    -------
    None.

    """
    if ((res is None) and (file is not None)):
        with open(file,"rb") as f:
            exp_id,res = pickle.load(f)
            if (title is None):
                _title = 'Line modulations. min_modulation={:f}%, min_mod_SNR={:f}   experiment: {:s}.'.format(min_modulation,min_mod_SNR,exp_id)
    else:
        _title = title
            
    good_wl = []
    if (trange is not None):
        _trange = trange
    else:
        _trange = [1e4,-1]
    for w in res.keys():
        for f in res[w].keys():
            ind = np.nonzero(np.logical_and(res[w][f]['line_mod'] > res[w][f]['line_mod_error'] * min_mod_SNR,
                                            res[w][f]['line_mod'] / res[w][f]['line_data'] * 100 > min_modulation
                                            )
                                            )[0]
            if (len(ind) > 0):
                good_wl.append(w)
                _trange[0] = min([_trange[0],min(res[w][f]['time'])])
                _trange[1] = max([_trange[1],max(res[w][f]['time'])])
                break        
    if (len(good_wl) == 0):
        print("No modulation was found.")
    n_wl = len(good_wl)
    fig = plt.figure(figsize=(20,15))
    
    if (n_wl <= 4):
        nc = n_wl
        nr = 1
    else:
        nc = int(round(math.sqrt(n_wl)))
        nr = n_wl // nc
        if (n_wl % nc != 0):
            nr += 1
    fibre_image_size = 2.5 # mm
    for i,w in enumerate(good_wl):    
        # This will collect the intensity and modulation as a function of R
        mod_amp = []
        amp = []
        err_amp = []
        plot_R = []
        for f in res[w].keys():
            append = True
            if (len(plot_R) > 0):
                if (np.min(np.abs(np.array(plot_R) - res[w][f]['R'])) < fibre_image_size * 1e-3 / 2):
                    ind = np.argmin(np.abs(np.array(plot_R) - res[w][f]['R']))
                    mod_amp[ind] += res[w][f]['line_mod']
                    amp[ind] += res[w][f]['line_data']
                    err_amp[ind] += res[w][f]['line_mod_error']
                    append = False
            if (append):
                plot_R.append(res[w][f]['R'])
                mod_amp.append(res[w][f]['line_mod'])
                amp.append(res[w][f]['line_data'])
                err_amp.append(res[w][f]['line_mod_error'])
                    
        ax = plt.subplot(nr,nc,i+1)
        plt.xlim(*_trange)
        plt.xlabel('Time [s]')
        plt.ylim(*R_range)
        plt.ylabel('R [cm]')
        plt.title("{:5.1f}nm".format(w))
        time = res[w][next(iter(res[w]))]['time']
        tres = res[w][next(iter(res[w]))]['tres']
        for iR,R in enumerate(plot_R):
            rel_mod = mod_amp[iR] / amp[iR] * 100
            ind = np.nonzero(np.logical_and(mod_amp[iR] > err_amp[iR] * min_mod_SNR,
                                            rel_mod > min_modulation
                                            )
                             )[0]
            c_ind = np.clip((rel_mod - mod_range[0]) /( mod_range[1] - mod_range[0]),0,1)
            for ii in ind:
                patch = Rectangle((time[ii] - tres / 2,R - fibre_image_size * 1e-3 / 2),
                                                     width = tres,
                                                     height = fibre_image_size * 1e-3,
                                                     color = (1,1-float(c_ind[ii]),1-float(c_ind[ii]))
                                                     )
                ax.add_patch(patch)
             
    
    if (_title is not None):
        plt.suptitle(_title)
    
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.85, hspace=0.1)
    #plt.tight_layout()    
    
    # Adding a custom colorscale
    ax_color = fig.add_subplot([0.1,0.9,0.84,0.02])
    plt.xlim(*mod_range)
    plt.xlabel('Modulation [%]')
    plt.ylim(0,1)
    plt.yticks([])
    plt.title('Modulation scale')
    nc = 100
    for i in range(nc):
        patch = Rectangle((i/nc*100, 0), width=100/nc , height = 1, color=(1,1-i/nc,1-i/nc))
        ax_color.add_patch(patch)
    
    plt.show()

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
#cxrs_summary('20250401.026',integration_time=1,linewidth=0.2,test=True)   # 529 nm
cxrs_summary('20250403.018',integration_time=2,linewidth=0.2,test=True)  # 530 nm
#cxrs_summary('20240926.028',integration_time=1,linewidth=0.1,test_spectrum=True,test_corr=True,tau_max=4)  # 529 nm
#cxrs_summary('20250402.028',integration_time=4,linewidth=0.2,test=False,min_spectral_noiselevel=150)   # 584 nm  Na line
#cxrs_summary('20250402.064',integration_time=4,linewidth=0.2,test=True,min_spectral_noiselevel=30)   # 585 nm

show_modulation(file='ABES_CXRS_save.dat')