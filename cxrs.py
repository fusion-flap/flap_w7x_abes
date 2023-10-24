# -*- coding: utf-8 -*-
"""
Created on Tue Aug 8 15:27:49 2023

@author: bcsillag

Data processing code for Wendelstein 7-X QSI CXRS spectra measured during OP2.1
"""

import requests
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy import interpolate
from scipy import sparse as sp
import flap
import flap_w7x_abes
    
def parab(x,a,b,c):
    return a*x**2 + b*x + c

def gauss(lambd,s,A,up,down):
    centr = (up + down) / 2
    return A*np.e**(-(((lambd-centr)**2)/s**2))

def doppler_broadening(kbt,a): #pontosítandó
    M = 19.9419 
    return a*np.sqrt(2*1.602176487*abs(kbt)/M) / 30000

def wavelength_grid_generator(grid,wavelength_setting,roi):
    calib_array = np.loadtxt("wavelength_calib_2023_"+grid+".txt") #loading the calibration coeffs
    c0 = calib_array[roi-1,0]
    c1 = calib_array[roi-1,2]
    c2 = calib_array[roi-1,4]
    c3 = calib_array[roi-1,6] #writing the out to clear variables
    pix_values = np.arange(0,1024,1) - 511 
    return c3*pix_values**3 + c2*pix_values**2 + c1*pix_values + c0 + wavelength_setting

def wavelength_calib(dataobj,grid,roi,wl):
    """
    It performs wavelength calibration on the loaded raw spectra.
    dataobj: data object containing the measured spectra
    grid: applied grid in the Isoplane spectrometer in the corresponding discharge.
            It can be "1200g_per_mm", "1800g_per_mm" or "2400g_per_mm"
    roi: region of interest - in other words, channel. It can be P01-P06
    wl: central wavelength setting during the corresponding discharge
    """
    wl_values = wavelength_grid_generator(grid,wl,roi) #calibration
    
    wl_coord_flap = flap.Coordinate(name='Wavelength', unit="nm", 
                                    mode=flap.CoordinateMode(equidistant=False), shape=(1024),
                                  values=wl_values, dimension_list=[2])
    dataobj.add_coordinate_object(wl_coord_flap)
    dataobj.del_coordinate("Coord 2") #exchanging coordinates: pixel -> wavelength
    return dataobj

def get_spectra(expe_id,grid,wavelength_set,roi):
    """
    It can load the QSI CXRS spectra that were measured during OP2.1.
    expe_id: shot ID
    grid: applied grid in the Isoplane spectrometer in the corresponding discharge.
            It can be "1200g_per_mm", "1800g_per_mm" or "2400g_per_mm"
    wavelength_set: central wavelength setting during the corresponding discharge
    roi: region of interest - in other words, channel. It can be P01-P06
    """
    #get data by shotID
    qsi_cxrs = flap.get_data('W7X_WEBAPI', name='Test/raw/W7X/QSI/cxrs_DATASTREAM/0/Images/',
                              exp_id=expe_id,
                              options={'Scale Time': True,
                                      'Cache Data': False},
                              object_name='QSI_CXRS_data')

    #changing to ROI coord
    roi_coord_flap = flap.Coordinate(name='ROI', unit="1", mode=flap.CoordinateMode(equidistant=False), shape=(6),
                                  values=["P01","P02","P03","P04","P05","P06"], dimension_list=[1])
    qsi_cxrs.add_coordinate_object(roi_coord_flap)
    qsi_cxrs.del_coordinate("Coord 1")
    
    #adding coordinate Device R
    R_coord_flap = flap.Coordinate(name='Device R', unit="m", mode=flap.CoordinateMode(equidistant=False), shape=(6),
                                  values=[6.175608381963992,6.18531258540187,6.196755395927777,
                                          6.207903188440903,6.22187706429289,6.241915857714211], dimension_list=[1])
    qsi_cxrs.add_coordinate_object(R_coord_flap)

    #correcting time coordinate
    c=qsi_cxrs.get_coordinate_object("Time")
    c.values = c.values + 1
    c=qsi_cxrs.get_coordinate_object("Time")

    #wavelength calibration
    qsi_cxrs = wavelength_calib(qsi_cxrs,grid,roi,wavelength_set)
    
    return qsi_cxrs

def slice_by_wl(spectra,w,expe_id,roi):
    """
    It can get and plot the temporal evolution of 1 pixel (closest to the given wavelength).
    spectra: dataobject containing the wavelength calibrated spectra
    w: the wavelength along wich the dataobject is sliced
    expe_id: shot ID
    roi: region of interest - in other words, channel. It can be P01-P06
    """
    spectra_1w = spectra.slice_data(slicing={"ROI" :"P0"+str(roi)})
    spectra_1w = spectra_1w.slice_data(slicing={"Wavelength" : w})
    plt.figure()
    spectra_1w.plot(axes="Time")
    plt.xlabel("Time [s]",fontsize = 15)
    plt.ylabel("Intensity",fontsize = 15)
    plt.title("Time evolution of a pixel at wavelength "+str(w)+" nm")
    
def interval_shift(expe_id):
    shift = 0
    if(expe_id[:8]=="20230314"):
        shift = -0.0509
    elif(expe_id[:8]=="20230315"):
        shift = -0.05087939698492462
    elif(expe_id=="20230316.043"):
        shift = -0.037
    elif(expe_id=="20230316.072"):
        shift = -0.05012562814070352
    elif(expe_id[:8]=="20230323"):
        shift = 0.04764529058116232
        # shift = -0.07229458917835671#-0.0592
    elif(expe_id[:8]=="20230328"):
        shift = 0.051633
    elif(expe_id[:8]=="20230330"):
        shift = 0.05275551102204408#0.051633
    
    return shift
    
def slice_by_wl_range(spectra,wstart,wstop,expe_id,roi,timerange):
    spectra_1w = spectra.slice_data(slicing={"ROI" :"P0"+str(roi)})
    spectra_1w = spectra_1w.slice_data(slicing={"Wavelength" : flap.Intervals(wstart, wstop)},summing = {"Wavelength":"Mean"})
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
        
    #correcting the timescales
    el = interval_shift(expe_id)
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
        
    # plt.figure()
    # legend = []
    # legend.append('Beam on')
    # legend.append('Beam off')
    
    spectra_1w.plot(axes="Time")
    # d_beam_off.plot(axes=['Time',650],plot_type='scatter',options={'Force':True})
    # d_beam_on.plot(axes=['Time',660],plot_type='scatter',options={'Force':True})
    plt.xlabel("Time [s]",fontsize = 15)
    plt.ylabel("Intensity",fontsize = 15)
    plt.title("Time evolution of a pixel at wavelength ["+str(wstart)+", "+str(wstop)+"] nm")
    # legend.append('avg. spectral intensity')
    
    # c = spectra_1w.get_coordinate_object("Time")
    # timeres = np.mean(c.values[1:(len(c.values)-1)]-c.values[0:(len(c.values)-2)])
    # t0 = c.values[0]
    # c.mode.equidistant = True
    # c.shape = spectra_1w.data.shape
    # c.start = t0
    # c.step = timeres
    # c.dimension_list = [0]
    
    # ap = spectra_1w.apsd(options={'Interval':1, 'Range':[1,100],'Res':0.1}) 
    # plt.figure()
    # ap.plot(axes="Frequency",options={'Log x': True, 'Log y':True, 'Error':True, 'X range':[1,100]})
    # plt.title(" ")
    
def active_passive(spectra,roi,t_start,t_stop,expe_id,timerange):
    """
    A function that calculates the average spectra at beam-on and beam-off,
    then it subtracts  the two to gain the active spectrum
    spectra: dataobject containing the wavelength calibrated spectra
    roi: region of interest - in other words, channel. It can be P01-P06
    t_start: start of the non-transient plasma interval
    t_stop: end of the non-transient plasma interval
    el: temporal shift between the CCD and the beam clocks
    expe_id: shot ID
    timerange: temporal range of beam intervals that one aims to utilize
    """
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    #correcting the timescales
    el = interval_shift(expe_id)
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    s_on_intervals_full=ROI1.slice_data(slicing={'Time':d_beam_on})
    s_off_intervals_full=ROI1.slice_data(slicing={'Time':d_beam_off})
    s_on_intervals = s_on_intervals_full.data[:,1:,].mean(axis = 1)
    s_off_intervals = s_off_intervals_full.data[:,1:,].mean(axis = 1)
    s_on_data = s_on_intervals.mean(axis = 1)
    s_off_data = s_off_intervals.mean(axis = 1)
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
    s_on.data = s_on_data
    s_off.data = s_off_data
    
    # print(s_on.data.shape)
    # raise ValueError("Stop")
    plt.figure()
    s_on.plot(axes = "Wavelength")
    s_off.plot(axes = "Wavelength")
    legend = ["beam on","beam off"]
    plt.legend(legend)
    plt.title(expe_id,fontsize = 15)
    plt.title(expe_id+", averaged spectra, ROI = P0"+str(roi),fontsize = 15)
    plt.grid()
    
    s_subs = s_on
    s_subs.data = s_on.data-s_off.data #gaining the active spectrum
    plt.figure()
    s_subs.plot(axes = "Wavelength")
    plt.title(expe_id+", active spectrum, ROI = P0"+str(roi),fontsize = 15)
    plt.grid()
    
def spectral_error_calc(spec):
    spec_perint = np.zeros((spec.data.shape[0],spec.data.shape[2]))
    # avg_int_length = 0
    for i in range(spec.data.shape[2]):
        ind = np.nonzero(abs(spec.data[10,:,i]) > 1e-10)
        # avg_int_length += len(ind[0]) / spec.data.shape[2]
        for j in ind[0]:
            spec_perint[:,i] = spec_perint[:,i] + spec.data[:,j,i] / len(ind[0])
            
    # linear = lambda x,a,b : a*x + b
    # for i in range(spec.data.shape[0]):
    #     popt,param = curve_fit(linear, np.arange(0,spec.data.shape[2],1), spec_perint[i,:])
    #     spec_perint[i,:] = spec_perint[i,:] - linear(np.arange(0,spec.data.shape[2],1), *popt)
    return np.sqrt(spec_perint.var(axis = 1)/ (spec.data.shape[2]))

def get_line_intensity(spectra,roi,t_start,t_stop,expe_id,timerange,lstart,lstop,bg_wls=[0],plots=False):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    el = interval_shift(expe_id)
    #correcting the timescales
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    ROI1 = ROI1.slice_data(slicing={"Wavelength" :flap.Intervals(lstart, lstop)})
    
    
    if(bg_wls != [0]):
        ROI1_witbg = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),
                                        "Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                        summing = {"Wavelength":"Mean","Time":"Mean"})
        ROI1.data[:,:] = ROI1.data[:,:] - ROI1_witbg.data
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
    
    lambd = s_on.coordinate("Wavelength")[0]
    gaus = lambda x,A,s,mu : A*np.e**(-(((x-mu)**2)/s**2))
    
    popton, pcovon = curve_fit(gaus,lambd, s_on.data,p0 = [max(s_on.data),0.1,lambd.mean()])
    poptoff, pcovoff = curve_fit(gaus,lambd, s_off.data,p0 = [max(s_off.data),0.1,lambd.mean()])
    if(plots == True):
        fs = 15
        plt.figure()
        plt.plot(lambd,s_on.data,"+")
        plt.plot(lambd,gaus(lambd,*popton))
        plt.xlabel("Wavelength [nm]",fontsize = fs)
        plt.ylabel("Spectral intensity",fontsize = fs)
        plt.legend(["Data","Fit"],fontsize = fs-2)
        plt.title(expe_id+", Beam on line intensity fit")
        
        plt.figure()
        plt.plot(lambd,s_off.data,"+")
        plt.plot(lambd,gaus(lambd,*poptoff))
        plt.xlabel("Wavelength [nm]",fontsize = fs)
        plt.ylabel("Spectral intensity",fontsize = fs)
        plt.legend(["Data","Fit"],fontsize = fs-2)
        plt.title(expe_id+", Beam off line intensity fit")
        
    return popton[0], poptoff[0]

def indep_spectral_error_calc(spec):
    spec_perint = np.zeros((spec.data.shape[0],spec.data.shape[2]))
    for i in range(spec.data.shape[2]):
        ind = np.nonzero(abs(spec.data[10,:,i]) > 1e-10)
        for j in ind[0]:
            spec_perint[:,i] = spec_perint[:,i] + spec.data[:,j,i] / len(ind[0])
    # for i in range(spec_perint.shape[0]):
    #       spec_perint[i,:] = spec_perint[i,:] - spec_perint[i,:].mean()
         
    HR = np.zeros((spec.data.shape[0]))
    HL = np.zeros((spec.data.shape[0]))
    for i in range(1,spec_perint.shape[0]-1):
        M1 = spec_perint[i,:]
        HR[i]=np.mean(M1**2)-np.mean(M1*spec_perint[i+1,:])*M1.mean()/spec_perint[i+1,:].mean()
        HL[i]=np.mean(M1**2)-np.mean(M1*spec_perint[i-1,:])*M1.mean()/spec_perint[i-1,:].mean()
    H = (HR + HL) / 2
    M1 = spec_perint[0,:]
    H[0] = np.mean(M1**2)-np.mean(M1*spec_perint[1,:])*M1.mean()/spec_perint[1,:].mean()
    M1 = spec_perint[-1,:]
    H[-1] = np.mean(M1**2)-np.mean(M1*spec_perint[-2,:])*M1.mean()/spec_perint[-2,:].mean()
    err = np.sqrt( abs(H) / spec_perint.shape[1])
    # print(err)
    # raise ValueError("STOP")
    # plt.figure()
    # plt.plot(err,marker = "+")
    return err
    
def active_passive_with_error(spectra,roi,t_start,t_stop,expe_id,timerange,bg_wls=[0],plots=False):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    #correcting the timescales
    el = interval_shift(expe_id)
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    s_on_sliced=ROI1.slice_data(slicing={'Time':d_beam_on}) #for error calculation
    s_off_sliced=ROI1.slice_data(slicing={'Time':d_beam_off})
    
    #the intervals are taken as independent measurements
    error_on = spectral_error_calc(s_on_sliced)
    error_off = spectral_error_calc(s_off_sliced)
    
    if(bg_wls != [0]):
        ROI1_witbg = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),
                                        "Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                        summing = {"Wavelength":"Mean","Time":"Mean"})
        ROI1.data[:,:] = ROI1.data[:,:] - ROI1_witbg.data
        
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
    s_on_intervals_full=ROI1.slice_data(slicing={'Time':d_beam_on})
    s_off_intervals_full=ROI1.slice_data(slicing={'Time':d_beam_off})
    s_on_intervals = s_on_intervals_full.data[:,1:,].mean(axis = 1)
    s_off_intervals = s_off_intervals_full.data[:,1:,].mean(axis = 1)
    s_on_data = s_on_intervals.mean(axis = 1)
    s_off_data = s_off_intervals.mean(axis = 1)
    s_on.data = s_on_data
    s_off.data = s_off_data
    
    if(plots == True):
        plt.figure()
        s_on.plot(axes = "Wavelength")
        s_off.plot(axes = "Wavelength")
        legend = ["beam on","beam off"]
        plt.legend(legend)
        plt.title(expe_id,fontsize = 15)
        plt.title(expe_id+", averaged spectra, ROI = P0"+str(roi),fontsize = 15)
        plt.grid()
    
    s_subs = s_on
    s_subs.data = s_on.data-s_off.data #gaining the active spectrum
    s_subs.error = np.sqrt(error_on**2 + error_off**2)
    
    if(plots == True):
        plt.figure()
        s_subs.plot(axes = "Wavelength")
        plt.title(expe_id+", active spectrum, ROI = P0"+str(roi),fontsize = 15)
        plt.grid()
    
    return s_subs

def error_distr(spectra,roi,t_start,t_stop,expe_id,minint,timerange,lstart,lstop,bg_wls=[0],plots=False):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    #correcting the timescales
    el = interval_shift(expe_id)
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    ROI1 = ROI1.slice_data(slicing={"Wavelength" :flap.Intervals(lstart, lstop)})
    s_on_sliced=ROI1.slice_data(slicing={'Time':d_beam_on}) #for error calculation
    s_off_sliced=ROI1.slice_data(slicing={'Time':d_beam_off})
    
    
    
    #the intervals are taken as independent measurements
    error_on = list(indep_spectral_error_calc(s_on_sliced))
    error_off = list(indep_spectral_error_calc(s_off_sliced))
    if(bg_wls != [0]):
        ROI1_witbg = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),
                                        "Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                        summing = {"Wavelength":"Mean","Time":"Mean"})
        ROI1.data[:,:] = ROI1.data[:,:] - ROI1_witbg.data
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
    int_on = list(s_on.data)
    int_off = list(s_off.data)
    
    err = np.array(error_on + error_off)
    intensity = np.array(int_on + int_off)
    
    err = err[intensity.argsort()]
    intensity = intensity[intensity.argsort()]
    
    err = err[intensity > minint]
    intensity = intensity[intensity >minint]    
    
    sq = lambda x,a,b : a*np.sqrt(x) + b
    if(plots == True):
        fs = 15
        plt.figure()
        plt.plot(intensity,err,"+")
        popt, pcov = curve_fit(sq, intensity, err,p0=[0.03,0.3])
        plt.plot(intensity,sq(intensity,*popt))
        plt.xlabel("Spectral intensity",fontsize = fs)
        plt.ylabel("Error",fontsize = fs)
        plt.legend(["Data","Fit"],fontsize = fs-2)
        plt.grid()
    return popt

def passive(qsi_cxrs,roi,t_start,t_stop,expe_id):
    
    plt.figure()
    ROI1 = qsi_cxrs.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    s_on=ROI1.slice_data(summing = {"Time":"Mean"})
    
    s_on.plot(axes = "Wavelength")
    plt.title(expe_id,fontsize = 15)
    plt.title(expe_id+", averaged spectra, ROI = P0"+str(roi),fontsize = 15)
    plt.grid()
    
def tshif(qsi_cxrs,roi,t_start,t_stop,expe_id,timerange,wstart,wstop,N,bg_wls = [0]):
    tsh = np.linspace(-0.075,0.075,N)
    lineint = np.zeros((N))
    for i in range(len(tsh)):
        d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                                 options={'State':{'Chop': 0, 'Defl': 0}},\
                                 object_name='Beam_on',
                                 coordinates={'Time': timerange})
    
        d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                                 options={'State':{'Chop': 1, 'Defl': 0}},\
                                 object_name='Beam_off',
                                 coordinates={'Time': timerange})
        
        c=d_beam_on.get_coordinate_object("Time")
        c.start = c.start + tsh[i]
        c=d_beam_off.get_coordinate_object("Time")
        c.start = c.start + tsh[i]
        
        ROI1 = qsi_cxrs.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
        if(bg_wls != [0]):
            ROI1_witbg = qsi_cxrs.slice_data(slicing={"ROI" :"P0"+str(roi),
                                            "Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                            summing = {"Wavelength":"Mean","Time":"Mean"})
            ROI1.data[:,:] = ROI1.data[:,:] - ROI1_witbg.data
            
        s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
        s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
        s_on_intervals_full=ROI1.slice_data(slicing={'Time':d_beam_on})
        s_off_intervals_full=ROI1.slice_data(slicing={'Time':d_beam_off})
        s_on_intervals = s_on_intervals_full.data[:,1:,].mean(axis = 1)
        s_off_intervals = s_off_intervals_full.data[:,1:,].mean(axis = 1)
        s_on_data = s_on_intervals.mean(axis = 1)
        s_off_data = s_off_intervals.mean(axis = 1)
        s_on.data = s_on_data
        s_off.data = s_off_data
        
        s_subs = s_on
        s_subs.data = s_on.data-s_off.data
        lineint[i] = s_subs.slice_data(slicing={"Wavelength" : flap.Intervals(wstart, wstop)},
                                       summing = {"Wavelength":"Mean"}).data
    
    fs = 15
    plt.figure()
    plt.plot(tsh,lineint,marker = "o")
    plt.xlabel("$t_{shift}$ [s]",fontsize = fs)
    plt.ylabel("Spectral intensity", fontsize = fs)
    plt.title("Na 3s->3p, $\Delta m = 0$, ["+str(t_start)+", "+str(t_stop)+"] s",fontsize = fs)
    
    print("The best temporal shift length in second:")
    print(tsh[np.argmax(lineint)])
    
def autocorr(qsi_cxrs,roi,t_start,t_stop,lstart,lstop,expe_id):
    ROI1 = qsi_cxrs.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    dt = np.mean(ROI1.coordinate("Time")[0][1:,0]-ROI1.coordinate("Time")[0][:-1,0])
    line = ROI1.slice_data(slicing={"Wavelength" :flap.Intervals(lstart, lstop)},
                           summing ={"Wavelength":"Sum"})
    flap.list_data_objects(line)
    c=line.get_coordinate_object("Time")
    c.mode.equidistant = True
    c.shape = line.data.shape
    c.step = dt
    c.start = 0.0
    c.dimension_list = [0]
    t = ROI1.coordinate("Time")[0][:,0]
    popt, pcov = curve_fit(parab, t, line.data)
    std_mea = line.data - parab(t,*popt)
    v = std_mea.var()
    std_mea = std_mea/v
    tcorr = t-t.mean()
    corr = np.correlate(std_mea, std_mea,mode = "same")
    
    fs = 15
    plt.figure()
    plt.plot(tcorr,corr,"bo-")
    plt.xlabel(r'$\tau [s]$',fontsize = fs)
    plt.title(expe_id+", Auto-correlation, $\lambda = ["+str(lstart)+", "+str(lstop)+"]$ nm")
    plt.grid()
    
    
def normed_intensity(qsi_cxrs,grid,t_start,t_stop,expe_id,timerange,lstart,lstop,bg_wls):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    el = interval_shift(expe_id)
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    qsi_cxrs = qsi_cxrs.slice_data(slicing={"Time":flap.Intervals(t_start, t_stop)})
    if(bg_wls != [0]):
        qsi_cxrs_bg = qsi_cxrs.slice_data(slicing={"Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                        summing = {"Wavelength":"Mean","Time":"Mean"})
        for i in range(qsi_cxrs_bg.data.shape[0]):
            qsi_cxrs.data[:,i,:] = qsi_cxrs.data[:,i,:] - qsi_cxrs_bg.data[i]
    s_on=qsi_cxrs.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
    s_off=qsi_cxrs.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
    s_subs = s_on
    s_subs.data = s_on.data-s_off.data
    
    s_subs = s_subs.slice_data(slicing={"Wavelength" :flap.Intervals(lstart, lstop)})
    
    gridfac = 0
    if(grid == "1200g_per_mm"):
        gridfac = 1
    elif(grid == "1800g_per_mm"):
        gridfac = 12572/12974
    elif(grid == "2400g_per_mm"):
        gridfac = 7620/12572
        
    lambd = s_on.coordinate("Wavelength")[0]
    gaus = lambda x,A,s,mu,C : C + A*np.e**(-(((x-mu)**2)/s**2))/s
    
    lineint = []
    int_calib = np.array([5961,5667,1209,3367])/5961
    for i in range(4):
        y = s_subs.data[i,:].ravel()
        lambd = s_subs.coordinate("Wavelength")[0][i,:]
        popton, pcovon = curve_fit(gaus,lambd, y,p0 = [max(y),0.1,lambd.mean(),0])
        # plt.figure()
        # plt.plot(lambd,y,"b+")
        # plt.plot(lambd,gaus(lambd,*popton))
        lineint.append(max(gaus(lambd,*popton)))
        print(popton[3])
    lineint = np.array(lineint)
    print(lineint)
    return gridfac*lineint/int_calib

def Li_spectrum(B,E,v,n,isotopemass):
    # input parameters: (they are lists, so they are not vectors or matrices)
    # B = [Bx,By,Bz]    Cartesian components of the magneric field vector (in Tesla) in the laboratory frame
    # E = [Ex,Ey,Ez]    Cartesian components of the electric field vector (in V/m) in the laboratory frame 
    # v = [vx,vy,vz]    Cartesian components of the atomic velocity (in m/s) in the laboratory frame
    # n = [nx,ny,nz]    Cartesian components of the viewing direction in the laboratory frame 
    # (direction TO detector, i.e. propagation direction of light)
        
    # output: structure spectrum with the following fields:
    # lambdas           the wavelengths of the spectral lines in nm
    # intensities       the relative intensities 
    
    # isotopemass       6 or 7, the mass of the isotope considered.
    
    B=np.array(B)
    E=np.array(E)
    v=np.array(v)
    n=np.array(n)

    # n is normalized, just in case:
    
    n = n/np.sqrt(sum(n**2))
    
    # physical constants:
    c = 299792458; # vacuum light speed in m/s
    muB = 9.27400915e-24; # Bohr magneton in J/T
    hbar = 1.05457172600000e-34;
    muB_per_h = 1e10*9.27400915/(2 * np.pi*1.054571726)
    h = hbar * 2 * np.pi;
    
    # relevant Li parameters
    gS = 2.0023193043622; # electron spin g-factor  
    gL = 0.99997613; # electron orbital g-factor 
    
    if(isotopemass == 6):
        E_D1 = 446.789598e12; # 6Li 3P1/2 excited state energy relative to ground state (E/h) in Hz
        E_D2 = 446.799650e12; # 6Li 3P3/2 excited state energy relative to ground state (E/h) in Hz
    elif(isotopemass == 7):
        E_D1 = 446.800132e12; # 7Li 3P1/2 excited state energy relative to ground state (E/h) in Hz
        E_D2 = 446.810184e12; # 7Li 3P3/2 excited state energy relative to ground state (E/h) in Hz
    else:
        raise ValueError("Isotope index invalid!")
        
    # calculating the fields in a reference frame moving with v:
    vabs = np.sqrt(sum(v**2))
    gam = 1/np.sqrt(1-(vabs**2)/c**2)
    if (vabs > 0):
        Lambda = np.matrix([[ gam , -gam*v[0]/c , -gam*v[1]/c , -gam*v[2]/c], # Lorentz boost tensor
                [-gam*v[0]/c , 1+(gam-1)*v[0]**2 / vabs**2 , v[0]*v[1]*(gam-1) / vabs**2 , v[0]*v[2]*(gam-1) / vabs**2], 
                [-gam*v[1]/c , v[0]*v[1]*(gam-1) / vabs**2 , 1+(gam-1)*v[1]**2/ vabs**2 , v[1]*v[2]*(gam-1) / vabs**2],
                [-gam*v[2]/c , v[0]*v[2]*(gam-1) / vabs**2 , v[1]*v[2]*(gam-1) / vabs**2 , 1+(gam-1)*v[2]**2 / vabs**2]])
    else:
        Lambda = sp.dia_matrix((np.ones((4)),0),shape=(4,4)).todense()
        
    F_lab = np.matrix([[0, -E[0]/c, -E[1]/c, -E[2]/c],      # laboratry frame electromagnetic tensor
                   [E[0]/c, 0, -B[2], B[1]],
                   [E[1]/c, B[2], 0, -B[0]],
                   [E[2]/c, -B[1], B[0], 0]])
    
    
    # the moving frame field electromagnetic tensor:
    F_mf = np.matmul(Lambda,(np.conjugate(np.matmul(Lambda,(np.conjugate(F_lab).T))).T))
    B_mf = np.array([F_mf[3,2],-F_mf[3,1],F_mf[2,1]])        # moving frame B field
    B_abs = np.sqrt(sum(B_mf**2))                          # magnitude of B in moving frame
    if (B_abs > 0):
        e_Bmf = B_mf/B_abs                     # unit vector along moving frame B
    else: 
        e_Bmf = np.array([0,0,1])
    
    # Calculating the Zeeman shifted energies:
    # we are working in an LS basis, the labeling of the states is as follows:
    # | g/e , m_L , 2*m_S >
    # e.g. g_m1_p1   means ground state, m_L=-1, m_S=1/2

    # Zeeman shifted ground-state energies:
    E_g_m1 = gS * (-1/2) * B_abs*muB_per_h; 
    E_g_p1 = gS * (1/2) * B_abs*muB_per_h;
    
    # for the excited state shitfs we diagonalize the fine-structure + Zeeman Hamiltonian
    H0 = np.diag(E_D1*np.ones(6));
    H_hf2 = (E_D2-E_D1)*np.array([[1,0,0,0,0,0], [0, 1/3, np.sqrt(2)/3, 0,0,0],
    [0, np.sqrt(2)/3,2/3,0,0,0],[0,0,0,2/3,np.sqrt(2)/3,0], [0,0,0,np.sqrt(2)/3,1/3,0], [0,0,0,0,0,1]])
    H_B = (muB * B_abs / h) * np.diag([ gL*(-1)+gS*(-1/2), gL*(-1)+gS*(1/2),
    gL*(0)+gS*(-1/2), gL*(0)+gS*(1/2), gL*(1)+gS*(-1/2), gL*(1)+gS*(1/2)]);
    H0 = np.matrix(H0)
    H_hf2 = np.matrix(H_hf2)
    H_B = np.matrix(H_B)
    
    ax1 = (H_B+H_hf2); # E/h now in units of Hz
    [vals,vecs] = np.linalg.eig(ax1);
    idx = vals.argsort()#[::-1]   
    vals = vals[idx]
    vecs = vecs[:,idx]
    vecs=np.matrix(vecs)
    energies = vals; # the energies of the excited states, E/h in Hz
    
    freqs = np.zeros(12)
    freqs[0:6] = E_D1 + energies - E_g_m1;
    freqs[6:12] = E_D1 + energies - E_g_p1;
    
    dipoleblock = np.matrix([[1,0,1,0,1,0], [0,1,0,1,0,1]])
    #the nondiagonal block of the dipole operator matrix on the LS basis
    
    matrixelements = np.zeros(12);  # transition matrix elements:
    
    matrixelements[:6] = np.matrix(dipoleblock[0,:])*vecs;
    matrixelements[6:12] = np.matrix(dipoleblock[1,:])*vecs;
    
    # to get the parity of the transitions, we calculate the expectation value of the J_z = L_z + S_z operator and
    # compute delta m_J
    J_z = np.matrix(np.diag([-3/2,-1/2,-1/2,1/2,1/2,3/2])); # the m_J operator is diagonal on our base
    deltam = np.zeros(12);
        
    for i in range(6):
        deltam[i] = (-1/2) - np.transpose(vecs[:,i])*(J_z*vecs[:,i]); # decaying to state |g,0,-1>
        deltam[i+6] = (1/2) - np.transpose(vecs[:,i])*(J_z*vecs[:,i]); # decaying to state |g,0,1>
   
    # and we exclude the forbidden transitions from the lists (which would have delta m_J=+-2)

    indices = np.nonzero(np.abs((np.round(deltam)).astype(int))<=1)[0]
        
    freqs = freqs[indices];
    matrixelements2 = matrixelements[indices];
    deltam = np.round(deltam[indices])
    
    # indices for sigmaplus, pi and sigmaminus polarizations:
    index_sp = np.nonzero(deltam==-1)[0]
    index_pi = np.nonzero(deltam==0)[0]
    index_sm = np.nonzero(deltam==1)[0]
       
    lambdas = c/freqs;
    wavevectors = 2*np.pi/lambdas
    Dopplershifts = sum(n*v)*wavevectors
    freqs = freqs + Dopplershifts;
    
    # finally the detected wavelengths in nm
    flambdas = (1e9/freqs)*c;
    
    # Calculation of the relative intensities:
    costheta = sum(n*e_Bmf);   # angle between B and detector direction 
    sigmapdist = (1/4)*(1+costheta**2); # the angular distributions
    sigmamdist = (1/4)*(1+costheta**2);
    pidist = (1/2)*(1-costheta**2);
    angular_dist = np.zeros(10);
    angular_dist[index_sp] = sigmapdist;
    angular_dist[index_pi] = pidist;
    angular_dist[index_sm] = sigmamdist;
    
    intensities = angular_dist*(matrixelements2**2)/2; 
    # we divide by 2 so that the overall sum is 1
    return flambdas,intensities

def mintavetel(R,M):
    pontokR = [0]
    pontokphi = [0]
    for i in range(M+1):
        if(i > 0):
            for j in range(i):
                pontokR.append(i*R/(M))
                pontokphi.append(2*np.pi*(j/i) + np.pi*i/2)
    return np.array(pontokR),np.array(pontokphi)

def polardescart3d(R,phi):
    return np.array([R*np.cos(phi),R*np.sin(phi),0])
    
def save_spectral_config(lis):
    with open('spectral_fit_config.cfg', 'w') as f:
        for line in range(len(lis)):
            if(type(lis[line]) != str):
                f.write(str(lis[line]))
            else:
                f.write(lis[line])
            f.write('\n')
    return
def load_spectral_config():
    with open('spectral_fit_config.cfg', 'r') as file:
        lis = file.readlines()
        lis = [line.strip() for line in lis]
    return lis
        

def CVI_529_line_generator(grid,roi,B,theta,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,kbt,A,dslit):
    # location where the web service is hosted
    pc_location = 'http://sv-coda-wsvc-28.ipp-hgw.mpg.de:6055'

    # fetching the fine structure of the predefined line
    fine_structure_query = '/getZeeman.json?name=C-VI-5291&B='+str(B)+'&theta1='+str(theta)
    fine_structure = requests.get(pc_location + fine_structure_query).json()
    wl_values = wavelength_grid_generator(grid,wavelength_setting,roi)#loading the wavelength grid
    wl_grid0 = wl_values[wl_values > lower_wl_lim] #slicing the wavelength grid
    wl_grid = wl_grid0[upper_wl_lim > wl_grid0]
    
    #projecting the loaded spectrum into the grid
    projection = np.zeros((wl_grid.shape[0]))
    locations = np.array(fine_structure['wavelengths'])/10
    intensities = np.array(fine_structure['amplitude'])
    locations = locations + mu_add #Doppler shift + calibration uncertainty
    mu = np.dot(locations,intensities)/sum(intensities)
    for i in range(len(intensities)):
        diff = abs(wl_grid - locations[i])
        sor = np.argsort(diff)
        closest_ind = sor[:2]
        distance = diff[sor[:2]]
        I2 = intensities[i] / (1 + distance[1]/distance[0])
        projection[closest_ind[1]] += I2
        projection[closest_ind[0]] += intensities[i] - I2
        
    #addition of the Doppler broadening
    s = doppler_broadening(kbt,mu)
    gaussian = gauss(wl_grid,s,A,upper_wl_lim,lower_wl_lim)
    doppler_spectrum = np.convolve(projection, gaussian, mode = "same")
    
    #convolution with instrumental function
    instr = np.load("instr_funcs/"+grid+"_P0"+str(roi)+"_"+str(int(dslit))+"micron_slit.npy").ravel()
    complete_spectrum = np.convolve(doppler_spectrum, instr, mode = "same")
    
    return complete_spectrum

def Li_line_generator(B,roi,grid,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,A):
    E = np.array([0,0,0])
    centre_of_lens = np.array([1.305624,6.094843,-3.013095])
    rois_positions = np.array([[1.91089568, 5.87253068, 0],
                                [1.91391106, 5.88175452, 0],
                                [1.91353862, 5.89390765, 0],
                                [1.91983473, 5.90358337, 0]]).T
    isomass = 7
    v_direction = rois_positions[:,0] - rois_positions[:,1]
    v_direction = v_direction / np.sqrt(np.sum(v_direction**2))

    m_iso = 1.0
    if(isomass == 6):
        m_iso = 6.015122
    elif(isomass==7):
        m_iso = 7.016003
    v_abs = np.sqrt(16.02177/(1.660539*m_iso))*1e6
    v = v_direction*v_abs
    n = rois_positions[:,roi-1] - centre_of_lens
    spectrum=Li_spectrum (B,E,v,n,isomass)

    #loading the wavelength grid
    wl_values = flap_w7x_abes.wavelength_grid_generator(grid,wavelength_setting,roi)

    #slicing the wavelength grid
    wl_grid0 = wl_values[wl_values > lower_wl_lim]
    wl_grid = wl_grid0[upper_wl_lim > wl_grid0]

    #projecting the loaded spectrum into the grid
    projection = np.zeros((wl_grid.shape[0]))

    locations = spectrum[0]
    intensities = spectrum[1]

    locations = locations + mu_add

    for i in range(len(intensities)):
        diff = abs(wl_grid - locations[i])
        sor = np.argsort(diff)
        closest_ind = sor[:2]
        distance = diff[sor[:2]]
        I2 = intensities[i] / (1 + distance[1]/distance[0])
        projection[closest_ind[1]] += I2
        projection[closest_ind[0]] += intensities[i] - I2

    instr = np.load("instr_funcs/"+grid+"_P0"+str(roi)+".npy").ravel()
    return A*np.convolve(projection, instr, mode = "same")

def Li_line_generator2(B,roi,grid,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,A,R,M):
    E = np.array([0,0,0])
    centre_of_lens = np.array([1.305624,6.094843,-3.013095])
    rois_positions = np.array([[1.91089568, 5.87253068, 0],
                                [1.91391106, 5.88175452, 0],
                                [1.91353862, 5.89390765, 0],
                                [1.91983473, 5.90358337, 0]]).T
    isomass = 7
    v_direction = rois_positions[:,0] - rois_positions[:,1]
    v_direction = v_direction / np.sqrt(np.sum(v_direction**2))

    m_iso = 1.0
    if(isomass == 6):
        m_iso = 6.015122
    elif(isomass==7):
        m_iso = 7.016003
    v_abs = np.sqrt((3/11)*16.02177/(1.660539*m_iso))*1e6
    v = v_direction*v_abs
    n = rois_positions[:,roi-1] - centre_of_lens
    
    K = np.array([-(n[0]+n[1])/n[2],1,1])
    K = K/np.sqrt(sum(K**2))
    L = np.cross(n,K)
    L = L/np.sqrt(sum(L**2))
    minta = mintavetel(R,M)
    szegmensek = np.array(np.zeros((minta[0].shape[0],3)))
    for i in range(minta[0].shape[0]):
        popt = polardescart3d(minta[0][i],minta[1][i])
        szegmensek[i,:] = n + popt[0]*K + popt[1]*L
        
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(minta[1],minta[0],marker = "o",linestyle="")
    raise ValueError("stop")
    
    #loading the wavelength grid
    wl_values = flap_w7x_abes.wavelength_grid_generator(grid,wavelength_setting,roi)

    #slicing the wavelength grid
    wl_grid0 = wl_values[wl_values > lower_wl_lim]
    wl_grid = wl_grid0[upper_wl_lim > wl_grid0]
    projection = np.zeros((wl_grid.shape[0]))
    
    for k in range(szegmensek.shape[0]):
        spectrum=Li_spectrum (B,E,v,szegmensek[k,:],isomass)
        locations = spectrum[0]
        intensities = spectrum[1]
        locations = locations + mu_add
        for i in range(len(intensities)):
            diff = abs(wl_grid - locations[i])
            sor = np.argsort(diff)
            closest_ind = sor[:2]
            distance = diff[sor[:2]]
            I2 = intensities[i] / (1 + distance[1]/distance[0])
            projection[closest_ind[1]] += I2
            projection[closest_ind[0]] += intensities[i] - I2
        
    instr = np.load("instr_funcs/1800g_per_mm_P03_670nm.npy").ravel()
    return A*np.convolve(projection, instr, mode = "same")

def CVI_fitfunc(esti):
    mu_add = esti[0]
    kbt = esti[1]
    A = esti[2]
    
    param = load_spectral_config()
    grid = param[0]
    roi = int(param[1])
    B = float(param[2])
    theta = float(param[3])
    wavelength_setting = float(param[4])
    lower_wl_lim = float(param[5])
    upper_wl_lim = float(param[6])
    dslit = int(param[7])
    
    modelled=CVI_529_line_generator(grid,roi,B,theta,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,kbt,A,dslit)
    # measured = flap.load("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement.dat")
    measured=np.load("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement.npy")
    error=np.load("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement_error.npy")
    C = (modelled - measured)/error
    return (np.dot(C,C) - 3) / measured.shape[0]

def Li_fitfunc(esti):
    mu_add = esti[0]
    A = esti[1]
    B = np.array([esti[2],esti[3],esti[4]])#maybe we should leave that const
    
    param = load_spectral_config()
    grid = param[0]
    roi = int(param[1])
    wavelength_setting = float(param[2])
    lower_wl_lim = float(param[3])
    upper_wl_lim = float(param[4])
    R = float(param[5])
    M = int(param[6])
    
    modelled = Li_line_generator2(B,roi,grid,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,A,R,M)
    measured=np.load("Li_670nm_P0"+str(roi)+"_"+grid+"_measurement.npy")
    error=np.load("Li_670nm_P0"+str(roi)+"_"+grid+"_measurement_error.npy")
    C = (modelled - measured)/error
    return (np.dot(C,C) - 5) / measured.shape[0]

def CVI_fitfunc_plot(measured,measured_err,lambd,mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                     tstop,bg,B,theta,lvl,uvl,tr,dslit,save=False):
    
    modelled = CVI_529_line_generator(grid,roi,B,theta,ws,lvl,uvl,mu_add,kbt,A,dslit)
    fs = 15
    plt.figure()
    plt.errorbar(lambd,measured, measured_err,color = "blue")
    plt.plot(lambd,modelled,color = "red")
    plt.xlabel("Wavelength [nm]",fontsize = fs)
    plt.ylabel("Intensity [a.u.]",fontsize = fs)
    plt.grid()
    plt.legend(["Calculated","Measured"],loc = "best",fontsize = (fs-2))
    return sum(((modelled - measured)/measured_err)**2) / measured.shape[0]

def Li_fitfunc_plot(measured,measured_err,lambd,mu_add,A,expe_id,grid,ws,roi,tstart,
                     tstop,bg,B,lvl,uvl,tr,R,M):
    
    modelled = Li_line_generator2(B,roi,grid,ws,lvl,uvl,mu_add,A,R,M)
    fs = 15
    plt.figure()
    plt.errorbar(lambd,measured, measured_err,color = "blue")
    plt.plot(lambd,modelled,color = "red")
    plt.xlabel("Wavelength [nm]",fontsize = fs)
    plt.ylabel("Intensity [a.u.]",fontsize = fs)
    plt.grid()
    plt.legend(["Calculated","Measured"],loc = "best",fontsize = (fs-2))
    return sum(((modelled - measured)/measured_err)**2) / measured.shape[0]

def Li_Bfit(spectra,mu_add,A,expe_id,grid,ws,roi,tstart,tstop,bg,B,lvl,uvl,tr,R,M):
    
    met="Powell"
    esti = np.array([mu_add,A,B[0],B[1],B[2]])
    measured = active_passive_with_error(spectra,roi,tstart,tstop,expe_id,
                                         tr,bg_wls=bg,plots=False)
    measured = measured.slice_data(slicing={"Wavelength":flap.Intervals(lvl, uvl)})
    save_spectral_config([grid,roi,ws,lvl,uvl,R,M])
    np.save("Li_670nm_P0"+str(roi)+"_"+grid+"_measurement",measured.data)
    np.save("Li_670nm_P0"+str(roi)+"_"+grid+"_measurement_error",measured.error)
        
    lambd = measured.coordinate("Wavelength")[0]
    es_chisq=Li_fitfunc_plot(measured.data,measured.error,lambd,mu_add,A,expe_id,grid,ws,roi,tstart,
                         tstop,bg,B,lvl,uvl,tr,R,M)
    plt.title("$\chi^2 = $"+str(round(es_chisq,6)))
    # raise ValueError("stop")
    # solution=minimize(Li_fitfunc,esti,method=met,tol=1e-8,options={"maxiter":1000})
    # print(solution)
    # if(solution.success == True):
    #     sol=solution.x
    #     # print(solution.hess_inv.matmat(np.eye(3)))
    #     Li_fitfunc_plot(measured.data,measured.error,lambd,sol[0],sol[1],expe_id,grid,ws,roi,tstart,
    #                          tstop,bg,np.array([sol[2],sol[3],sol[4]]),lvl,uvl,tr)

def tempfit_error(sol,fmin,stepsize):
    s = sol[1]
    sl = [s]
    #f = [fmin]
    while(CVI_fitfunc(np.array([sol[0],s,sol[2]])) < 1.0625*fmin):
        s += stepsize
        sl.append(s)
        #f.append(CVI_fitfunc(np.array([sol[0],s,sol[2]])))
        
    sb = abs(s - sol[1])
    s = sol[1]
    while(CVI_fitfunc(np.array([sol[0],s,sol[2]])) < 1.0625*fmin):
        s += -stepsize
        sl.append(s)
        #f.append(CVI_fitfunc(np.array([sol[0],s,sol[2]])))
    sj = abs(s - sol[1])
    return (sj+sb)*2

def tempfit_error_fitfunc(Tfit):
    mu_add,kbt,A,fmin = np.loadtxt("Terror_fit_parameters.txt")
    return abs(CVI_fitfunc(np.array([mu_add,kbt+Tfit[0],A])) - 1.0625*fmin)

def tempfit_error_fit(sol,fmin,met):
    np.savetxt("Terror_fit_parameters.txt",np.array([sol[0],sol[1], sol[2],fmin]))
    Terr_per4_pos = minimize(tempfit_error_fitfunc,np.array([100.0]),method = met,
                         tol=1e-4, options={"maxiter":2000},bounds=((1,None),(None,None)))
    Terr_per4_neg = minimize(tempfit_error_fitfunc,np.array([-100.0]),method = met,
                         tol=1e-4, options={"maxiter":2000},bounds=((None,-1),(None,None)))
    print(Terr_per4_pos)
    print(Terr_per4_neg)
    if(Terr_per4_pos.success == True and Terr_per4_neg.success == True):
        return abs(2*Terr_per4_pos.x[0]) + abs(2*Terr_per4_neg.x[0])

def tempfit_error_curve(sol,stepsize,N):
    res = []
    for i in range(1,N+1):
        res.append(CVI_fitfunc(np.array([sol[0],i*stepsize,sol[2]])))
    plt.figure()
    plt.plot(np.arange(1,N+1,stepsize),np.array(res),"+")
    
def derivative(sol,fmin,h,i): #centered
    sol_p = sol.copy()
    sol_n = sol.copy()
    sol_p[i] = sol_p[i] + h
    sol_n[i] = sol_n[i] - h
    fp1 = CVI_fitfunc(sol_p)
    fm1 = CVI_fitfunc(sol_n)
    return (fp1 - fm1)/(2*h)
    
def second_deriv(sol,fmin,h,i):
    sol_p = sol.copy()
    sol_n = sol.copy()
    sol_p[i] = sol_p[i] + h
    sol_n[i] = sol_n[i] - h
    fp1 = CVI_fitfunc(sol_p)
    fm1 = CVI_fitfunc(sol_n)
    return (fp1 + fm1 - 2*fmin)/h**2

def partial_cross_deriv(sol,fmin,i,j,hi,hj):
    sol_pp = sol.copy()
    sol_pn = sol.copy()
    sol_np = sol.copy()
    sol_nn = sol.copy()
    sol_pp[i] = sol_pp[i] + hi
    sol_pp[j] = sol_pp[j] + hj
    
    sol_pn[i] = sol_pn[i] + hi
    sol_pn[j] = sol_pn[j] - hj
    
    sol_np[i] = sol_np[i] - hi
    sol_np[j] = sol_np[j] + hj
    
    sol_nn[i] = sol_nn[i] - hi
    sol_nn[j] = sol_nn[j] - hj
    
    fpp = CVI_fitfunc(sol_pp)
    fpn = CVI_fitfunc(sol_pn)
    fnp = CVI_fitfunc(sol_np)
    fnn = CVI_fitfunc(sol_nn)
    
    return (fpp-fpn-fnp+fnn) / (4*hi*hj)

    
    
def hesse(sol,fmin,h):
    hess = np.matrix(np.zeros((3,3)))
    for i in range(3):
        for j in range(3):
            if(i == j):
                hess[i,i] = second_deriv(sol,fmin,h[i],i)
            else:
                hess[i,j] = partial_cross_deriv(sol,fmin,i,j,h[i],h[j])
                
    return hess

def error_from_hesse(sol,fmin,h):
    alpha = hesse(sol,fmin,h)#/2
    I = np.linalg.inv(alpha)
    return np.sqrt(abs(I))
    
def CVI_tempfit(spectra,mu_add,kbt,A,expe_id,grid,ws,roi,tstart,tstop,bg,B,theta,lvl,uvl,tr,dslit):
    
    met="Powell"
    h = np.array([1e-5,1e-3,1e-8])
    esti = np.array([mu_add,kbt,A])
    measured = active_passive_with_error(spectra,roi,tstart,tstop,expe_id,
                                         tr,bg_wls=bg,plots=False)
    measured = measured.slice_data(slicing={"Wavelength":flap.Intervals(lvl, uvl)})
    save_spectral_config([grid,roi,B,theta,ws,lvl,uvl,dslit])
    np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement",measured.data)
    np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement_error",measured.error)
        
    lambd = measured.coordinate("Wavelength")[0]
    es_chisq=CVI_fitfunc_plot(measured.data,measured.error,lambd,mu_add,kbt,A,
                              expe_id,grid,ws,roi,tstart,tstop,bg,B,theta,lvl,uvl,tr,dslit,save=True)
    plt.title("$\chi^2 = $"+str(round(es_chisq,6)))
    # raise ValueError("stop")
    solution=minimize(CVI_fitfunc,esti,method=met,bounds = ((None,None),(0.1,None),(None,None)),tol=1e-12,
                                                            options={"maxiter":2000})
    print(solution)
    if(solution.success == True):
        sol=solution.x
        # print(solution.hess_inv.matmat(np.eye(3)))
        CVI_fitfunc_plot(measured.data,measured.error,lambd,sol[0],sol[1],sol[2],
                         expe_id,grid,ws,roi,tstart,tstop,bg,B,theta,lvl,uvl,tr,dslit)
        err = error_from_hesse(sol,solution.fun,h)
        print(err)
        # stepsize = 1
        # err = tempfit_error_fit(sol,solution.fun,met)#tempfit_error(sol,solution.fun,stepsize)
        R_plot = round(spectra.coordinate("Device R")[0][0,(roi-1),0],4)
        # plt.title("R = "+str(R_plot)+" m, $\chi^2 = $"+str(round(solution.fun,6))+", $T_C$ = "+str(round(sol[1],2)))
        plt.title("R = "+str(R_plot)+" m, $\chi^2 = $"+str(round(solution.fun,6))+", $T_C$ = "+str(round(sol[1],2))+" $\pm$ "+str(round(err[1,1],2))+" ev")
        # N = 1000
        # tempfit_error_curve(sol,stepsize,N)
        
def grid_slit_intensity(grid,dslit):
    baseint = 11230
    if(grid == "1200g_per_mm" and dslit == 100):
        return 1
    elif(grid == "1800g_per_mm" and dslit == 100):
        return 16922/baseint
    elif(grid == "2400g_per_mm" and dslit == 100):
        return 10688/baseint
    
    elif(grid == "1200g_per_mm" and dslit == 70):
        return 12810/(baseint*1.5)
    elif(grid == "1800g_per_mm" and dslit == 70):
        return 20885/(baseint*1.5)
    elif(grid == "2400g_per_mm" and dslit == 70):
        return 13971/(baseint*1.5)
    
    elif(grid == "1200g_per_mm" and dslit == 50):
        return 14715/(baseint*2)
    elif(grid == "1800g_per_mm" and dslit == 50):
        return 21279/(baseint*2)
    elif(grid == "2400g_per_mm" and dslit == 50):
        return 14274/(baseint*2)
    
    elif(grid == "1200g_per_mm" and dslit == 35):
        return 12470/(baseint*2.8)
    elif(grid == "1800g_per_mm" and dslit == 35):
        return 21396/(baseint*2.8)
    elif(grid == "2400g_per_mm" and dslit == 35):
        return 16699/(baseint*2.8)
    
    else:
        raise ValueError("Wrong grid or slit size.")
        
def CVI_line_simulator(mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                       tstop,bg,B,theta,lvl,uvl,tr,scalef,errparam,dslit,plots=False):
    measured=flap.load("CVI_529nm_P0"+str(roi)+"_1200g_per_mm_measurement.dat")
    lamb = measured.coordinate("Wavelength")[0]
    gaus = lambda x,A,s,mu : A*np.e**(-(((x-mu)**2)/s**2))
    
    popt, pcov = curve_fit(gaus,lamb, measured.data,p0 = [max(measured.data),0.1,lamb.mean()])
    if(plots == True):
        wl_values = wavelength_grid_generator(grid,ws,roi)#loading the wavelength grid
        wl_grid0 = wl_values[wl_values > lvl] #slicing the wavelength grid
        lambd = wl_grid0[uvl > wl_grid0]
        fs = 15
        plt.figure()
        plt.plot(lamb,measured.data,"+")
        plt.plot(lambd,gaus(lambd,*popt))
        plt.xlabel("Wavelength [nm]",fontsize = fs)
        plt.ylabel("Spectral intensity",fontsize = fs)
        plt.legend(["Data","Fit"],fontsize = fs-2)
        plt.title(expe_id+", Beam on line intensity fit")
    
    gridfac = grid_slit_intensity(grid,dslit)
        
    calculated=CVI_529_line_generator(grid,roi,B,theta,ws,lvl,uvl,mu_add,kbt,A,dslit)
    calculated=calculated/max(calculated)
    calculated = gridfac*scalef*popt[0]*calculated
    
    if(plots == True):
        plt.figure()
        plt.plot(lamb,measured.data,"+")
        plt.plot(lambd,calculated, marker = "o", color = "black")
        plt.grid()
        
    sq = lambda x,a,b : a*np.sqrt(x) + b
    err = sq(calculated*10,errparam[0],errparam[1]) #assuming 10% modulation
    
    if(plots == True):
        plt.figure()
        plt.plot(lamb,measured.data,"o")
        plt.errorbar(lambd,calculated,err)
        plt.grid()
        
    calculated = np.random.normal(loc = calculated, scale=err)
    
    if(plots == True):
        plt.figure()
        plt.plot(lamb,measured.data,"+")
        plt.plot(lambd,calculated, marker = "o", color = "black")
        plt.grid()
        
    return calculated,err

def CVI_line_simulator_me(mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                       tstop,bg,B,theta,lvl,uvl,tr,scalef,dslit,plots=False):
    measured=flap.load("CVI_529nm_P0"+str(roi)+"_1200g_per_mm_measurement.dat")
    lamb = measured.coordinate("Wavelength")[0]
    gaus = lambda x,A,s,mu : A*np.e**(-(((x-mu)**2)/s**2))
    
    popt, pcov = curve_fit(gaus,lamb, measured.data,p0 = [max(measured.data),0.1,lamb.mean()])
    if(plots == True):
        wl_values = wavelength_grid_generator(grid,ws,roi)#loading the wavelength grid
        wl_grid0 = wl_values[wl_values > lvl] #slicing the wavelength grid
        lambd = wl_grid0[uvl > wl_grid0]
        fs = 15
        plt.figure()
        plt.plot(lamb,measured.data,"+")
        plt.plot(lambd,gaus(lambd,*popt))
        plt.xlabel("Wavelength [nm]",fontsize = fs)
        plt.ylabel("Spectral intensity",fontsize = fs)
        plt.legend(["Data","Fit"],fontsize = fs-2)
        plt.title(expe_id+", Beam on line intensity fit")
    
    gridfac = grid_slit_intensity(grid,dslit)
        
    calculated=CVI_529_line_generator(grid,roi,B,theta,ws,lvl,uvl,mu_add,kbt,A,100)
    calculated = gridfac*scalef*popt[0]*calculated/max(calculated)
    
    if(plots == True):
        plt.figure()
        plt.plot(lamb,measured.data,"+")
        plt.plot(lambd,calculated, marker = "o", color = "black")
        plt.grid()
        
    err = measured.error
    
    if(plots == True):
        plt.figure()
        plt.plot(lamb,measured.data,"o")
        plt.errorbar(lambd,calculated,err)
        plt.grid()
        
    calculated = np.random.normal(loc = calculated, scale=err)
    
    if(plots == True):
        plt.figure()
        plt.plot(lamb,measured.data,"+")
        plt.plot(lambd,calculated, marker = "o", color = "black")
        plt.grid()
        
    return calculated,err

def CVI_Ti_error_sim(mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                     tstop,bg,B,theta,lvl,uvl,tr,scalef,errparam,iter_num,dslit,plots = False):
    spectra = get_spectra(expe_id, "1200g_per_mm", ws, roi)
    met="Powell"
    h = np.array([1e-5,1e-3,1e-8])
    line_param = np.array([mu_add,kbt,A])
    measured = active_passive_with_error(spectra,roi,tstart,tstop,expe_id,
                                         tr,bg_wls=bg,plots=False)
    measured = measured.slice_data(slicing={"Wavelength":flap.Intervals(lvl, uvl)})
    save_spectral_config([grid,roi,B,theta,ws,lvl,uvl,dslit])
    np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement",measured.data)
    np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement_error",measured.error)
    
    T_i = np.zeros((iter_num))
    T_i_err = np.zeros((iter_num))
    chisq = np.zeros((iter_num))
    
    for i in range(iter_num):
        print("Iteration "+str(i))
        sim,sim_err = CVI_line_simulator(mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                               tstop,bg,B,theta,lvl,uvl,tr,scalef,errparam,dslit,plots=False)
        np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement",sim)
        np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement_error",sim_err)
        if(plots == True):
            wl_values = wavelength_grid_generator(grid,ws,roi)#loading the wavelength grid
            wl_grid0 = wl_values[wl_values > lvl] #slicing the wavelength grid
            lambd = wl_grid0[uvl > wl_grid0]
            es_chisq=CVI_fitfunc_plot(sim,sim_err,lambd,mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                                      tstop,bg,B,theta,lvl,uvl,tr,dslit,save=True)
            plt.title("$\chi^2 = $"+str(round(es_chisq,6)))
        solution=minimize(CVI_fitfunc,line_param,method=met,bounds = ((None,None),(0.1,None),(None,None)),tol=1e-8,
                                                                options={"maxiter":2000})
        if(solution.success == False):
            raise ValueError("Failed T_i fit")
        # print(solution)
        sol=solution.x
        H = error_from_hesse(sol,solution.fun,h)
        err = H[1,1]
        print(H)
        print(sol[1])
        T_i[i] = sol[1]
        T_i_err[i] = err
        chisq[i] = solution.fun
        if(plots == True):
            CVI_fitfunc_plot(sim,sim_err,lambd,sol[0],sol[1],sol[2],expe_id,grid,ws,roi,
                             tstart,tstop,bg,B,theta,lvl,uvl,tr,dslit)
            R_plot = round(spectra.coordinate("Device R")[0][0,(roi-1),0],4)
            plt.title("R = "+str(R_plot)+" m, $\chi^2 = $"+str(round(solution.fun,6))+", $T_C$ = "+str(round(sol[1],2))+" $\pm$ "+str(round(err,2))+" ev")
            
    print("Average T_i:")
    print(T_i.mean())
    print("STD of T_i:")
    print(np.std(T_i))
    
    print("Average T_i error:")
    print(T_i_err.mean())
    print("STD of T_i error:")
    print(np.std(T_i_err))
    
    print("Average chi square:")
    print(chisq.mean())
    print("STD of chi square:")
    print(np.std(chisq))
    
def CVI_Ti_error_sim_me(mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                     tstop,bg,B,theta,lvl,uvl,tr,scalef,errparam,iter_num,dslit,plots = False):
    spectra = get_spectra(expe_id, "1200g_per_mm", ws, roi)
    met="Powell"
    h = np.array([1e-5,1e-3,1e-8])
    line_param = np.array([mu_add,kbt,A])
    measured = active_passive_with_error(spectra,roi,tstart,tstop,expe_id,
                                         tr,bg_wls=bg,plots=False)
    measured = measured.slice_data(slicing={"Wavelength":flap.Intervals(lvl, uvl)})
    save_spectral_config([grid,roi,B,theta,ws,lvl,uvl,dslit])
    np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement",measured.data)
    np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement_error",measured.error)
    
    T_i = np.zeros((iter_num))
    T_i_err = np.zeros((iter_num))
    chisq = np.zeros((iter_num))
    
    for i in range(iter_num):
        print("Iteration "+str(i))
        sim,sim_err = CVI_line_simulator_me(mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                               tstop,bg,B,theta,lvl,uvl,tr,scalef,dslit,plots=False)
        np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement",sim)
        np.save("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement_error",sim_err)
        if(plots == True):
            wl_values = wavelength_grid_generator(grid,ws,roi)#loading the wavelength grid
            wl_grid0 = wl_values[wl_values > lvl] #slicing the wavelength grid
            lambd = wl_grid0[uvl > wl_grid0]
            es_chisq=CVI_fitfunc_plot(sim,sim_err,lambd,mu_add,kbt,A,expe_id,grid,ws,roi,tstart,
                                      tstop,bg,B,theta,lvl,uvl,tr,dslit,save=True)
            plt.title("$\chi^2 = $"+str(round(es_chisq,6)))
        solution=minimize(CVI_fitfunc,line_param,method=met,bounds = ((None,None),(0.1,None),(None,None)),tol=1e-8,
                                                                options={"maxiter":2000})
        if(solution.success == False):
            raise ValueError("Failed T_i fit")
        # print(solution)
        sol=solution.x
        H = error_from_hesse(sol,solution.fun,h)
        err = H[1,1]
        print(H)
        print(sol[1])
        T_i[i] = sol[1]
        T_i_err[i] = err
        chisq[i] = solution.fun
        if(plots == True):
            CVI_fitfunc_plot(sim,sim_err,lambd,sol[0],sol[1],sol[2],expe_id,grid,ws,roi,
                            tstart,tstop,bg,B,theta,lvl,uvl,tr,dslit)
            R_plot = round(spectra.coordinate("Device R")[0][0,(roi-1),0],4)
            plt.title("R = "+str(R_plot)+" m, $\chi^2 = $"+str(round(solution.fun,6))+", $T_C$ = "+str(round(sol[1],2))+" $\pm$ "+str(round(err,2))+" ev")
            
    print("Average T_i:")
    print(T_i.mean())
    print("STD of T_i:")
    print(np.std(T_i))
    
    print("Average T_i error:")
    print(T_i_err.mean())
    print("STD of T_i error:")
    print(np.std(T_i_err))
    
    print("Average chi square:")
    print(chisq.mean())
    print("STD of chi square:")
    print(np.std(chisq))
    
def minim(K0,S,sigma,M,H): #additional function for the following deconv fun
    K=K0
    M=np.array(M)
    N=np.zeros((H.shape[1],H.shape[1]))
    for k in range(H.shape[1]):
        for i in range(H.shape[1]):
            N[k,i] = H[k,i] + K*sum(M[:,k]*M[:,i]/sigma**2)
    L=K*M.T/sigma**2
    N=np.matrix(N)
    L=np.matrix(L)
    perm=L*S
    M=np.matrix(M)
    P=np.linalg.inv(N)*perm
    X=(np.array(S-M*P))**2/sigma**2
    print(sum(X[:,0])/X.shape[0])
    print(K)
    return P,M
    
def spectrum_deconv(spectra,grid,roi,lvl,uvl,n,K):
    instr = np.load("instr_funcs/"+grid+"_P0"+str(roi)+".npy").ravel()
    spectra = spectra.slice_data(slicing={"Wavelength":flap.Intervals(lvl, uvl)})
    S=spectra.data
    sigma=spectra.error
    instr=np.load("instr_funcs/"+grid+"_P0"+str(roi)+".npy").ravel()
    instr=instr-min(instr)
    instr=instr/max(instr) 
    res = spectra.coordinate("Wavelength")[0][1]-spectra.coordinate("Wavelength")[0][0]
    
    lambd0=np.linspace(0,instr.shape[0],instr.shape[0])
    lambd02=np.linspace(0,instr.shape[0],instr.shape[0]*n)
    f=interpolate.interp1d(lambd0,instr)
    A=f(lambd02)
    M0=np.zeros((n*S.shape[0],n*S.shape[0]))
    for i in range(n*S.shape[0]):
        M0[i,i] = 1
        M0[i,:]=np.convolve(M0[i,:],A,"same")
    M=np.zeros((S.shape[0],n*S.shape[0]))
    
    for i in range(S.shape[0]):
        for k in range(n*S.shape[0]):
            if(abs(i - k/n) < abs(i + 1 - k/n) and abs(i - k/n) < abs(i - 1 - k/n)):
                M[i,:] += M0[k,:]
            elif(abs(i - k/n) == abs(i + 1 - k/n) and i+1 < S.shape[0]):
                M[i,:] += M0[k,:]/2
                M[i+1,:] += M0[k,:]/2
    M[0,:] = 2*M[0,:] 
    
    H0=np.zeros((n*S.shape[0],n*S.shape[0]))
    for i in range(n*S.shape[0]):
        for j in range(n*S.shape[0]):
            if(i==j):
                H0[i,j]=n/res
                if(j < n*S.shape[0]-1):
                    H0[i,j+1] = -n/res
    H = np.dot(H0,H0.T)
    S=np.matrix(S).T
    mtx=minim(K,S,sigma,M,H)
    
    lambd=np.linspace(lvl,uvl,S.shape[0])
    lambd2=np.linspace(lvl,uvl,H.shape[0])
    
    plt.figure()
    plt.plot(lambd,S,color="blue")
    plt.plot(lambd,np.array(mtx[1]*mtx[0]),color="red")
    plt.xlabel("$\lambda [nm]$",fontsize=15)
    plt.ylabel("$Intensity$",fontsize=15)
    plt.title("20181018_010, channel 43, m/n = "+str(n)+",K="+str(K),fontsize=15)
    plt.grid()
    
    plt.figure()
    plt.plot(lambd2,np.array(mtx[0]),color="green")
    plt.xlabel("$\lambda [nm]$",fontsize=15)
    plt.ylabel("$Intensity$",fontsize=15)
    plt.title("20181018_010, channel 43, m/n = "+str(n)+",K="+str(K),fontsize=15)
    plt.grid()