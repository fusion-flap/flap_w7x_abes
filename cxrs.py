# -*- coding: utf-8 -*-
"""
Created on Tue Aug 8 15:27:49 2018

@author: bcsillag

Data processing code for Wendelstein 7-X QSI CXRS spectra measured during OP2.1
"""

import datetime
import requests
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import flap
import flap_w7x_abes
    
def parab(x,a,b,c):
    return a*x**2 + b*x + c

def gauss(lambd,s,A,up,down):
    centr = (up + down) / 2
    return A*np.e**(-(((lambd-centr)**2)/s**2)) / s

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
    
def slice_by_wl_range(spectra,wstart,wstop,expe_id,roi):
    spectra_1w = spectra.slice_data(slicing={"ROI" :"P0"+str(roi)})
    spectra_1w = spectra_1w.slice_data(slicing={"Wavelength" : flap.Intervals(wstart, wstop)},summing = {"Wavelength":"Mean"})
    plt.figure()
    spectra_1w.plot(axes="Time")
    plt.xlabel("Time [s]",fontsize = 15)
    plt.ylabel("Intensity",fontsize = 15)
    plt.title("Time evolution of a pixel at wavelength ["+str(wstart)+", "+str(wstop)+"] nm")
    
    c = spectra_1w.get_coordinate_object("Time")
    timeres = np.mean(c.values[1:(len(c.values)-1)]-c.values[0:(len(c.values)-2)])
    t0 = c.values[0]
    c.mode.equidistant = True
    c.shape = spectra_1w.data.shape
    c.start = t0
    c.step = timeres
    c.dimension_list = [0]
    
    plt.figure()
    ap = spectra_1w.apsd(options={'Interval':1, 'Range':[1,100],'Res':0.1})   
    ap.plot(axes="Frequency",options={'Log x': True, 'Log y':True, 'Error':True, 'X range':[1,100]})
    plt.title(" ")
    
def active_passive(spectra,roi,t_start,t_stop,el,expe_id,timerange):
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
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Sum"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Sum"})
    
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
    return np.sqrt(spec_perint.var(axis = 1) / (spec.data.shape[2]))

def get_line_intensity(spectra,roi,t_start,t_stop,el,expe_id,timerange,lstart,lstop,bg_wls=[0],plots=False):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    #correcting the timescales
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    ROI1 = ROI1.slice_data(slicing={"Wavelength" :flap.Intervals(lstart, lstop)})
    s_on_sliced=ROI1.slice_data(slicing={'Time':d_beam_on}) #for error calculation
    s_off_sliced=ROI1.slice_data(slicing={'Time':d_beam_off})
    
    
    if(bg_wls != [0]):
        ROI1_witbg = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),
                                        "Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                        summing = {"Wavelength":"Mean","Time":"Mean"})
        ROI1.data[:,:] = ROI1.data[:,:] - ROI1_witbg.data
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Sum"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Sum"})
    
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
    
def active_passive_with_error(spectra,roi,t_start,t_stop,el,expe_id,timerange,bg_wls=[0],plots=False):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    #correcting the timescales
    c=d_beam_on.get_coordinate_object("Time")
    c.start = c.start + el
    c=d_beam_off.get_coordinate_object("Time")
    c.start = c.start + el
    
    #slicing the data
    ROI1 = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),"Time":flap.Intervals(t_start, t_stop)})
    s_on_sliced=ROI1.slice_data(slicing={'Time':d_beam_on}) #for error calculation
    s_off_sliced=ROI1.slice_data(slicing={'Time':d_beam_off})
    
    #the intervals are taken as independent measurements
    error_on = indep_spectral_error_calc(s_on_sliced)
    error_off = indep_spectral_error_calc(s_off_sliced)
    
    if(bg_wls != [0]):
        ROI1_witbg = spectra.slice_data(slicing={"ROI" :"P0"+str(roi),
                                        "Wavelength":flap.Intervals(bg_wls[0], bg_wls[1])},
                                        summing = {"Wavelength":"Mean","Time":"Mean"})
        ROI1.data[:,:] = ROI1.data[:,:] - ROI1_witbg.data
    s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
    s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
    
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

def error_distr(spectra,roi,t_start,t_stop,el,expe_id,minint,timerange,lstart,lstop,bg_wls=[0],plots=False):
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    #correcting the timescales
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
    
def tshif(qsi_cxrs,roi,t_start,t_stop,expe_id,timerange,wstart,wstop,N):
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
        s_on=ROI1.slice_data(slicing={'Time':d_beam_on},summing = {"Rel. Time in int(Time)":"Mean"})
        s_off=ROI1.slice_data(slicing={'Time':d_beam_off},summing = {"Rel. Time in int(Time)":"Mean"})
        
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
    
    
def normed_intensity(qsi_cxrs,t_start,t_stop,expe_id,timerange,el,bg_wls):
    print("***** Reading beam-on, beam-off intervals into DataObjects.")
    d_beam_on=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 0, 'Defl': 0}},\
                             object_name='Beam_on',
                             coordinates={'Time': timerange})

    d_beam_off=flap.get_data('W7X_ABES',exp_id=expe_id,name='Chopper_time',
                             options={'State':{'Chop': 1, 'Defl': 0}},\
                             object_name='Beam_off',
                             coordinates={'Time': timerange})
    
    flap.list_data_objects(d_beam_on)
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
    flap.list_data_objects(s_subs)
    
    lineint = []
    int_calib = np.array([5961,5667,1209,3367])/5961
    for i in range(4):
        line = max(s_subs.data[i,200:])
        lineint.append(line)
    lineint = np.array(lineint)
    print(lineint/int_calib)
    
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
        

def CVI_529_line_generator(grid,roi,B,theta,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,kbt,A):
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
    instr = np.load("instr_funcs/"+grid+"_P0"+str(roi)+".npy").ravel()
    complete_spectrum = np.convolve(doppler_spectrum, instr, mode = "same")
    
    return complete_spectrum

def CVI_fitfunc(esti):
    mu_add = esti[0]
    kbt = esti[1]
    A = esti[2]
    
    param = load_spectral_config()
    grid = param[0]
    roi = int(param[1])
    B = float(param[2])
    theta = int(param[3])
    wavelength_setting = float(param[4])
    lower_wl_lim = float(param[5])
    upper_wl_lim = float(param[6])
    
    modelled = CVI_529_line_generator(grid,roi,B,theta,wavelength_setting,lower_wl_lim,upper_wl_lim,mu_add,kbt,A)
    measured = flap.load("CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement.dat")
    C = (modelled - measured.data)/measured.error
    return (np.dot(C,C) - 3) / measured.data.shape[0]

def CVI_fitfunc_plot(spectra,mu_add,kbt,A,expe_id,grid,ws,roi,tshift,tstart,
                     tstop,bg,B,theta,lvl,uvl,tr,save=False):
    
    modelled = CVI_529_line_generator(grid,roi,B,theta,ws,lvl,uvl,mu_add,kbt,A)
    measured = active_passive_with_error(spectra,roi,tstart,tstop,tshift,expe_id,
                                         tr,bg_wls=bg,plots=False)
    measured = measured.slice_data(slicing={"Wavelength":flap.Intervals(lvl, uvl)})
    if(save == True):
        save_spectral_config([grid,roi,B,theta,ws,lvl,uvl])
        flap.save(measured,"CVI_529nm_P0"+str(roi)+"_"+grid+"_measurement.dat")
    lambd = measured.coordinate("Wavelength")[0]
    fs = 15
    plt.figure()
    plt.errorbar(lambd,measured.data, measured.error,color = "blue")
    plt.plot(lambd,modelled,color = "red")
    plt.xlabel("Wavelength [nm]",fontsize = fs)
    plt.ylabel("Intensity [a.u.]",fontsize = fs)
    plt.grid()
    plt.legend(["Calculated","Measured"],loc = "best",fontsize = (fs-2))
    return sum(((modelled - measured.data)/measured.error)**2) / measured.data.shape[0]

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
    
def CVI_tempfit(spectra,mu_add,kbt,A,expe_id,grid,ws,roi,tshift,tstart,tstop,bg,B,theta,lvl,uvl,tr):
    
    met="Powell"
    esti = np.array([mu_add,kbt,A])
    es_chisq=CVI_fitfunc_plot(spectra,mu_add,kbt,A,expe_id,grid,ws,roi,tshift,tstart,
                              tstop,bg,B,theta,lvl,uvl,tr,save=True)
    plt.title("$\chi^2 = $"+str(round(es_chisq,6)))
    solution=minimize(CVI_fitfunc,esti,method=met,bounds = ((None,None),(0.1,None),(None,None)),tol=1e-8,
                                                            options={"maxiter":2000})
    print(solution)
    if(solution.success == True):
        sol=solution.x
        CVI_fitfunc_plot(spectra,sol[0],sol[1],sol[2],expe_id,grid,ws,roi,
                         tshift,tstart,tstop,bg,B,theta,lvl,uvl,tr)
        # stepsize = 1
        err = tempfit_error_fit(sol,solution.fun,met)#tempfit_error(sol,solution.fun,stepsize)
        R_plot = round(spectra.coordinate("Device R")[0][0,(roi-1),0],4)
        plt.title("R = "+str(R_plot)+" m, $\chi^2 = $"+str(round(solution.fun,6))+", $T_C$ = "+str(round(sol[1],2))+" $\pm$ "+str(round(err,2))+" ev")
        # N = 1000
        # tempfit_error_curve(sol,stepsize,N)