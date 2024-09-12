# -*- coding: utf-8 -*-
"""
Created on Tue Aug 8 15:27:49 2023

@author: bcsillag

Data processing code for Wendelstein 7-X QSI CXRS spectra measured during OP2.1
It can be modified easely to process data from later campaigns as well
"""

import requests
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import pandas as pd
import flap


def wavelength_grid_generator_op21(grid, w_s, roi,datapath_base):
    """
    Generates a wavelength grid for spectra measured by IsoPlane using 
    coefficients that were fitted by calibration (which assumed third-order
    polynomial relation between the wavelength values and the pixels). This
    wavelength grid is only valid for the given ROI (channel).
    
    INPUT:
        grid: "1200g_per_mm", "1800g_per_mm" or "2400g_per_mm"
        w_s: central wavelength setting of the measurement (float)
        roi: Region Of Interest (ROI), which belongs to the channel in question (int)
        datapath_base : string, path of the supplementary data (like wavelength calibration)
        
    OUTPUT:
        1D numpy array with the calibrated wavelength values
    """
    
    calib_array = np.loadtxt(datapath_base+"wavelength_calib_2023_"+grid+".txt")
    # loading the calibration coeffs
    c0 = calib_array[roi-1, 0]
    c1 = calib_array[roi-1, 2]
    c2 = calib_array[roi-1, 4]
    c3 = calib_array[roi-1, 6]  # writing the out to clear variables
    pix_values = np.arange(0, 1024, 1) - 511 #centralizing the pixel range
    return c3*pix_values**3 + c2*pix_values**2 + c1*pix_values + c0 + w_s

def wavelength_grid_generator_op22(grid, w_s,datapath_base):
    """
    Generates a wavelength grid for spectra measured by IsoPlane using 
    coefficients that were fitted by calibration (which assumed second-order
    polynomial relation between the wavelength values and the pixels).
    
    INPUT:
        grid: "1200g_per_mm", "1800g_per_mm" or "2400g_per_mm"
        w_s: central wavelength setting of the measurement (float)
        roi: Region Of Interest (ROI), which belongs to the channel in question (int)
        datapath_base : string, path of the supplementary data (like wavelength calibration)
        
    OUTPUT:
        1D numpy array with the calibrated wavelength values
    """
    
    calib_array = np.loadtxt(datapath_base+"wavelength_calib_2024_"+grid+".txt")
    # loading the calibration coeffs
    c0 = calib_array[0]
    c1 = calib_array[2]
    c2 = calib_array[4] # writing the out to clear variables
    pix_values = np.arange(0, 1024, 1) - 511 #centralizing the pixel range
    return c2*pix_values**2 + c1*pix_values + c0 + w_s


def interval_shift(expe_id):
    """
    A function to correct the discrepancies between the clock of the CMOS and
    the spectrometer camera. It can be determined using the spectra.tshif()
    function.
    
    Parameters
    ----------
    expe_id : string, experiment ID or shot ID

    Returns
    -------
    shift: temporal shift between the two clocks (float).

    """
    shift = 0
    if(expe_id[:8] == "20230314"):
        shift = -0.0509
    elif(expe_id[:8] == "20230314.025"):
        shift = -0.04924242424242424
    elif(expe_id[:8] == "20230315"):
        shift = -0.05087939698492462
    if(expe_id[:8] == "20230316"):
        shift = -0.05012562814070352
    elif(expe_id == "20230316.043"):
        shift = -0.037
    elif(expe_id == "20230316.047"):
        shift = -0.058466933867735466
    elif(expe_id == "20230316.072"):
        shift = -0.05012562814070352
    elif(expe_id[:8] == "20230323"):
        shift = -0.075
        # shift = -0.07229458917835671#-0.0592
    elif(expe_id[:8] == "20230328"):
        shift = 0.051633
    elif(expe_id[:8] == "20230330"):
        shift = 0.05275551102204408  # 0.051633

    return shift

def spectral_error_calc_op21(spec):
    """
    Intensity error calculation function for spectra averaged over intervals.
    It calculates the uncertainty based on the temporal variance of the
    intensity, and it takes a beam interval as one measurement.

    Parameters
    ----------
    spec : spectra class object.

    Returns
    -------
    A 2D numpy array with the errors.

    """
    spec_perint = np.zeros((spec.data.shape[0], spec.data.shape[2]))
    for i in range(spec.data.shape[2]):
        ind = np.nonzero(abs(spec.data[10, :, i]))
        for j in range(1, max(ind[0]+1)):
            spec_perint[:, i] = spec_perint[:, i] + \
                (spec.data[:, j, i] / (len(ind[0])-1))
    for i in range(spec.data.shape[2]):
        if(sum(spec_perint[:, i]) == 0 and i > 0):
            spec_perint[:, i] = spec_perint[:, i-1]
        if(sum(spec_perint[:, i]) == 0 and i == 0):
            spec_perint[:, i] = spec_perint[:, i+1]

    return np.sqrt(spec_perint.var(axis=1) / (spec.data.shape[2]))


def indep_spectral_error_calc_op21(spec):
    """
    Intensity error calculation function for spectra averaged over intervals.
    It calculates the uncertainty based on the INDEPENDENT temporal variance of 
    the intensity - which means it does not consider the joint changing of the
    spectral line as error. It takes a beam interval as one measurement.

    Parameters
    ----------
    spec : spectra class object.

    Returns
    -------
    A 2D numpy array with the errors.
    """
    spec_perint = np.zeros((spec.data.shape[0], spec.data.shape[2]))
    for i in range(spec.data.shape[2]):
        ind = np.nonzero(abs(spec.data[0, :, i]) > 1e-10)
        for j in range(1, max(ind[0]+1)):
            spec_perint[:, i] = spec_perint[:, i] + \
                (spec.data[:, j, i] / (len(ind[0])-1))
    for i in range(spec.data.shape[2]):
        if(sum(spec_perint[:, i]) == 0 and i > 0):
            spec_perint[:, i] = spec_perint[:, i-1]
        if(sum(spec_perint[:, i]) == 0 and i == 0):
            spec_perint[:, i] = spec_perint[:, i+1]

    HR = np.zeros((spec.data.shape[0]))
    HL = np.zeros((spec.data.shape[0]))
    for i in range(1, spec_perint.shape[0]-1):
        M1 = spec_perint[i, :]
        HR[i] = np.mean(M1**2)-np.mean(M1*spec_perint[i+1, :]) * \
            M1.mean()/spec_perint[i+1, :].mean()
        HL[i] = np.mean(M1**2)-np.mean(M1*spec_perint[i-1, :]) * \
            M1.mean()/spec_perint[i-1, :].mean()
    H = (HR + HL) / 2
    M1 = spec_perint[0, :]
    H[0] = np.mean(M1**2)-np.mean(M1*spec_perint[1, :]) * \
        M1.mean()/spec_perint[1, :].mean()
    M1 = spec_perint[-1, :]
    H[-1] = np.mean(M1**2)-np.mean(M1*spec_perint[-2, :]) * \
        M1.mean()/spec_perint[-2, :].mean()
    err = np.sqrt(abs(H) / (spec_perint.shape[1]))
    
    return err

class spectra:
    """
    A class for storing and manipulating measured spectra. At the definition
    one has to add the following parameters:
        data_source: string, data source (usually 'W7X_WEBAPI')
        expe_id: string, experiment ID or shot ID
    The optional parameters are:
        get_data: string, the way to reach data for the object - default: "by shotID"
        campaign: string, label of the operation phase - default: "OP2.1"
        spatcal: Boolean, wether or not add spatial coordinates to the 
                 dataobject - default: False
        time_correction: Boolean, whether or not perform correction on the
                         time axis using interval_shift()
    """

    def __init__(self, data_source, expe_id, get_data="by shotID",
                 campaign="OP2.1", spatcal=False, time_correction=False):
        self.data_source = data_source
        self.expe_id = expe_id
        self.campaign = campaign
        self.d_beam_on = None
        self.d_beam_off = None
        self.grid = None
        self.wavelength_setting = None
        self.observation_directions = None
        self.B = None
        self.current_roi = None
        self.Babs = None
        self.current_theta = None
        self.zc_locations = None
        self.zc_intensities = None
        self.dslit = None
        self.wstart = None
        self.wstop = None
        self.active_spectrum = None
        self.simulated = None
        self.simulated_error = None
        self.simd = None
        self.simgrid = None
        self.errparam = None
        self.supl_data_path = None
        self.instr_funcs_datapath = None

        #loading data
        if(data_source == 'W7X_WEBAPI' and get_data == "by shotID" and
                campaign == "OP2.1"):
            self.dataobj = flap.get_data(data_source,
                            name='Test/raw/W7X/QSI/cxrs_DATASTREAM/0/Images/',
                            exp_id=expe_id,
                            options={'Scale Time': True,'Cache Data': False},
                            object_name='QSI_spectral_data')
        if(data_source == 'W7X_WEBAPI' and get_data == "by shotID" and
                campaign == "OP2.2"):
            self.dataobj = flap.get_data(data_source,
                            name='Test/raw/W7X/QSI/cxrs_DATASTREAM/0/Images/',
                            exp_id=expe_id,
                            options={'Scale Time': True,'Cache Data': False},
                            object_name='QSI_spectral_data')
        else:
            raise ValueError("Undefined data source or loading process.")

        if(self.campaign == "OP2.1"):
            self.APDCAM_timerange = [1, 40]

        if(spatcal == True and self.campaign == "OP2.1"):
            # changing to ROI coord
            roi_coord_flap = flap.Coordinate(name='ROI', unit="1",
                            mode=flap.CoordinateMode(equidistant=False),
                            shape=(6),values=["P01", "P02", "P03",
                            "P04", "P05", "P06"],dimension_list=[1])
            self.dataobj.add_coordinate_object(roi_coord_flap)
            self.dataobj.del_coordinate("Coord 1")

            # adding coordinate Device R
            R_coord_flap = flap.Coordinate(name='Device R', unit="m",
                        mode=flap.CoordinateMode(equidistant=False), shape=(6),
                        values=[6.175608381963992,6.18531258540187,6.196755395927777,
                        6.207903188440903,6.22187706429289,6.241915857714211],
                        dimension_list=[1])
            self.dataobj.add_coordinate_object(R_coord_flap)

            # correcting a systematic error on the time coordinate
            c = self.dataobj.get_coordinate_object("Time")
            c.values = c.values + 1
            c = self.dataobj.get_coordinate_object("Time")
            
        if(spatcal == True and self.campaign == "OP2.2"):
             # changing to ROI coord
             
             roi_coord_flap = flap.Coordinate(name='ROI', unit="1",
                             mode=flap.CoordinateMode(equidistant=False),
                             shape=(54),values=np.arange(1,55,1),dimension_list=[1])
             self.dataobj.add_coordinate_object(roi_coord_flap)
             self.dataobj.del_coordinate("Coord 1")

             # adding coordinate Device R
             R_coord_flap = flap.Coordinate(name='Device R', unit="m",
                         mode=flap.CoordinateMode(equidistant=False), shape=(54),
                         values=[6.29246947, 6.28930728, 6.28203192, 6.27886972, 6.27254534, 6.26907393,
                          6.26591175, 6.26210779, 6.2589456, 6.25547416, 6.25262121, 6.24945902,
                          6.24598759, 6.24313463, 6.23997244, 6.23585928, 6.23333883, 6.23048586,
                          6.2263727,  6.22321051, 6.22004831, 6.21657687, 6.21372392, 6.21056173,
                          6.2070903,  6.20423734, 6.20012419, 6.196962,   6.19475076, 6.19063761,
                          6.18716615, 6.18431322, 6.18146043, 6.17798883, 6.17482664, 6.16881164,
                          6.16534006, 6.16248724, 6.15837414, 6.15585347, 6.15204974, 6.14857813,
                          6.26805975, 6.27139163, 6.24371252, 6.24447035, 6.22062642, 6.22043246,
                          6.19658979], #must be updated!!!
                         dimension_list=[1])
             self.dataobj.add_coordinate_object(R_coord_flap)

             # correcting a systematic error on the time coordinate
             c = self.dataobj.get_coordinate_object("Time")
             c.values = c.values + 1
             c = self.dataobj.get_coordinate_object("Time")

        elif(spatcal == True and self.campaign != "OP2.1"):
            raise ValueError(
            "Spatial calibration for that campaign has not been implemented yet.")

        if(self.campaign == "OP2.1" and time_correction == True):
            #getting beam intervals
            d_beam_on = flap.get_data('W7X_ABES', exp_id=self.expe_id, 
                                      name='Chopper_time',options={
                                          'State': {'Chop': 0, 'Defl': 0}},
                                      object_name='Beam_on',
                                      coordinates={'Time': self.APDCAM_timerange})

            d_beam_off = flap.get_data('W7X_ABES', exp_id=self.expe_id, 
                                       name='Chopper_time', options={
                                           'State': {'Chop': 1, 'Defl': 0}},
                                       object_name='Beam_off',
                                       coordinates={'Time': self.APDCAM_timerange})

            # correcting the timescales
            el = interval_shift(self.expe_id)
            c = d_beam_on.get_coordinate_object("Time")
            c.start = c.start + el
            c = d_beam_off.get_coordinate_object("Time")
            c.start = c.start + el

            self.d_beam_on = d_beam_on
            self.d_beam_off = d_beam_off

    def wavelength_calibration(self,man=False,grid=None,wavelength_setting=None,
                               options=None):
        """
        Changes the pixel coordinate axis to wavelength.
        
        INPUT:
            man: Boolean, whether or not one wants to give the calibration settings
                 manually
            grid: "1200g_per_mm", "1800g_per_mm" or "2400g_per_mm", but it only
                  matters when man = True - the default is None
            wavelength_setting: floot, central wavelength setting of the measurement,
                  but it only matters when man = True - the default is None
            options : dictionary
                'Supplementary_data_path': string, path of the supplementary data
                (like wavelength calibration)
        OUTPUT:
            1D numpy array with the calibrated wavelength values
        """
        default_options = {"Supplementary_data_path":"data"}
        _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES_CXRS')
        try:
            datapath_base = _options['Supplementary_data_path']
        except (KeyError, TypeError):
            datapath_base = 'data'
        if(self.campaign == "OP2.1" and man == False): #for OP2.1 a table contains the values
            settings_df = pd.read_excel(
                datapath_base+"discharge_table_for_spectral_measurements.xlsx")
            settings = settings_df.loc[settings_df["Discharge"] == float(
                self.expe_id)]
            self.grid = str(settings["Grid [g/mm]"].to_numpy()[0])+"g_per_mm"
            self.wavelength_setting = float(
                settings["Wavelength setting [nm]"])

        if(man == True): #these can also be added manually for other measurements
            self.grid = grid
            self.wavelength_setting = wavelength_setting
            
        if(self.campaign == "OP2.1"):
            wl_values = np.zeros((6, 1024))
            for i in range(1, 7):
                wl_values[i-1, :] = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, i,datapath_base)

            wl_coord_flap = flap.Coordinate(name='Wavelength', unit="nm",
                            mode=flap.CoordinateMode(equidistant=False),
                            shape=(6,1024),values=wl_values,dimension_list=[1,2])
            self.dataobj.add_coordinate_object(wl_coord_flap)
            self.dataobj.del_coordinate("Coord 2")
            
        if(self.campaign == "OP2.2"):
            wl_values = np.zeros((54, 1024))
            for i in range(1, 54):
                wl_values[i-1, :] = wavelength_grid_generator_op22(
                    self.grid, self.wavelength_setting, datapath_base)

            wl_coord_flap = flap.Coordinate(name='Wavelength', unit="nm",
                            mode=flap.CoordinateMode(equidistant=False),
                            shape=(54,1024),values=wl_values,dimension_list=[1,2])
            self.dataobj.add_coordinate_object(wl_coord_flap)
            self.dataobj.del_coordinate("Coord 2")
            
        elif(self.campaign != "OP2.1"):
            raise ValueError(
            "Wavelength calibration for that campaign has not been implemented yet.")
                

    def slice_by_wl(self, roi, wavelength):
        """
        It slices the temporal signal out in the corresponding wavelength, averages
        it if an interval was given, then it plots the result.
        
        INPUT:
            roi: Region Of Interest (ROI), in other words spectral channel (int)
            wavelength: at which the changes over time are important (float).
                        If it is a list two floats, then the function first slices,
                        then averages the signal between the two wavelengths.
        """
        spectra_1w = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi)})
        if(type(wavelength) == float or type(wavelength) == int):
            spectra_1w = spectra_1w.slice_data(
                slicing={"Wavelength": wavelength})
        elif(type(wavelength) == list):
            spectra_1w = spectra_1w.slice_data(
                slicing={"Wavelength": flap.Intervals(
                    wavelength[0], wavelength[1])},
                summing={"Wavelength": "Mean"})
        plt.figure()
        spectra_1w.plot(axes="Time")
        plt.xlabel("Time [s]", fontsize=15)
        plt.ylabel("Intensity", fontsize=15)
        if(type(wavelength) == float or type(wavelength) == int):
            plt.title(
                "Temporal evolution of the pixel at wavelength "+str(wavelength)+" nm")
        elif(type(wavelength) == list):
            tit = "Avg. temporal evolution of pixels between wavelength ["
            plt.title(tit+str(wavelength[0])+", "+str(wavelength[1])+"] nm")

    def passive(self, roi, t_start, t_stop, error=False):
        """
        It averages the spectra in the given time interval, then it plots the result.
        
        INPUT:
            roi: Region Of Interest (ROI), in other words spectral channel (int)
            t_start: beginning of the interval (float)
            t_stop: end of the interval (float)
            
        """
        ROI1 = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                                "Time": flap.Intervals(t_start, t_stop)})
        errors = None
            
        avg_spectrum = ROI1.slice_data(summing={"Time": "Mean"})
        errors0 = ROI1.data.copy()
        if(error == True):
            for i in range(ROI1.data.shape[1]):
                errors0[:,i] = (errors0[:,i] - avg_spectrum.data[i])**2
                errors = np.sqrt(errors0.mean(axis = 0)) / np.sqrt(errors0.shape[0])
            avg_spectrum.error = errors
        plt.figure()
        avg_spectrum.plot(axes="Wavelength")
        plt.title(self.expe_id, fontsize=15)
        plt.title(self.expe_id+", averaged spectra, ROI = P0" +
                  str(roi), fontsize=15)
        plt.grid()
        return avg_spectrum

    def autocorr(self, roi, t_start, t_stop, lstart, lstop):
        """
        It averages the spectra in the given wavelength interval, slices them into 
        the given time interval, then calculates the temporal autocorrelation of the 
        signal. Finally, it plots the result.
        
        INPUT:
            roi: Region Of Interest (ROI), in other words spectral channel (int)
            t_start: beginning of the time interval (float)
            t_stop: end of the time interval (float)
            lstart: beginning of the wavelength interval (float)
            lstop: end of the wavelength interval (float)
        """
        
        #slicing
        ROI1 = self.dataobj.slice_data(
            slicing={"ROI": "P0"+str(roi), "Time": flap.Intervals(t_start, t_stop)})
        line = ROI1.slice_data(slicing={"Wavelength": flap.Intervals(lstart, lstop)},
                               summing={"Wavelength": "Sum"})

        #calculating the correlation after normalization
        norm_sig=(line.data-np.mean(line.data))/(np.std(line.data)*np.sqrt(len(line.data)))
        corr = np.correlate(norm_sig, norm_sig, mode="same")
        
        #time axis for the correlation
        t = ROI1.coordinate("Time")[0][:, 0]
        tcorr = t-t.mean()

        #plotting
        fs = 15
        plt.figure()
        plt.plot(tcorr, corr, "bo-")
        plt.xlabel(r'$\tau [s]$', fontsize=fs)
        plt.title(
            self.expe_id+", Auto-correlation, $\lambda = ["+str(lstart)+", "+str(lstop)+"]$ nm")
        plt.grid()

    def tshif(self, roi, t_start, t_stop, wstart, wstop, N, background_interval=[0]):
        """
        It averages the spectra in the given wavelength interval, then it slices and
        averages it along the time axis in the intervals when the beam was on and off,
        and it subtracts the two. This process is done multiple times while the temporal
        shift between the spectral data time axis and the beam interval times is varying.
        Finally, it plots the average spectral intensity / frame in the given intervals
        in the subtracted spectra vs the temporalshifts.

        Parameters
        ----------
        roi : Region Of Interest (ROI), in other words spectral channel (int)
        t_start : beginning of the time interval (for APDCAM and CXRS data both) (float)
        t_stop : end of the time interval (for APDCAM and CXRS data both) (float)
        wstart : Lower end of the wavelength interval of the spectral line (float)
        wstop : Higher end of the wavelength interval of the spectral line (float)
        N : number of iterations for the temporal shift varying
        background_interval : A spectral interval that does not contain any notable
        spectra lines (list with two elements). The default is [0].

        """
        if(self.campaign == "OP2.1"):
            tsh = np.linspace(-0.075, 0.075, N) #temporal shifts
            lineint = np.zeros((N))
            for i in range(len(tsh)):
                
                #getting the beam intervals
                d_beam_on = flap.get_data('W7X_ABES', exp_id=self.expe_id, name='Chopper_time',
                                          options={
                                              'State': {'Chop': 0, 'Defl': 0}},
                                          object_name='Beam_on',
                                          coordinates={'Time': self.APDCAM_timerange})

                d_beam_off = flap.get_data('W7X_ABES', exp_id=self.expe_id, name='Chopper_time',
                                           options={
                                               'State': {'Chop': 1, 'Defl': 0}},
                                           object_name='Beam_off',
                                           coordinates={'Time': self.APDCAM_timerange})

                c = d_beam_on.get_coordinate_object("Time") #shifting them
                c.start = c.start + tsh[i]
                c = d_beam_off.get_coordinate_object("Time")
                c.start = c.start + tsh[i]
                
                #slicing
                ROI1 = self.dataobj.slice_data(
                    slicing={"ROI": "P0"+str(roi), "Time": flap.Intervals(t_start, t_stop)})
                if(background_interval != [0]):
                    ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                "Wavelength": flap.Intervals(background_interval[0],
                                background_interval[1])},summing={"Wavelength": "Mean",
                                                                  "Time": "Mean"})
                    ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data

                s_on = ROI1.slice_data(slicing={'Time': d_beam_on}, summing={
                                       "Rel. Time in int(Time)": "Mean"})
                s_off = ROI1.slice_data(slicing={'Time': d_beam_off}, summing={
                                        "Rel. Time in int(Time)": "Mean"})
                s_on_intervals_full = ROI1.slice_data(
                    slicing={'Time': d_beam_on})
                s_off_intervals_full = ROI1.slice_data(
                    slicing={'Time': d_beam_off})
                
                #averaging over the intervals (first elements are taken out)
                s_on_intervals = s_on_intervals_full.data[:, 1:, ].mean(axis=1)
                s_off_intervals = s_off_intervals_full.data[:, 1:, ].mean(
                    axis=1)
                s_on_data = s_on_intervals.mean(axis=1)
                s_off_data = s_off_intervals.mean(axis=1)
                s_on.data = s_on_data
                s_off.data = s_off_data

                s_subs = s_on
                s_subs.data = s_on.data-s_off.data
                lineint[i]=s_subs.slice_data(slicing={"Wavelength":flap.Intervals(wstart,wstop)},
                                               summing={"Wavelength": "Mean"}).data

            fs = 15 #plotting
            plt.figure()
            plt.plot(tsh, lineint, marker="o")
            plt.xlabel("$t_{shift}$ [s]", fontsize=fs)
            plt.ylabel("Spectral intensity", fontsize=fs)
            plt.title("["+str(wstart)+", "+str(wstop)+"] nm, [" +
                      str(t_start)+", "+str(t_stop)+"] s", fontsize=fs)

            print("The best temporal shift length in second:")
            print(tsh[np.argmax(lineint)])
        else:
            raise ValueError(
                "For that campaign this method is not implemented yet.")

    def active_passive(self, roi, t_start, t_stop, background_interval=[0],
                       error=False, plotting=True, plotting_wavelength = []):
        """
        It slices and averages it along the time axis in the intervals when the beam
        was on and off, and it subtracts the two.

        Parameters
        ----------
        roi : Region Of Interest (ROI), in other words spectral channel (int)
        t_start : beginning of the time interval (for APDCAM and CXRS data both) (float)
        t_stop : end of the time interval (for APDCAM and CXRS data both) (float)
        background_interval : A spectral interval that does not contain any notable
        spectra lines (list with two elements). The default is [0], which means no
        interval is given.
        error : wether one would like to perform error calculation or not (Boolean).
        The default is False.
        plotting : wether it should plot the results or not (Boolean or string). 
                   It can also be "on-off", when it plots the spectra measured and 
                   averaged during the beam-on, beam-off intervals, or it can be
                   "active", when it plots the difference of the two. The default is True.
        plotting_wavelength: a wavelength interval in which the spectra will be plotted
        (list of two floats). The default is an empty list.

        Returns
        -------
        s_subs : subtracted spectrum

        """
        if(self.campaign == "OP2.1"):

            # slicing the data
            ROI1 = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                                    "Time": flap.Intervals(t_start, t_stop)})
            error_on = None
            error_off = None
            if(error == True):
                s_on_sliced = ROI1.slice_data(
                    slicing={'Time': self.d_beam_on})  # for error calculation
                s_off_sliced = ROI1.slice_data(
                    slicing={'Time': self.d_beam_off})

                # the intervals are taken as independent measurements
                error_on = spectral_error_calc_op21(s_on_sliced)
                error_off = spectral_error_calc_op21(s_off_sliced)

            if(background_interval != [0]):
                ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                            "Wavelength": flap.Intervals(background_interval[0],
                             background_interval[1])},
                             summing={"Wavelength": "Mean", "Time": "Mean"})
                ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data

            #averaging over the intervals (first elements are taken out)
            s_on_intervals_full = ROI1.slice_data(
                slicing={'Time': self.d_beam_on})
            s_off_intervals_full = ROI1.slice_data(
                slicing={'Time': self.d_beam_off})
            s_on_intervals = s_on_intervals_full.data[:, 1:, ].mean(axis=1)
            s_off_intervals = s_off_intervals_full.data[:, 1:, ].mean(axis=1)
            s_on_data = s_on_intervals.mean(axis=1)
            s_off_data = s_off_intervals.mean(axis=1)
            s_on = ROI1.slice_data(slicing={'Time': self.d_beam_on},
                                   summing={"Rel. Time in int(Time)": "Mean"})
            s_off = ROI1.slice_data(slicing={'Time': self.d_beam_off},
                                    summing={"Rel. Time in int(Time)": "Mean"})
            s_on.data = s_on_data
            s_off.data = s_off_data

            if(plotting == True or plotting == "on-off"):
                
                if(plotting_wavelength != []):
                    s_on=s_on.slice_data(slicing={'Wavelength':flap.Intervals(
                        plotting_wavelength[0],plotting_wavelength[1])})
                    s_off=s_off.slice_data(slicing={'Wavelength':flap.Intervals(
                        plotting_wavelength[0],plotting_wavelength[1])})
                
                plt.figure()
                s_on.plot(axes="Wavelength")
                s_off.plot(axes="Wavelength")
                legend = ["beam on", "beam off"]
                plt.legend(legend,fontsize = 15)
                plt.title(self.expe_id+", beam on and off spectra, ROI = "+str(roi)+
                          ", ["+str(t_start)+","+str(t_stop)+"] s",fontsize=15)
                plt.ylabel("Intensity [a.u.]")
                plt.grid()

            s_subs = s_on
            s_subs.data = s_on.data-s_off.data  # gaining the active spectrum
            if(error == True):
                s_subs.error = np.sqrt(error_on**2 + error_off**2)

            if(plotting == True or plotting == "active"):
                if(plotting_wavelength != []):
                    s_subs=s_subs.slice_data(slicing={'Wavelength':flap.Intervals(
                        plotting_wavelength[0],plotting_wavelength[1])})
                plt.figure()
                s_subs.plot(axes="Wavelength")
                plt.title(self.expe_id+", active spectrum, ROI = P0"+str(roi)+
                          ", ["+str(t_start)+","+str(t_stop)+"] s",fontsize=15)
                plt.ylabel("Intensity [a.u.]")
                plt.grid()

            return s_subs
        else:
            raise ValueError(
                "For that campaign this method is not implemented yet.")

    def get_line_intensity(self, roi, t_start, t_stop, lstart, lstop,
                           background_interval=[0], plotting=False):
        """
        It slices and averages it along the time axis in the intervals when the beam
        was on and off, and it subtracts the two. Finally, it calculates the line
        intensity in the given wavelength interval for both by fitting a Gaussian.

        Parameters
        ----------
        roi : Region Of Interest (ROI), in other words spectral channel (int)
        t_start : beginning of the time interval (for APDCAM and CXRS data both) (float)
        t_stop : end of the time interval (for APDCAM and CXRS data both) (float)
        lstart : Lower end of the wavelength interval of the spectral line (float)
        lstop : higher end of the wavelength interval of the spectral line (float)
        background_interval : A spectral interval that does not contain any notable
        spectra lines (list with two elements). The default is [0], which means no
        interval is given.
        plotting : wether it should plot the results or not (Boolean). The default is True.

        Returns
        -------
        Average maximum line intensity at beam on and off

        """
        if(self.campaign == "OP2.1"):
            # slicing the data
            ROI1 = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                        "Time": flap.Intervals(t_start, t_stop),
                                        "Wavelength": flap.Intervals(lstart,lstop)})

            if(background_interval != [0]): #continuous background subtraction
                ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                            "Wavelength": flap.Intervals(background_interval[0],
                             background_interval[1])},
                             summing={"Wavelength": "Mean", "Time": "Mean"})
                ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data
            s_on_intervals_full = ROI1.slice_data(
                slicing={'Time': self.d_beam_on})
            s_off_intervals_full = ROI1.slice_data(
                slicing={'Time': self.d_beam_off})
            s_on_intervals = s_on_intervals_full.data[:, 1:, ].mean(axis=1)
            s_off_intervals = s_off_intervals_full.data[:, 1:, ].mean(axis=1)
            s_on_data = s_on_intervals.mean(axis=1)
            s_off_data = s_off_intervals.mean(axis=1)
            s_on = ROI1.slice_data(slicing={'Time': self.d_beam_on},
                                   summing={"Rel. Time in int(Time)": "Mean"})
            s_off = ROI1.slice_data(slicing={'Time': self.d_beam_off},
                                    summing={"Rel. Time in int(Time)": "Mean"})
            s_on.data = s_on_data
            s_off.data = s_off_data

            #fitting Gaussians to the lines
            lambd = s_on.coordinate("Wavelength")[0]
            gaus = lambda x, A, s, mu: A*np.e**(-(((x-mu)**2)/s**2))

            popton, pcovon = curve_fit(gaus, lambd, s_on.data, p0=[
                                       max(s_on.data), 0.1, lambd.mean()])
            print(popton[0])
            poptoff, pcovoff = curve_fit(gaus, lambd, s_off.data, p0=[
                                         max(s_off.data), 0.1, lambd.mean()])
            if(plotting == True): #plotting the results
                fs = 15
                plt.figure()
                plt.plot(lambd, s_on.data, "+")
                plt.plot(lambd, gaus(lambd, *popton))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")

                plt.figure()
                plt.plot(lambd, s_off.data, "+")
                plt.plot(lambd, gaus(lambd, *poptoff))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam off line intensity fit")

            return popton[0], poptoff[0]
        else:
            raise ValueError(
                "For that campaign this method is not implemented yet.")

    def error_distr(self, t_start, t_stop, lstart, lstop,
                    background_interval=[0], plotting=False, ROI="ALL"):
        """
        Error(Intensity) distribution calculator function. It assumes that the
        fluctuations in spectra which affect more than one point are not errors. It
        also fits a root square function to the results, assuming only electronic and
        photonic error.
        
        Parameters
        ----------
        t_start : beginning of the time interval (float)
        t_stop : end of the time interval (float)
        lstart : lower end of the wavelength interval of the spectral line (float)
        lstop : higher end of the wavelength interval of the spectral line (float)
        background_interval : A spectral interval that does not contain any notable
        spectra lines (list with two elements). The default is [0], which means no
        interval is given.
        plotting : wether it should plot the results or not (Boolean). The default is False.
        ROI : Region Of Interest, in other words spectral channel (int or string).
        The default is "ALL". Besides that, it can be 1,2,3 or 4.
        
        Returns
        -------
        Fitted average coefficients (a,b) of the E(I) = a*sqrt(I)+b function

        """
        if(self.campaign == "OP2.1" and ROI == "ALL"):
            errorparam = np.zeros((4, 2))
            for roi in range(1, 5):

                # slicing the data
                ROI1 = self.dataobj.slice_data(
                    slicing={"ROI": "P0"+str(roi), "Time": flap.Intervals(t_start, t_stop)})
                ROI1 = ROI1.slice_data(
                    slicing={"Wavelength": flap.Intervals(lstart, lstop)})
                s_on_sliced = ROI1.slice_data(
                    slicing={'Time': self.d_beam_on})  # for error calculation
                s_off_sliced = ROI1.slice_data(
                    slicing={'Time': self.d_beam_off})

                # the intervals are taken as independent measurements
                error_on = list(indep_spectral_error_calc_op21(s_on_sliced))
                error_off = list(indep_spectral_error_calc_op21(s_off_sliced))
                
                if(background_interval != [0]): #constant background removal
                    ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                "Wavelength": flap.Intervals(background_interval[0],
                                 background_interval[1])},summing={"Wavelength":"Mean",
                                 "Time":"Mean"})
                    ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data
                    
                #averaging
                s_on = ROI1.slice_data(slicing={'Time': self.d_beam_on}, summing={
                                       "Rel. Time in int(Time)": "Mean"})
                s_off = ROI1.slice_data(slicing={'Time': self.d_beam_off}, summing={
                                        "Rel. Time in int(Time)": "Mean"})
                
                #adding the errors and creating the distribution
                int_on = list(s_on.data)
                int_off = list(s_off.data)
                err = np.array(error_on + error_off)
                intensity = np.array(int_on + int_off)
                err = err[intensity.argsort()]
                intensity = intensity[intensity.argsort()]

                #fitting the funcion
                sq = lambda x, a, b: a*np.sqrt(x) + b
                if(plotting == True):
                    fs = 15
                    plt.figure()
                    plt.plot(intensity, err, "+")
                    popt, pcov = curve_fit(sq, intensity, err, p0=[0.03, 0.3])
                    plt.plot(intensity, sq(intensity, *popt))
                    plt.xlabel("Spectral intensity", fontsize=fs)
                    plt.ylabel("Error", fontsize=fs)
                    plt.legend(["Data", "Fit"], fontsize=fs-2)
                    plt.grid()
                    plt.title(self.expe_id+", ROI"+str(roi) +
                              ", t = ["+str(t_start)+","+str(t_stop)+"] s", fontsize=fs)
                errorparam[roi-1, 0] = popt[0]
                errorparam[roi-1, 1] = popt[1]
            print(errorparam)
            print(errorparam.mean(axis=0))
            return errorparam.mean(axis=0)

        elif(self.campaign == "OP2.1"):
            # slicing the data
            ROI1 = self.dataobj.slice_data(
                slicing={"ROI": "P0"+str(ROI), "Time": flap.Intervals(t_start, t_stop)})
            ROI1 = ROI1.slice_data(
                slicing={"Wavelength": flap.Intervals(lstart, lstop)})
            s_on_sliced = ROI1.slice_data(
                slicing={'Time': self.d_beam_on})  # for error calculation
            s_off_sliced = ROI1.slice_data(slicing={'Time': self.d_beam_off})

            # the intervals are taken as independent measurements
            error_on = list(indep_spectral_error_calc_op21(s_on_sliced))
            error_off = list(indep_spectral_error_calc_op21(s_off_sliced))
            
            if(background_interval != [0]): #constant background removal
                ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(ROI),
                             "Wavelength": flap.Intervals(background_interval[0],
                             background_interval[1])},summing={"Wavelength": "Mean",
                             "Time": "Mean"})
                ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data
            
            #averaging
            s_on = ROI1.slice_data(slicing={'Time': self.d_beam_on}, summing={
                                   "Rel. Time in int(Time)": "Mean"})
            s_off = ROI1.slice_data(slicing={'Time': self.d_beam_off}, summing={
                                    "Rel. Time in int(Time)": "Mean"})
            
            #adding the errors and creating the distribution
            int_on = list(s_on.data)
            int_off = list(s_off.data)
            err = np.array(error_on + error_off)
            intensity = np.array(int_on + int_off)
            err = err[intensity.argsort()]
            intensity = intensity[intensity.argsort()]
            
            #fitting the funcion
            sq = lambda x, a, b: a*np.sqrt(x) + b
            if(plotting == True):
                fs = 15
                plt.figure()
                plt.plot(intensity, err, "+")
                popt, pcov = curve_fit(sq, intensity, err, p0=[0.03, 0.3])
                plt.plot(intensity, sq(intensity, *popt))
                plt.xlabel("Spectral intensity", fontsize=fs)
                plt.ylabel("Error", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.grid()
                plt.title(self.expe_id+", ROI"+str(ROI) +
                          ", t = ["+str(t_start)+","+str(t_stop)+"] s", fontsize=fs)
            print(popt[0])
            print(popt[1])
            return popt[0], popt[1]
        else:
            raise ValueError(
                "For that campaign this method is not implemented yet.")

    def C_line_generator(self, mu, kbt, A,sim=False):
        """
        A function to generate theoretical carbon ion line shapes based on the Zeeman
        components of the line which are defined elsewhere.

        Parameters
        ----------
        mu : term that constitutes from the Doppler shift of the line, the
        wavelength uncertainty of the spectrometer (float), and the line central wavelength
        kbt : temperature times Boltzmann constant (float)
        A : Line intensity factor (float)
        sim : (Boolean) Wether the function is used for simulating experiments, or fitting
        measurements. The code reads the instrumental functions from different places
        in these two cases, and generates different wavelength grids. The default is
        False.

        Returns
        -------
        complete_spectrum : Theoretical spectrum in case of the given parameters

        """
        wl_grid = None
        if(sim==False):
            wl_values = wavelength_grid_generator_op21(
                self.grid, self.wavelength_setting, self.current_roi,self.supl_data_path)
            # slicing the wavelength grid
            wl_grid0 = wl_values[wl_values > self.wstart]
            wl_grid = wl_grid0[self.wstop > wl_grid0]
        elif(sim==True):
            wl_values = wavelength_grid_generator_op21(
                self.simgrid, self.wavelength_setting, self.current_roi,self.supl_data_path)
            # slicing the wavelength grid
            wl_grid0 = wl_values[wl_values > self.wstart]
            wl_grid = wl_grid0[self.wstop > wl_grid0] 
        # Doppler shift + calibration uncertainty
        locations = self.zc_locations
        projection = np.zeros((wl_grid.shape[0]))
        for i in range(len(self.zc_intensities)):
            diff = abs(wl_grid - locations[i])
            sor = np.argsort(diff)
            closest_ind = sor[:2]
            distance = diff[sor[:2]]
            I2 = self.zc_intensities[i] / (1 + distance[1]/distance[0])
            projection[closest_ind[1]] += I2
            projection[closest_ind[0]] += self.zc_intensities[i] - I2

        # addition of the Doppler broadening
        doppler_broadening=lambda kbt,a:a*np.sqrt(2*1.602176487*abs(kbt)/19.9419)/30000
        s = doppler_broadening(kbt, mu)
        gauss=lambda lambd,s,A,mu:A*np.e**(-(((lambd-mu)**2)/s**2))
        gaussian = gauss(wl_grid, s, A, mu)
        doppler_spectrum = np.convolve(projection, gaussian, mode="same")
        instr = None
        if(sim==False):
            # convolution with instrumental function
            instr = np.load(self.instr_funcs_datapath+self.grid+"_P0"+
                            str(self.current_roi) +
                            "_"+str(int(self.dslit))+"micron_slit.npy").ravel()
        elif(sim==True):
            # convolution with instrumental function
            if(self.simd == 100):
                instr = np.load(self.instr_funcs_datapath+self.simgrid+"_P0"+
                                str(self.current_roi) +
                            "_"+str(int(self.simd))+"micron_slit.npy").ravel()
                print(instr)
                plt.figure()
                plt.plot(instr)
                raise ValueError("stop")
            else:
                instr = np.load(self.instr_funcs_datapath+self.simgrid+"_P03" +
                            "_"+str(int(self.simd))+"micron_slit.npy").ravel()
                print(instr)
                plt.figure()
                plt.plot(instr)
                raise ValueError("stop")
        complete_spectrum = np.convolve(doppler_spectrum, instr, mode="same")

        return complete_spectrum

    def C_fitfunc(self, esti):
        """
        It computes the weighted, squared, summed difference between the measured
        and the generated spectral line.

        Parameters
        ----------
        esti : (np.array) a vector that contains relevant parameters in the following 
        way: np.array([mu, kbt, A]) - see the variables above

        Returns
        -------
        chi^2 (normalized by # of parameters, data points)
        """
        mu = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu, kbt, A)
        C = (modelled - self.active.data)/self.active.error
        return (np.dot(C, C) - 3) / self.active.data.shape[0]
    
    def C_fitfunc_sim(self, esti):
        """
        It computes the weighted, squared, summed difference between the simulated
        and the generated spectral line.

        Parameters
        ----------
        esti : (np.array) a vector that contains relevant parameters in the following 
        way: np.array([mu, kbt, A]) - see the variables above

        Returns
        -------
        chi^2 (normalized by # of parameters, data points)
        """
        mu = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu, kbt, A,sim = True)
        C = (modelled - self.simulated)/self.simulated_error
        return (np.dot(C, C) - 3) / self.simulated.shape[0]

    def C_fitfunc_plot(self, esti):
        """
        It computes the weighted, squared, summed difference between the measured
        and the generated spectral line. It also plots the measured and the generated
        lines.

        Parameters
        ----------
        esti : (np.array) a vector that contains relevant parameters in the following 
        way: np.array([mu, kbt, A]) - see the variables above

        Returns
        -------
        chi^2 (normalized by # of parameters, data points)
        """
        mu = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu, kbt, A)
        C = (modelled - self.active.data)/self.active.error

        lambd = self.active.coordinate("Wavelength")[0]

        fs = 15
        plt.figure()
        plt.errorbar(lambd, self.active.data, self.active.error, color="blue")
        plt.plot(lambd, modelled, color="red")
        plt.xlabel("Wavelength [nm]", fontsize=fs)
        plt.ylabel("Intensity [a.u.]", fontsize=fs)
        plt.grid()
        plt.legend(["Calculated", "Measured"], loc="best", fontsize=(fs-2))
        return (np.dot(C, C) - 3) / self.active.data.shape[0]
    
    def second_deriv(self,sol, fmin, h, i):
        """
        It is part of an error calculation algorithm that is based on the Hessian matrix
        of chi^2. However, this algorithm gave unreasonably large errors, therefore it
        is not in use
        """
        sol_p = sol.copy()
        sol_n = sol.copy()
        sol_p[i] = sol_p[i] + h
        sol_n[i] = sol_n[i] - h
        fp1 = self.C_fitfunc(sol_p)
        fm1 = self.C_fitfunc(sol_n)
        return (fp1 + fm1 - 2*fmin)/h**2

    def partial_cross_deriv(self,sol, fmin, i, j, hi, hj):
        """
        It is part of an error calculation algorithm that is based on the Hessian matrix
        of chi^2. However, this algorithm gave unreasonably large errors, therefore it
        is not in use
        """
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

        fpp = self.C_fitfunc(sol_pp)
        fpn = self.C_fitfunc(sol_pn)
        fnp = self.C_fitfunc(sol_np)
        fnn = self.C_fitfunc(sol_nn)

        return (fpp-fpn-fnp+fnn) / (4*hi*hj)

    def hesse(self,sol, fmin, h):
        """
        It is part of an error calculation algorithm that is based on the Hessian matrix
        of chi^2. However, this algorithm gave unreasonably large errors, therefore it
        is not in use
        """
        hess = np.matrix(np.zeros((3, 3)))
        for i in range(3):
            for j in range(3):
                if(i == j):
                    hess[i, i] = self.second_deriv(sol, fmin, h[i], i)
                else:
                    hess[i, j] = self.partial_cross_deriv(sol, fmin, i, j, h[i], h[j])
        return hess

    def error_from_hesse(self,sol,fmin, h):
        """
        It is part of an error calculation algorithm that is based on the Hessian matrix
        of chi^2. However, this algorithm gave unreasonably large errors, therefore it
        is not in use
        """
        alpha = self.hesse(sol, fmin, h)
        I = np.linalg.inv(alpha)
        return np.sqrt(abs(I))
    
    def C_fitfunc_plot_sim(self, esti):
        """
        It computes the weighted, squared, summed difference between the simulated
        and the generated spectral line. It also plots the simulated and the generated
        lines.

        Parameters
        ----------
        esti : (np.array) a vector that contains relevant parameters in the following 
        way: np.array([mu, kbt, A]) - see the variables above

        Returns
        -------
        chi^2 (normalized by # of parameters, data points)
        """
        mu = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu, kbt, A,sim = True)
        C = (modelled - self.simulated)/self.simulated_error

        # loading the wavelength grid
        wl_values = wavelength_grid_generator_op21(
            self.simgrid, self.wavelength_setting, self.current_roi,self.supl_data_path)
        wl_grid0 = wl_values[wl_values > self.wstart]  # slicing the wavelength grid
        lambd = wl_grid0[self.wstop > wl_grid0]

        fs = 15
        plt.figure()
        plt.errorbar(lambd, self.simulated, self.simulated_error, color="blue")
        plt.plot(lambd, modelled, color="red")
        plt.xlabel("Wavelength [nm]", fontsize=fs)
        plt.ylabel("Intensity [a.u.]", fontsize=fs)
        plt.grid()
        plt.legend(["Calculated", "Measured"], loc="best", fontsize=(fs-2))
        return (np.dot(C, C) - 3) / self.simulated.shape[0]
    
    
    def C_line_simulator(self,mu,kbt,A,fittype,scalef=None,plots=False,sim=False):
        """
        A function that generates a realistic Carbon ion line spectrum including 
        noise. 

        Parameters
        ----------
        mu : term that constitutes from the Doppler shift of the line, the
        wavelength uncertainty of the spectrometer (float), and the line central wavelength
        kbt : temperature times Boltzmann constant (float)
        A : Line intensity factor (float)
        fittype : (string) The kind of line that is to be generated. At the moment, 
        it can be "CV","CVI", or None. The default is None.
        scalef : scaling number between the measured and to be generated line intensity
        (float). Only relevant when sim=True. The default is None.
        plots : (Boolean) wether plot the generated spectrum (during the generation
        stages and afterwards) or not. The default is False.
        sim : (Boolean) wether the function is used for simulating experiments, or fitting
        measurements. The code reads the instrumental functions from different places
        in these two cases, and generates different wavelength grids. The default is
        False.

        Returns
        -------
        calculated : generated realistic spectrum (np.array)
        err : uncertainties of the spectrum (np.array)

        """
        gaus = lambda x, A, s, mu: A*np.e**(-(((x-mu)**2)/s**2))
        if(sim==False):
            #slicing to the wavelength interval
            measured = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(self.wstart, self.wstop)})
            
            #calculating the measured line intensity maximum by fitting a Gaussian
            lamb = measured.coordinate("Wavelength")[0]
            popt, pcov = curve_fit(gaus, lamb, measured.data, p0=[
                                   max(measured.data), 0.1, lamb.mean()])
            if(plots == True):
                #loading the wavelength grid
                wl_values = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, self.current_roi)
                wl_grid0=wl_values[wl_values>self.wstart] #slicing the wavelength grid
                lambd = wl_grid0[self.wstop > wl_grid0]
                fs = 15
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, gaus(lambd, *popt))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")
    
            #generating a line with the same maximum intensity as the measured
            calculated = self.C_line_generator(mu, kbt, A)
            calculated = popt[0]*calculated/max(calculated)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            #the error of the generated spectra is the same as for the measured
            err = measured.error
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "o")
                plt.errorbar(lambd, calculated, err)
                plt.grid()
    
            #adding noise to the generated spectra based on its errors
            calculated = np.random.normal(loc=calculated, scale=err)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            return calculated, err
        
        if(sim==True):
            #slicing to the wavelength interval
            measured = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(self.wstart, self.wstop)})
            
            #calculating the measured line intensity maximum by fitting a Gaussian
            lamb = measured.coordinate("Wavelength")[0] 
            popt, pcov = curve_fit(gaus, lamb, measured.data, p0=[
                                   max(measured.data), 0.1, lamb.mean()])
            if(plots == True):
                #loading the wavelength grid
                wl_values = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, self.current_roi)
                wl_grid0=wl_values[wl_values>self.wstart] #slicing the wavelength grid
                lambd = wl_grid0[self.wstop > wl_grid0]
                fs = 15
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, gaus(lambd, *popt))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")
    
            calculated = self.C_line_generator(mu, kbt, A,sim=True)
            
            gridfac = None 
            #grid and slit width factor for cases when one would like to generate 
            #a spectrum by different grid or slit (or both) than the measured
            #it depends on the wavelength
            if(fittype == "CVI"):
                baseint = 11230
                if(self.simgrid == "1200g_per_mm" and self.simd == 100):
                    gridfac = 1
                elif(self.simgrid == "1800g_per_mm" and self.simd == 100):
                    gridfac = 16922/baseint
                elif(self.simgrid == "2400g_per_mm" and self.simd == 100):
                    gridfac = 10688/baseint
                elif(self.simgrid == "1200g_per_mm" and self.simd == 70):
                    gridfac = 12810/(baseint*1.5)
                elif(self.simgrid == "1800g_per_mm" and self.simd == 70):
                    gridfac = 20885/(baseint*1.5)
                elif(self.simgrid == "2400g_per_mm" and self.simd == 70):
                    gridfac = 13971/(baseint*1.5)
                elif(self.simgrid == "1200g_per_mm" and self.simd == 50):
                    gridfac = 14715/(baseint*2)
                elif(self.simgrid == "1800g_per_mm" and self.simd == 50):
                    gridfac = 21279/(baseint*2)
                elif(self.simgrid == "2400g_per_mm" and self.simd == 50):
                    gridfac = 14274/(baseint*2)
                elif(self.simgrid == "1200g_per_mm" and self.simd == 35):
                    gridfac = 12470/(baseint*2.8)
                elif(self.simgrid == "1800g_per_mm" and self.simd == 35):
                    gridfac = 21396/(baseint*2.8)
                elif(self.simgrid == "2400g_per_mm" and self.simd == 35):
                    gridfac = 16699/(baseint*2.8)
                else:
                    raise ValueError("Wrong grid or slit size.")
                    
            elif(fittype == "CV"):
                baseint = 3840
                if(self.simgrid == "1200g_per_mm"):
                    gridfac = 1
                elif(self.simgrid == "1800g_per_mm"):
                    gridfac = 2707/baseint
                elif(self.simgrid == "2400g_per_mm"):
                    gridfac = 1505/baseint
                    
                slitfac = None
                if(self.simd == 100):
                    slitfac = 1
                elif(self.simd == 70):
                    slitfac = 0.81823415
                elif(self.simd == 50):
                    slitfac = 0.65055357
                elif(self.simd == 35):
                    slitfac = 0.46871601   
                gridfac = gridfac*slitfac
            
            #generating the line spectrum with the given parameters
            calculated = gridfac*scalef*popt[0]*calculated/max(calculated)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
                
            sq = lambda x, a, b: a*np.sqrt(x) + b
            #assuming 10% modulation
            err = sq(calculated*10, self.errparam[0], self.errparam[1])

            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "o")
                plt.errorbar(lambd, calculated, err)
                plt.grid()

            #adding noise
            calculated = np.random.normal(loc=calculated, scale=err)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            return calculated, err
    
    def C_Ti_error_sim_me(self,mu,kbt,A,iter_num,fittype,plots=False):
        """
        Function for calculating uncertainty of kb*T_i for a spectra fitted on measured
        data. It uses Monte Carlo error calculation.

        Parameters
        ----------
        mu : a term that constitutes from the Doppler shift of the line, the
        wavelength uncertainty of the spectrometer (float), and the line central wavelength
        kbt : temperature times Boltzmann constant (float)
        A : Line intensity factor (float)
        iter_num : Number of iterations for the Monte Carlo error calculation process
        (int)
        fittype : (string) The kind of line that is to be generated. At the moment, 
        it can be "CV" or "CVI".
        plots : (Boolean) wether plot the generated spectra during the error calculation
        process or not. The default is False.

        Returns
        -------
        float
            Standard deviation of the fitted temperatures to the generated spectra.
        T_i : np.array with size iter_num
            Fitted temperatures to the generated spectra.
        chisq : np.array with size iter_num
            Chi^2 from every iteration.

        """
        met = "Powell" #as far as I know that is the best choice for such purpose
        line_param = np.array([mu, kbt, A])

        T_i = np.zeros((iter_num))
        chisq = np.zeros((iter_num))

        for i in range(iter_num):
            print("Iteration "+str(i))
            #simulating same spectrum as the measured based on its maximum 
            #intensity and its errors
            self.simulated, self.simulated_error = self.C_line_simulator(mu, kbt, A,
                                                                    fittype, plots=False)
            if(plots == True): #plot with initial parameters before fit
                es_chisq = self.C_fitfunc_plot_sim(line_param)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
            
            #minimizing the chi^2 using scipy.optimize.minimize
            solution = minimize(self.C_fitfunc_sim, line_param, method=met,
                                bounds=((None,None),(0.1,None),(None,None)),tol=1e-8,
                                options={"maxiter": 2000})
            if(solution.success == False):
                raise ValueError("Failed T_i fit")
            sol = solution.x
            print(sol[1])
            T_i[i] = sol[1]
            chisq[i] = solution.fun
            if(plots == True): #plotting the line shape with the fitted parameters
                self.C_fitfunc_plot_sim(sol)
                R_plot = round(self.dataobj.coordinate("Device R")[0][0,
                                                        (self.current_roi-1), 0], 4)
                plt.title("R = "+str(R_plot)+" m, $\chi^2 = $" +
                          str(round(solution.fun, 6))+", $T_C$ = "+str(round(sol[1], 2)))

        return np.std(T_i), T_i, chisq

    def tempfit(self,fittype,roi,wstart,wstop,mu,kbt,A,dslit,t_start,t_stop,bcg,N,
                plots=False,options = None):
        """
        A function for fitting the ion temperature (among A, mu) by modelling the
        line shape. The uncertainties are gained by Monte Carlo error calculation.
        The function prints the results, and also plots the fitted line upon request.

        Parameters
        ----------
        fittype : (string) The kind of line that is to be generated. At the moment, 
        it can be "CV" or "CVI".
        roi : Region Of Interest, in other words spectral channel (int).
        wstart : lower end of the wavelength interval of the spectral line (float)
        wstop : higher end of the wavelength interval of the spectral line (float)
        mu : initial guess for a term that constitutes from the Doppler shift of the line, the
        wavelength uncertainty of the spectrometer (float), and the line central wavelength
        kbt : initial guess for temperature times Boltzmann constant (float)
        A : initial guess for line intensity factor (float)
        dslit : slit width in the spectrometer in the discharge (int)
        t_start : lower end of the interval for which the ion temperature to be 
        calculated
        t_stop : higher end of the interval for which the ion temperature to be 
        calculated
        bcg : A wavelength interval that does not contain any notable
        spectral lines (list with two floats).
        N : Number of iterations for the Monte Carlo error calculation process (int)
        plots : (Boolean) wether plot the generated spectra during the error calculation
        process or not. The default is False.
        options : dictionary
            'Supplementary_data_path': string, path of the supplementary data
            (like wavelength calibration)
            'Instrumental_functions_datapath': string, path of the saved instrumental
            functions, which should be .npy files
        """
        default_options = {"Supplementary_data_path":"data"}
        _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES_CXRS')
        try:
            self.supl_data_path = _options['Supplementary_data_path']
            self.instr_funcs_datapath = _options['Instrumental_functions_datapath']
        except (KeyError, TypeError):
            self.supl_data_path = 'data'
            self.instr_funcs_datapath = 'data'
        if(self.campaign == "OP2.1"):
            
            #central position of the first lens
            centre_of_lens = np.array([1.305624, 6.094843, -3.013095])
            #viewed positions of the channels
            roi_pos = np.loadtxt(self.supl_data_path+"ROI_pos.txt")
            self.observation_directions = np.zeros(
                (roi_pos.shape[0], roi_pos.shape[1]))
            for i in range(roi_pos.shape[1]):
                self.observation_directions[:,i] = roi_pos[:, i]-centre_of_lens
            self.current_roi = roi
            self.dslit = dslit
            self.wstart = wstart
            self.wstop = wstop
            self.simgrid = self.grid
            self.simd = dslit

            #getting the measured spectra
            self.active = self.active_passive(roi, t_start, t_stop,
                                    background_interval=bcg, error=True, plotting=False)
            self.active = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(wstart, wstop)})
            if(fittype == "CV"):
                #magnetic field data from Gábor Cseh
                self.B = np.loadtxt(self.supl_data_path+self.expe_id+".txt")
                self.Babs = np.sqrt(np.sum(self.B[:, roi-1]**2))
                divider = np.sqrt(
                    np.sum(self.observation_directions[:, roi-1]**2))*self.Babs
                n = self.observation_directions[:, roi-1]
                self.current_theta = np.arccos(
                    np.sum(n*self.B[:, roi-1])/divider)*180/np.pi
                
                #no exact data - saved approximation proposed by Oliver Ford
                #for the Zeeman components, done by the webservice of Dorothea Gradic
                with open(self.supl_data_path+'getZeeman_CV_ROI'+str(self.current_roi)+'.json', 'r') as f:
                    datafile = json.load(f)
                
                self.zc_locations = np.array(datafile['wavelengths'])/10
                self.zc_intensities = np.array(datafile['amplitude'])
                
                #fitting the ion temperature, and other parameters
                esti = np.array([mu, kbt, A])
                es_chisq = self.C_fitfunc_plot(esti)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
                solution = minimize(self.C_fitfunc, esti, method="Powell",
                                    bounds=((None, None), (0.1, None), (None, None)), tol=1e-12,
                                    options={"maxiter": 2000})

                print(solution)
                if(solution.success == True):
                    #error calculation
                    sol = solution.x
                    Terr, T_iters, chi_iters = self.C_Ti_error_sim_me(sol[0],
                                        sol[1],sol[2],N,fittype,plots=plots)
                    # h = np.array([1e-5, 1e-3, 1e-8])
                    # hessian_error = self.error_from_hesse(sol, solution.fun, h)[1,1]
                    # print("Error based on Hessian matrix (in eV):")
                    # print(hessian_error)
                    #plotting
                    R_plot = round(self.dataobj.coordinate(
                        "Device R")[0][0, (roi-1), 0], 4)
                    self.C_fitfunc_plot(sol)
                    plt.title(self.expe_id+", channel "+str(roi)+
                              ", R = "+str(round(R_plot,4))+" m, $\chi^2 = $" +
                              str(round(solution.fun, 2))+", $T_C$ = "+str(round(sol[1], 2))+
                              "$\pm$ "+str(round(Terr, 2))+" eV",fontsize=10)
                    plt.figure()
                    plt.subplot(211)
                    plt.ylabel("$T_i [eV]$",fontsize = 15)
                    plt.grid()
                    plt.title(self.expe_id+", channel "+str(roi)+
                              ", R = "+str(round(R_plot,4))+" m, $\chi^2 = $" +
                              str(round(solution.fun, 2))+", $T_C$ = "+str(round(sol[1], 2))+
                              "$\pm$ "+str(round(Terr, 2))+" eV",fontsize=10)
                    plt.hlines(sol[1], 0, N, color = "red")
                    plt.fill_between(np.arange(N), sol[1] - Terr, sol[1] + Terr,color='red',
                                     alpha=0.2)
                    plt.plot(np.arange(N),T_iters,marker = "o",linestyle="",color="blue")
                    plt.xlim(0,N-1)
                    plt.subplot(212)
                    plt.plot(np.arange(N),chi_iters,marker = "o",linestyle="",color="green")
                    plt.ylabel("$\chi^2$",fontsize = 15)
                    plt.xlabel("number of iterations",fontsize = 15)
                    plt.xlim(0,N-1)
                    plt.grid()
                
            elif(fittype == "CVI"):
                self.B = np.loadtxt(self.supl_data_path+self.expe_id+".txt")
                self.Babs = np.sqrt(np.sum(self.B[:, roi-1]**2))
                divider = np.sqrt(
                    np.sum(self.observation_directions[:, roi-1]**2))*self.Babs
                n = self.observation_directions[:, roi-1]
                self.current_theta = np.arccos(
                    np.sum(n*self.B[:, roi-1])/divider)*180/np.pi
                
                # location where the web service is hosted
                pc_location = 'http://sv-coda-wsvc-28.ipp-hgw.mpg.de:6055'

                # fetching the fine structure of the predefined line
                fine_structure_query = '/getZeeman.json?name=C-VI-5291&B=' + \
                    str(self.Babs)+'&theta1='+str(self.current_theta)
                fine_structure = requests.get(
                    pc_location + fine_structure_query).json()

                self.zc_locations = np.array(fine_structure['wavelengths'])/10
                self.zc_intensities = np.array(fine_structure['amplitude'])
                
                #fitting the ion temperature, and other parameters
                esti = np.array([mu, kbt, A])
                es_chisq = self.C_fitfunc_plot(esti)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
                solution = minimize(self.C_fitfunc, esti, method="Powell",
                                    bounds=((None, None), (0.1, None), (None, None)),
                                    tol=1e-12,options={"maxiter": 2000})

                print(solution)
                if(solution.success == True):
                    #error calculation
                    sol = solution.x
                    Terr, T_iters, chi_iters = self.C_Ti_error_sim_me(sol[0],
                                        sol[1],sol[2],N,fittype,plots=plots)
                    # h = np.array([1e-5, 1e-3, 1e-8])
                    # hessian_error = self.error_from_hesse(sol, solution.fun, h)[1,1]
                    # print("Error based on Hessian matrix (in eV):")
                    # print(hessian_error)
                    R_plot = round(self.dataobj.coordinate( #plotting
                        "Device R")[0][0, (roi-1), 0], 4)
                    self.C_fitfunc_plot(sol)
                    plt.title(self.expe_id+", channel "+str(roi)+
                              ", R = "+str(round(R_plot,4))+" m, $\chi^2 = $" +
                              str(round(solution.fun, 2))+", $T_C$ = "+str(round(sol[1], 2))+
                              "$\pm$ "+str(round(Terr, 2))+" eV",fontsize=10)
                    plt.figure()
                    plt.subplot(211)
                    plt.ylabel("$T_i [eV]$",fontsize = 15)
                    plt.grid()
                    plt.title(self.expe_id+", channel "+str(roi)+
                              ", R = "+str(round(R_plot,4))+" m, $\chi^2 = $" +
                              str(round(solution.fun, 2))+", $T_C$ = "+str(round(sol[1], 2))+
                              "$\pm$ "+str(round(Terr, 2))+" eV",fontsize=10)
                    plt.hlines(sol[1], 0, N, color = "red")
                    plt.fill_between(np.arange(N), sol[1] - Terr, sol[1] + Terr,color='red',
                                     alpha=0.2)
                    plt.plot(np.arange(N),T_iters,marker = "o",linestyle="",color="blue")
                    plt.xlim(0,N-1)
                    plt.subplot(212)
                    plt.plot(np.arange(N),chi_iters,marker = "o",linestyle="",color="green")
                    plt.ylabel("$\chi^2$",fontsize = 15)
                    plt.xlabel("number of iterations",fontsize = 15)
                    plt.xlim(0,N-1)
                    plt.grid()

    def C_Ti_error_sim(self,mu,kbt,A,tstart,tstop,iter_num,scalef,fittype,plots=False):
        """
        A function for simulating uncertainties of possible future ion temperature
        measurements when certain parameters are defined previously (like ion type,
        Zeeman spectra, setting of the spectrometer, etc.). A previous measurement
        is utilized.

        Parameters
        ----------
        mu : expected value for term that constitutes from the Doppler shift of the line, the
        wavelength uncertainty of the spectrometer (float), and the line central wavelength
        kbt : expected value for temperature times Boltzmann constant (float)
        A : expected value for line intensity factor (float)
        tstart : lower end of the interval for which the ion temperature to be 
        calculated in the measurement
        tstop : higher end of the interval for which the ion temperature to be 
        calculated in the measurement
        iter_num : Number of iterations for the Monte Carlo error calculation process (int)
        scalef : scaling number between the measured and to be generated line intensity
        (float)
        fittype : (string) The kind of line that is to be generated. At the moment, 
        it can be "CV" or "CVI".
        plots : (Boolean) wether plot the generated spectra during the error calculation
        process or not. The default is False.
        """
        met = "Powell" #as far as I know that is the best choice for such purpose
        line_param = np.array([mu, kbt, A])

        T_i = np.zeros((iter_num))
        chisq = np.zeros((iter_num))
        
        R_plot=round(self.dataobj.coordinate("Device R")[0][0,(self.current_roi-1),0],4)

        for i in range(iter_num):
            #measurement simulation, then fit iter_num times
            print("Iteration "+str(i))
            self.simulated, self.simulated_error = self.C_line_simulator(mu,
                                    kbt, A,fittype,scalef=scalef, plots=False,sim=True)
            if(plots == True):
                es_chisq = self.C_fitfunc_plot_sim(line_param)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
            solution = minimize(self.C_fitfunc_sim, line_param, method=met,
                                bounds=((None, None), (0.1, None), (None, None)), tol=1e-8,
                                options={"maxiter": 2000})
            if(solution.success == False):
                raise ValueError("Failed T_i fit")
            sol = solution.x
            print(sol[1])
            T_i[i] = sol[1]
            chisq[i] = solution.fun
            if(plots == True):
                self.C_fitfunc_plot_sim(sol)
                plt.title("R = "+str(R_plot)+" m, $\chi^2 = $" +
                          str(round(solution.fun, 6))+", $T_C$ = "+str(round(sol[1], 2)))
        Terr = np.std(T_i)
        plt.figure()
        plt.subplot(211)
        plt.ylabel("$T_i [eV]$",fontsize = 15)
        plt.grid()
        plt.title(self.expe_id+", channel "+str(self.current_roi)+
                  ", R = "+str(round(R_plot,4))+" m, $\chi^2 = $" +
                  str(round(chisq.mean(), 2))+", $T_C$ = "+str(round(np.mean(T_i), 2))+
                  "$\pm$ "+str(round(Terr, 2))+" eV",fontsize=10)
        plt.hlines(np.mean(T_i), 0, iter_num, color = "red")
        plt.fill_between(np.arange(iter_num),np.mean(T_i)-Terr,np.mean(T_i)+Terr,
                         color='red',alpha=0.2)
        plt.plot(np.arange(iter_num),T_i,marker = "o",linestyle="",color="blue")
        plt.xlim(0,iter_num-1)
        plt.subplot(212)
        plt.plot(np.arange(iter_num),chisq,marker = "o",linestyle="",color="green")
        plt.ylabel("$\chi^2$",fontsize = 15)
        plt.xlabel("number of iterations",fontsize = 15)
        plt.xlim(0,iter_num-1)
        plt.grid()
        return Terr
    
    def Ti_error_simulation(self,fittype,roi,wstart,wstop,mu,kbt,A,dslit,
                            t_start,t_stop,bcg,N,simd,simgrid,scalef,plots=False,
                            options = None):
        """
        A function for simulating uncertainties of possible future ion temperature
        measurements using C_Ti_error_sim() for which the specific details are defined
        here. A previous measurement is utilized.

        Parameters
        ----------
        fittype : (string) The kind of line that is to be generated. At the moment, 
        it can be "CV" or "CVI".
        roi : Region Of Interest, in other words spectral channel (int).
        wstart : lower end of the wavelength interval of the spectral line (float)
        wstop : higher end of the wavelength interval of the spectral line (float)
        mu : expected value for term that constitutes from the Doppler shift of the line, the
        wavelength uncertainty of the spectrometer (float), and the line central wavelength
        kbt : expected value for temperature times Boltzmann constant (float)
        A : expected value for line intensity factor (float)
        dslit: slit size in the spectrometer in the measurement
        t_start : lower end of the interval for which the ion temperature to be 
        calculated in the measurement
        t_stop : higher end of the interval for which the ion temperature to be 
        calculated in the measurement
        bcg : A wavelength interval that does not contain any notable spectral lines 
        in the measurement (list with two floats).
        N : Number of iterations for the Monte Carlo error calculation process (int)
        simd : spectrometer slit size in the simulation
        simgrid : spectrometer grid in the simulation
        scalef : scaling number between the measured and to be generated line intensity
        (float)
        plots : (Boolean) wether plot the generated spectra during the error calculation
        process or not. The default is False.
        options : dictionary
            'Supplementary_data_path': string, path of the supplementary data
            (like wavelength calibration)
            'Instrumental_functions_datapath': string, path of the saved instrumental
            functions, which should be .npy files
        """
        default_options = {"Supplementary_data_path":"data"}
        _options = flap.config.merge_options(default_options,options,data_source='W7X_ABES_CXRS')
        try:
            self.supl_data_path = _options['Supplementary_data_path']
            self.instr_funcs_datapath = _options['Instrumental_functions_datapath']
        except (KeyError, TypeError):
            self.supl_data_path = 'data'
            self.instr_funcs_datapath = 'data'
        if(self.campaign == "OP2.1"):
            centre_of_lens = np.array([1.305624, 6.094843, -3.013095])
            roi_pos = np.loadtxt(self.supl_data_path+"ROI_pos.txt")
            self.observation_directions = np.zeros(
                (roi_pos.shape[0], roi_pos.shape[1]))
            for i in range(roi_pos.shape[1]):
                self.observation_directions[:,i] = roi_pos[:, i]-centre_of_lens
            self.current_roi = roi
            self.dslit = dslit
            self.wstart = wstart
            self.wstop = wstop
            self.simd = simd
            self.simgrid = simgrid

            self.active = self.active_passive(roi, t_start, t_stop,
                                    background_interval=bcg, error=True, plotting=False)
            self.active = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(wstart, wstop)})
            if(fittype == "CV"):
                self.B = np.loadtxt(self.supl_data_path+self.expe_id+".txt")
                self.Babs = np.sqrt(np.sum(self.B[:, roi-1]**2))
                divider = np.sqrt(
                    np.sum(self.observation_directions[:, roi-1]**2))*self.Babs
                n = self.observation_directions[:, roi-1]
                self.current_theta = np.arccos(
                    np.sum(n*self.B[:, roi-1])/divider)*180/np.pi
                
                self.errparam = np.array([0.07875014474033012,0.4501855052541737]) 
                #on a H line (486 nm), and on the background, respectively, for 9s
                
                with open(self.supl_data_path+'getZeeman_CV_ROI'+str(self.current_roi)+
                          '.json', 'r') as f:
                    datafile = json.load(f)
                
                self.zc_locations = np.array(datafile['wavelengths'])/10
                self.zc_intensities = np.array(datafile['amplitude'])
                
                Terr = self.C_Ti_error_sim(mu,kbt,A,t_start,t_stop,N,scalef,fittype,plots=plots)
                print("T_i error with the given parameters:")
                print(Terr)
                return Terr
            
            elif(fittype == "CVI"):
                self.B = np.loadtxt(self.supl_data_path+self.expe_id+".txt")
                self.Babs = np.sqrt(np.sum(self.B[:, roi-1]**2))
                divider = np.sqrt(
                    np.sum(self.observation_directions[:, roi-1]**2))*self.Babs
                n = self.observation_directions[:, roi-1]
                self.current_theta = np.arccos(
                    np.sum(n*self.B[:, roi-1])/divider)*180/np.pi
                self.errparam = np.array([0.08126764,0.44721794]) 
                #average ROI1-4 error parameters on Sodium lines for 6s
                
                # location where the web service is hosted
                pc_location = 'http://sv-coda-wsvc-28.ipp-hgw.mpg.de:6055'

                # fetching the fine structure of the predefined line
                fine_structure_query = '/getZeeman.json?name=C-VI-5291&B=' + \
                    str(self.Babs)+'&theta1='+str(self.current_theta)
                fine_structure = requests.get(
                    pc_location + fine_structure_query).json()

                self.zc_locations = np.array(fine_structure['wavelengths'])/10
                self.zc_intensities = np.array(fine_structure['amplitude'])

                Terr = self.C_Ti_error_sim(mu,kbt,A,t_start,t_stop,N,scalef,fittype,plots=plots)
                print("T_i error with the given parameters:")
                print(Terr)
                return Terr