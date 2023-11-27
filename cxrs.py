# -*- coding: utf-8 -*-
"""
Created on Tue Aug 8 15:27:49 2023

@author: bcsillag

Data processing code for Wendelstein 7-X QSI CXRS spectra measured during OP2.1
"""

import requests
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import pandas as pd
import flap


def wavelength_grid_generator_op21(grid, wavelength_setting, roi):
    # loading the calibration coeffs
    calib_array = np.loadtxt("wavelength_calib_2023_"+grid+".txt")
    c0 = calib_array[roi-1, 0]
    c1 = calib_array[roi-1, 2]
    c2 = calib_array[roi-1, 4]
    c3 = calib_array[roi-1, 6]  # writing the out to clear variables
    pix_values = np.arange(0, 1024, 1) - 511
    return c3*pix_values**3 + c2*pix_values**2 + c1*pix_values + c0 + wavelength_setting


def interval_shift(expe_id):
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
        shift = 0.04764529058116232
        # shift = -0.07229458917835671#-0.0592
    elif(expe_id[:8] == "20230328"):
        shift = 0.051633
    elif(expe_id[:8] == "20230330"):
        shift = 0.05275551102204408  # 0.051633

    return shift

def spectral_error_calc_op21(spec):
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

    def __init__(self, data_source, expe_id, get_data="by shotID",
                 campaign="OP2.1", spatcal=True, time_correction=True):
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

        if(data_source == 'W7X_WEBAPI' and get_data == "by shotID" and
                campaign == "OP2.1"):
            self.dataobj = flap.get_data(data_source,
                                         name='Test/raw/W7X/QSI/cxrs_DATASTREAM/0/Images/',
                                         exp_id=expe_id,
                                         options={'Scale Time': True,
                                                  'Cache Data': False},
                                         object_name='QSI_spectral_data')
        else:
            raise ValueError("Undefined data source or loading process.")

        if(self.campaign == "OP2.1"):
            self.APDCAM_timerange = [1, 40]

        if(spatcal == True and self.campaign == "OP2.1"):
            # changing to ROI coord
            roi_coord_flap = flap.Coordinate(name='ROI', unit="1",
                                             mode=flap.CoordinateMode(equidistant=False), shape=(6),
                                             values=["P01", "P02", "P03",
                                                     "P04", "P05", "P06"],
                                             dimension_list=[1])
            self.dataobj.add_coordinate_object(roi_coord_flap)
            self.dataobj.del_coordinate("Coord 1")

            # adding coordinate Device R
            R_coord_flap = flap.Coordinate(name='Device R', unit="m",
                                           mode=flap.CoordinateMode(equidistant=False), shape=(6),
                                           values=[6.175608381963992, 6.18531258540187, 6.196755395927777,
                                                   6.207903188440903, 6.22187706429289, 6.241915857714211],
                                           dimension_list=[1])
            self.dataobj.add_coordinate_object(R_coord_flap)

            # correcting time coordinate
            c = self.dataobj.get_coordinate_object("Time")
            c.values = c.values + 1
            c = self.dataobj.get_coordinate_object("Time")

        elif(spatcal == True and self.campaign != "OP2.1"):
            raise ValueError(
                "Spatial for that campaign has not been implemented yet.")

        if(self.campaign == "OP2.1" and time_correction == True):
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

            # correcting the timescales
            el = interval_shift(self.expe_id)
            c = d_beam_on.get_coordinate_object("Time")
            c.start = c.start + el
            c = d_beam_off.get_coordinate_object("Time")
            c.start = c.start + el

            self.d_beam_on = d_beam_on
            self.d_beam_off = d_beam_off

    def wavelength_calibration(self, grid=None, wavelength_setting=None):
        if(self.campaign == "OP2.1"):
            settings_df = pd.read_excel(
                "OP21/discharge_table_for_spectral_measurements.xlsx")
            settings = settings_df.loc[settings_df["Discharge"] == float(
                self.expe_id)]
            self.grid = str(settings["Grid [g/mm]"].to_numpy()[0])+"g_per_mm"
            self.wavelength_setting = float(
                settings["Wavelength setting [nm]"])
            if(self.campaign == "OP2.1"):
                wl_values = np.zeros((6, 1024))
                for i in range(1, 7):
                    wl_values[i-1, :] = wavelength_grid_generator_op21(
                        self.grid, self.wavelength_setting, i)

                wl_coord_flap = flap.Coordinate(name='Wavelength', unit="nm",
                                                mode=flap.CoordinateMode(equidistant=False), shape=(6, 1024),
                                                values=wl_values, dimension_list=[1, 2])
                self.dataobj.add_coordinate_object(wl_coord_flap)
                self.dataobj.del_coordinate("Coord 2")

        else:
            self.grid = grid
            self.wavelength_setting = wavelength_setting

    def slice_by_wl(self, roi, wavelength):

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
                "Temporal evolution of the pixel() at wavelength "+str(wavelength)+" nm")
        elif(type(wavelength) == list):
            tit = "Avg. temporal evolution of pixels between wavelength ["
            plt.title(tit+str(wavelength[0])+", "+str(wavelength[1])+"] nm")

    def passive(self, roi, t_start, t_stop):

        ROI1 = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                                "Time": flap.Intervals(t_start, t_stop)})
        avg_spectrum = ROI1.slice_data(summing={"Time": "Mean"})

        plt.figure()
        avg_spectrum.plot(axes="Wavelength")
        plt.title(self.expe_id, fontsize=15)
        plt.title(self.expe_id+", averaged spectra, ROI = P0" +
                  str(roi), fontsize=15)
        plt.grid()

    def autocorr(self, roi, t_start, t_stop, lstart, lstop):
        ROI1 = self.dataobj.slice_data(
            slicing={"ROI": "P0"+str(roi), "Time": flap.Intervals(t_start, t_stop)})
        dt = np.mean(ROI1.coordinate("Time")[
                     0][1:, 0]-ROI1.coordinate("Time")[0][:-1, 0])
        line = ROI1.slice_data(slicing={"Wavelength": flap.Intervals(lstart, lstop)},
                               summing={"Wavelength": "Sum"})
        c = line.get_coordinate_object("Time")
        c.mode.equidistant = True
        c.shape = line.data.shape
        c.step = dt
        c.start = 0.0
        c.dimension_list = [0]
        t = ROI1.coordinate("Time")[0][:, 0]
        norm_sig=(line.data-np.mean(line.data))/(np.std(line.data)*np.sqrt(len(line.data)))
        tcorr = t-t.mean()
        corr = np.correlate(norm_sig, norm_sig, mode="same")

        fs = 15
        plt.figure()
        plt.plot(tcorr, corr, "bo-")
        plt.xlabel(r'$\tau [s]$', fontsize=fs)
        plt.title(
            self.expe_id+", Auto-correlation, $\lambda = ["+str(lstart)+", "+str(lstop)+"]$ nm")
        plt.grid()

    def tshif(self, roi, t_start, t_stop, wstart, wstop, N, background_interval=[0]):
        if(self.campaign == "OP2.1"):
            tsh = np.linspace(-0.075, 0.075, N)
            lineint = np.zeros((N))
            for i in range(len(tsh)):
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

                c = d_beam_on.get_coordinate_object("Time")
                c.start = c.start + tsh[i]
                c = d_beam_off.get_coordinate_object("Time")
                c.start = c.start + tsh[i]

                ROI1 = self.dataobj.slice_data(
                    slicing={"ROI": "P0"+str(roi), "Time": flap.Intervals(t_start, t_stop)})
                if(background_interval != [0]):
                    ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                                                  "Wavelength": flap.Intervals(background_interval[0],
                                                                                               background_interval[1])},
                                                         summing={"Wavelength": "Mean", "Time": "Mean"})
                    ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data

                s_on = ROI1.slice_data(slicing={'Time': d_beam_on}, summing={
                                       "Rel. Time in int(Time)": "Mean"})
                s_off = ROI1.slice_data(slicing={'Time': d_beam_off}, summing={
                                        "Rel. Time in int(Time)": "Mean"})
                s_on_intervals_full = ROI1.slice_data(
                    slicing={'Time': d_beam_on})
                s_off_intervals_full = ROI1.slice_data(
                    slicing={'Time': d_beam_off})
                s_on_intervals = s_on_intervals_full.data[:, 1:, ].mean(axis=1)
                s_off_intervals = s_off_intervals_full.data[:, 1:, ].mean(
                    axis=1)
                s_on_data = s_on_intervals.mean(axis=1)
                s_off_data = s_off_intervals.mean(axis=1)
                s_on.data = s_on_data
                s_off.data = s_off_data

                s_subs = s_on
                s_subs.data = s_on.data-s_off.data
                lineint[i] = s_subs.slice_data(slicing={"Wavelength": flap.Intervals(wstart, wstop)},
                                               summing={"Wavelength": "Mean"}).data

            fs = 15
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
                       error=False, plotting=True):
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

            if(plotting == True):
                plt.figure()
                s_on.plot(axes="Wavelength")
                s_off.plot(axes="Wavelength")
                legend = ["beam on", "beam off"]
                plt.legend(legend)
                plt.title(self.expe_id+", active spectrum, ROI = P0"+str(roi)+
                          ", ["+str(t_start)+","+str(t_stop)+"] s",fontsize=15)
                plt.ylabel("Intensity [a.u.]")
                plt.grid()

            s_subs = s_on
            s_subs.data = s_on.data-s_off.data  # gaining the active spectrum
            if(error == True):
                s_subs.error = np.sqrt(error_on**2 + error_off**2)

            if(plotting == True):
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
        if(self.campaign == "OP2.1"):
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
            el = interval_shift(self.expe_id)
            # correcting the timescales
            c = d_beam_on.get_coordinate_object("Time")
            c.start = c.start + el
            c = d_beam_off.get_coordinate_object("Time")
            c.start = c.start + el

            # slicing the data
            ROI1 = self.dataobj.slice_data(
                slicing={"ROI": "P0"+str(roi), "Time": flap.Intervals(t_start, t_stop)})
            ROI1 = ROI1.slice_data(
                slicing={"Wavelength": flap.Intervals(lstart, lstop)})

            if(background_interval != [0]):
                ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                                              "Wavelength": flap.Intervals(background_interval[0], background_interval[1])},
                                                     summing={"Wavelength": "Mean", "Time": "Mean"})
                ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data
            s_on = ROI1.slice_data(slicing={'Time': d_beam_on}, summing={
                                   "Rel. Time in int(Time)": "Mean"})
            s_off = ROI1.slice_data(slicing={'Time': d_beam_off}, summing={
                                    "Rel. Time in int(Time)": "Mean"})

            lambd = s_on.coordinate("Wavelength")[0]
            def gaus(x, A, s, mu): return A*np.e**(-(((x-mu)**2)/s**2))

            popton, pcovon = curve_fit(gaus, lambd, s_on.data, p0=[
                                       max(s_on.data), 0.1, lambd.mean()])
            poptoff, pcovoff = curve_fit(gaus, lambd, s_off.data, p0=[
                                         max(s_off.data), 0.1, lambd.mean()])
            if(plotting == True):
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
                if(background_interval != [0]):
                    ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(roi),
                                                                  "Wavelength": flap.Intervals(background_interval[0], background_interval[1])},
                                                         summing={"Wavelength": "Mean", "Time": "Mean"})
                    ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data
                s_on = ROI1.slice_data(slicing={'Time': self.d_beam_on}, summing={
                                       "Rel. Time in int(Time)": "Mean"})
                s_off = ROI1.slice_data(slicing={'Time': self.d_beam_off}, summing={
                                        "Rel. Time in int(Time)": "Mean"})
                int_on = list(s_on.data)
                int_off = list(s_off.data)

                err = np.array(error_on + error_off)
                intensity = np.array(int_on + int_off)

                err = err[intensity.argsort()]
                intensity = intensity[intensity.argsort()]

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
            if(background_interval != [0]):
                ROI1_witbg = self.dataobj.slice_data(slicing={"ROI": "P0"+str(ROI),
                                                              "Wavelength": flap.Intervals(background_interval[0], background_interval[1])},
                                                     summing={"Wavelength": "Mean", "Time": "Mean"})
                ROI1.data[:, :] = ROI1.data[:, :] - ROI1_witbg.data
            s_on = ROI1.slice_data(slicing={'Time': self.d_beam_on}, summing={
                                   "Rel. Time in int(Time)": "Mean"})
            s_off = ROI1.slice_data(slicing={'Time': self.d_beam_off}, summing={
                                    "Rel. Time in int(Time)": "Mean"})
            int_on = list(s_on.data)
            int_off = list(s_off.data)

            err = np.array(error_on + error_off)
            intensity = np.array(int_on + int_off)

            err = err[intensity.argsort()]
            intensity = intensity[intensity.argsort()]

            def sq(x, a, b): return a*np.sqrt(x) + b
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

    def C_line_generator(self, mu_add, kbt, A,sim=False):
        wl_grid = None
        if(sim==False):
            wl_values = wavelength_grid_generator_op21(
                self.grid, self.wavelength_setting, self.current_roi)
            # slicing the wavelength grid
            wl_grid0 = wl_values[wl_values > self.wstart]
            wl_grid = wl_grid0[self.wstop > wl_grid0]
        elif(sim==True):
            wl_values = wavelength_grid_generator_op21(
                self.simgrid, self.wavelength_setting, self.current_roi)
            # slicing the wavelength grid
            wl_grid0 = wl_values[wl_values > self.wstart]
            wl_grid = wl_grid0[self.wstop > wl_grid0] 
        # Doppler shift + calibration uncertainty
        locations = self.zc_locations + mu_add
        projection = np.zeros((wl_grid.shape[0]))
        mu = np.dot(locations, self.zc_intensities)/sum(self.zc_intensities)
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
        gauss=lambda lambd,s,A,up,down:A*np.e**(-(((lambd-(up+down)/2)**2)/s**2))
        gaussian = gauss(wl_grid, s, A, self.wstop, self.wstart)
        doppler_spectrum = np.convolve(projection, gaussian, mode="same")
        instr = None
        if(sim==False):
            # convolution with instrumental function
            instr = np.load("instr_funcs/"+self.grid+"_P0"+str(self.current_roi) +
                            "_"+str(int(self.dslit))+"micron_slit.npy").ravel()
        elif(sim==True):
            # convolution with instrumental function
            if(self.simd == 100):
                instr = np.load("instr_funcs/"+self.simgrid+"_P0"+str(self.current_roi) +
                            "_"+str(int(self.simd))+"micron_slit.npy").ravel()
            else:
                instr = np.load("instr_funcs/"+self.simgrid+"_P03" +
                            "_"+str(int(self.simd))+"micron_slit.npy").ravel()
                
        # instr = np.load("instr_funcs/"+self.grid+"_P03_"+str(int(self.dslit))+"micron_slit.npy").ravel()
        complete_spectrum = np.convolve(doppler_spectrum, instr, mode="same")

        return complete_spectrum

    def CVI_fitfunc(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A)
        C = (modelled - self.active.data)/self.active.error
        return (np.dot(C, C) - 3) / self.active.data.shape[0]
    
    def CVI_fitfunc_sim(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A,sim = True)
        C = (modelled - self.simulated)/self.simulated_error
        return (np.dot(C, C) - 3) / self.simulated.shape[0]

    def CVI_fitfunc_plot(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A)
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
    
    def CV_fitfunc(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A)
        C = (modelled - self.active.data)/self.active.error
        return (np.dot(C, C) - 3) / self.active.data.shape[0]
    
    def CV_fitfunc_sim(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A,sim = True)
        C = (modelled - self.simulated)/self.simulated_error
        return (np.dot(C, C) - 3) / self.simulated.shape[0]

    def CV_fitfunc_plot(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A)
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
        sol_p = sol.copy()
        sol_n = sol.copy()
        sol_p[i] = sol_p[i] + h
        sol_n[i] = sol_n[i] - h
        fp1 = self.CVI_fitfunc(sol_p)
        fm1 = self.CVI_fitfunc(sol_n)
        return (fp1 + fm1 - 2*fmin)/h**2

    def partial_cross_deriv(self,sol, fmin, i, j, hi, hj):
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

        fpp = self.CVI_fitfunc(sol_pp)
        fpn = self.CVI_fitfunc(sol_pn)
        fnp = self.CVI_fitfunc(sol_np)
        fnn = self.CVI_fitfunc(sol_nn)

        return (fpp-fpn-fnp+fnn) / (4*hi*hj)

    def hesse(self,sol, fmin, h):
        hess = np.matrix(np.zeros((3, 3)))
        for i in range(3):
            for j in range(3):
                if(i == j):
                    hess[i, i] = self.second_deriv(sol, fmin, h[i], i)
                else:
                    hess[i, j] = self.partial_cross_deriv(sol, fmin, i, j, h[i], h[j])
        return hess

    def error_from_hesse(self,sol,fmin, h):
        alpha = self.hesse(sol, fmin, h)
        I = np.linalg.inv(alpha)
        return np.sqrt(abs(I))
    
    def CVI_fitfunc_plot_sim(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A,sim = True)
        C = (modelled - self.simulated)/self.simulated_error

        wl_values = wavelength_grid_generator_op21(
            self.simgrid, self.wavelength_setting, self.current_roi)  # loading the wavelength grid
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
    
    def CV_fitfunc_plot_sim(self, esti):
        mu_add = esti[0]
        kbt = esti[1]
        A = esti[2]
        modelled = self.C_line_generator(mu_add, kbt, A,sim = True)
        C = (modelled - self.simulated)/self.simulated_error

        wl_values = wavelength_grid_generator_op21(
            self.simgrid, self.wavelength_setting, self.current_roi)  # loading the wavelength grid
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
    
    def CVI_line_simulator(self,mu_add,kbt,A,tstart,tstop,scalef=None,plots=False,sim=False):
        gaus = lambda x, A, s, mu: A*np.e**(-(((x-mu)**2)/s**2))
        if(sim==False):
            measured = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(self.wstart, self.wstop)})
            lamb = measured.coordinate("Wavelength")[0]
    
            popt, pcov = curve_fit(gaus, lamb, measured.data, p0=[
                                   max(measured.data), 0.1, lamb.mean()])
            if(plots == True):
                wl_values = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, self.current_roi)  # loading the wavelength grid
                wl_grid0 = wl_values[wl_values > self.wstart]  # slicing the wavelength grid
                lambd = wl_grid0[self.wstop > wl_grid0]
                fs = 15
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, gaus(lambd, *popt))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")
    
            calculated = self.C_line_generator(mu_add, kbt, A)
            calculated = popt[0]*calculated/max(calculated)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            err = measured.error
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "o")
                plt.errorbar(lambd, calculated, err)
                plt.grid()
    
            calculated = np.random.normal(loc=calculated, scale=err)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            return calculated, err
        
        if(sim==True):
            measured = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(self.wstart, self.wstop)})
            lamb = measured.coordinate("Wavelength")[0] 
    
            popt, pcov = curve_fit(gaus, lamb, measured.data, p0=[
                                   max(measured.data), 0.1, lamb.mean()])
            if(plots == True):
                wl_values = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, self.current_roi)  # loading the wavelength grid
                wl_grid0 = wl_values[wl_values > self.wstart]  # slicing the wavelength grid
                lambd = wl_grid0[self.wstop > wl_grid0]
                fs = 15
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, gaus(lambd, *popt))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")
    
            calculated = self.C_line_generator(mu_add, kbt, A,sim=True)
            gridfac = None
            
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
            
            calculated = gridfac*scalef*popt[0]*calculated/max(calculated)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
            sq = lambda x, a, b: a*np.sqrt(x) + b
            err = sq(calculated*15, self.errparam[0], self.errparam[1]) #assuming 15% modulation

            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "o")
                plt.errorbar(lambd, calculated, err)
                plt.grid()

            calculated = np.random.normal(loc=calculated, scale=err)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            return calculated, err
        
    def CV_line_simulator(self,mu_add,kbt,A,tstart,tstop,scalef=None,plots=False,sim=False):
        gaus = lambda x, A, s, mu: A*np.e**(-(((x-mu)**2)/s**2))
        if(sim==False):
            measured = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(self.wstart, self.wstop)})
            lamb = measured.coordinate("Wavelength")[0]
    
            popt, pcov = curve_fit(gaus, lamb, measured.data, p0=[
                                   max(measured.data), 0.1, lamb.mean()])
            if(plots == True):
                wl_values = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, self.current_roi)  # loading the wavelength grid
                wl_grid0 = wl_values[wl_values > self.wstart]  # slicing the wavelength grid
                lambd = wl_grid0[self.wstop > wl_grid0]
                fs = 15
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, gaus(lambd, *popt))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")
    
            calculated = self.C_line_generator(mu_add, kbt, A)
            calculated = popt[0]*calculated/max(calculated)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            err = measured.error
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "o")
                plt.errorbar(lambd, calculated, err)
                plt.grid()
    
            calculated = np.random.normal(loc=calculated, scale=err)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            return calculated, err
        
        if(sim==True):
            measured = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(self.wstart, self.wstop)})
            lamb = measured.coordinate("Wavelength")[0] 
    
            popt, pcov = curve_fit(gaus, lamb, measured.data, p0=[
                                   max(measured.data), 0.1, lamb.mean()])
            if(plots == True):
                wl_values = wavelength_grid_generator_op21(
                    self.grid, self.wavelength_setting, self.current_roi)  # loading the wavelength grid
                wl_grid0 = wl_values[wl_values > self.wstart]  # slicing the wavelength grid
                lambd = wl_grid0[self.wstop > wl_grid0]
                fs = 15
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, gaus(lambd, *popt))
                plt.xlabel("Wavelength [nm]", fontsize=fs)
                plt.ylabel("Spectral intensity", fontsize=fs)
                plt.legend(["Data", "Fit"], fontsize=fs-2)
                plt.title(self.expe_id+", Beam on line intensity fit")
    
            calculated = self.C_line_generator(mu_add, kbt, A,sim=True)
            gridfac = 1#None
            
            # baseint = 11230
            # if(self.simgrid == "1200g_per_mm" and self.simd == 100):
            #     gridfac = 1
            # elif(self.simgrid == "1800g_per_mm" and self.simd == 100):
            #     gridfac = 16922/baseint
            # elif(self.simgrid == "2400g_per_mm" and self.simd == 100):
            #     gridfac = 10688/baseint
            # elif(self.simgrid == "1200g_per_mm" and self.simd == 70):
            #     gridfac = 12810/(baseint*1.5)
            # elif(self.simgrid == "1800g_per_mm" and self.simd == 70):
            #     gridfac = 20885/(baseint*1.5)
            # elif(self.simgrid == "2400g_per_mm" and self.simd == 70):
            #     gridfac = 13971/(baseint*1.5)
            # elif(self.simgrid == "1200g_per_mm" and self.simd == 50):
            #     gridfac = 14715/(baseint*2)
            # elif(self.simgrid == "1800g_per_mm" and self.simd == 50):
            #     gridfac = 21279/(baseint*2)
            # elif(self.simgrid == "2400g_per_mm" and self.simd == 50):
            #     gridfac = 14274/(baseint*2)
            # elif(self.simgrid == "1200g_per_mm" and self.simd == 35):
            #     gridfac = 12470/(baseint*2.8)
            # elif(self.simgrid == "1800g_per_mm" and self.simd == 35):
            #     gridfac = 21396/(baseint*2.8)
            # elif(self.simgrid == "2400g_per_mm" and self.simd == 35):
            #     gridfac = 16699/(baseint*2.8)
            # else:
            #     raise ValueError("Wrong grid or slit size.")
            
            calculated = gridfac*scalef*popt[0]*calculated/max(calculated)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
            sq = lambda x, a, b: a*np.sqrt(x) + b
            err = sq(calculated*15, self.errparam[0], self.errparam[1]) #assuming 15% modulation

            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "o")
                plt.errorbar(lambd, calculated, err)
                plt.grid()

            calculated = np.random.normal(loc=calculated, scale=err)
    
            if(plots == True):
                plt.figure()
                plt.plot(lamb, measured.data, "+")
                plt.plot(lambd, calculated, marker="o", color="black")
                plt.grid()
    
            return calculated, err
    
    def CVI_Ti_error_sim_me(self,mu_add,kbt,A,tstart,tstop,iter_num,plots=False):
        met = "Powell"
        line_param = np.array([mu_add, kbt, A])

        T_i = np.zeros((iter_num))
        chisq = np.zeros((iter_num))

        for i in range(iter_num):
            print("Iteration "+str(i))
            self.simulated, self.simulated_error = self.CVI_line_simulator(mu_add, kbt, A, tstart,tstop, plots=False)
            # raise ValueError("stop")
            if(plots == True):
                es_chisq = self.CVI_fitfunc_plot_sim(line_param)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
            solution = minimize(self.CVI_fitfunc_sim, line_param, method=met, bounds=((None, None), (0.1, None), (None, None)), tol=1e-8,
                                options={"maxiter": 2000})
            if(solution.success == False):
                raise ValueError("Failed T_i fit")
            sol = solution.x
            print(sol[1])
            T_i[i] = sol[1]
            chisq[i] = solution.fun
            if(plots == True):
                self.CVI_fitfunc_plot_sim(sol)
                R_plot = round(self.dataobj.coordinate("Device R")[0][0, (self.current_roi-1), 0], 4)
                plt.title("R = "+str(R_plot)+" m, $\chi^2 = $" +
                          str(round(solution.fun, 6))+", $T_C$ = "+str(round(sol[1], 2)))

        return np.std(T_i), T_i, chisq
    
    def CV_Ti_error_sim_me(self,mu_add,kbt,A,tstart,tstop,iter_num,plots=False):
        met = "Powell"
        line_param = np.array([mu_add, kbt, A])

        T_i = np.zeros((iter_num))
        chisq = np.zeros((iter_num))

        for i in range(iter_num):
            print("Iteration "+str(i))
            self.simulated, self.simulated_error = self.CV_line_simulator(mu_add, kbt, A, tstart,tstop, plots=False)
            # raise ValueError("stop")
            if(plots == True):
                es_chisq = self.CV_fitfunc_plot_sim(line_param)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
            solution = minimize(self.CV_fitfunc_sim, line_param, method=met, bounds=((None, None), (0.1, None), (None, None)), tol=1e-8,
                                options={"maxiter": 2000})
            if(solution.success == False):
                raise ValueError("Failed T_i fit")
            sol = solution.x
            print(sol[1])
            T_i[i] = sol[1]
            chisq[i] = solution.fun
            if(plots == True):
                self.CV_fitfunc_plot_sim(sol)
                R_plot = round(self.dataobj.coordinate("Device R")[0][0, (self.current_roi-1), 0], 4)
                plt.title("R = "+str(R_plot)+" m, $\chi^2 = $" +
                          str(round(solution.fun, 6))+", $T_C$ = "+str(round(sol[1], 2)))

        return np.std(T_i), T_i, chisq

    def tempfit(self,fittype,roi,wstart,wstop,mu_add,kbt,A,dslit,t_start,t_stop,bcg,N,plots=False):
        if(self.campaign == "OP2.1"):
            centre_of_lens = np.array([1.305624, 6.094843, -3.013095])
            roi_pos = np.loadtxt("OP21/ROI_pos.txt")
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

            self.active = self.active_passive(roi, t_start, t_stop,
                                              background_interval=bcg, error=True, plotting=False)
            self.active = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(wstart, wstop)})
            if(fittype == "CV"):
                self.B = np.loadtxt("OP21/"+self.expe_id+".txt")
                self.Babs = np.sqrt(np.sum(self.B[:, roi-1]**2))
                divider = np.sqrt(
                    np.sum(self.observation_directions[:, roi-1]**2))*self.Babs
                n = self.observation_directions[:, roi-1]
                self.current_theta = np.arccos(
                    np.sum(n*self.B[:, roi-1])/divider)*180/np.pi
                
                with open('OP21/getZeeman_CV_ROI'+str(self.current_roi)+'.json', 'r') as f:
                    datafile = json.load(f)
                
                self.zc_locations = np.array(datafile['wavelengths'])/10
                self.zc_intensities = np.array(datafile['amplitude'])
                # plt.figure()
                # plt.plot(self.zc_locations,self.zc_intensities,"+")
                
                esti = np.array([mu_add, kbt, A])
                es_chisq = self.CV_fitfunc_plot(esti)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
                solution = minimize(self.CV_fitfunc, esti, method="Powell",
                                    bounds=((None, None), (0.1, None), (None, None)), tol=1e-12,
                                    options={"maxiter": 2000})

                print(solution)
                if(solution.success == True):
                    sol = solution.x
                    R_plot = round(self.dataobj.coordinate(
                        "Device R")[0][0, (roi-1), 0], 4)
                    Terr, T_iters, chi_iters = self.CV_Ti_error_sim_me(sol[0],
                                        sol[1],sol[2],t_start,t_stop,N,plots=plots)
                    h = np.array([1e-5, 1e-3, 1e-8])
                    hessian_error = self.error_from_hesse(sol, solution.fun, h)[1,1]
                    print("Error based on Hessian matrix (in eV):")
                    print(hessian_error)
                    self.CV_fitfunc_plot(sol)
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
                    plt.fill_between(np.arange(N), sol[1] - Terr, sol[1] + Terr,color='red', alpha=0.2)
                    plt.plot(np.arange(N),T_iters,marker = "o",linestyle="",color="blue")
                    plt.xlim(0,N-1)
                    plt.subplot(212)
                    plt.plot(np.arange(N),chi_iters,marker = "o",linestyle="",color="green")
                    plt.ylabel("$\chi^2$",fontsize = 15)
                    plt.xlabel("number of iterations",fontsize = 15)
                    plt.xlim(0,N-1)
                    plt.grid()
                
            elif(fittype == "CVI"):
                self.B = np.loadtxt("OP21/"+self.expe_id+".txt")
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
                
                esti = np.array([mu_add, kbt, A])
                es_chisq = self.CVI_fitfunc_plot(esti)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
                solution = minimize(self.CVI_fitfunc, esti, method="Powell",
                                    bounds=((None, None), (0.1, None), (None, None)), tol=1e-12,
                                    options={"maxiter": 2000})

                print(solution)
                if(solution.success == True):
                    sol = solution.x
                    R_plot = round(self.dataobj.coordinate(
                        "Device R")[0][0, (roi-1), 0], 4)
                    Terr, T_iters, chi_iters = self.CVI_Ti_error_sim_me(sol[0],
                                        sol[1],sol[2],t_start,t_stop,N,plots=plots)
                    h = np.array([1e-5, 1e-3, 1e-8])
                    hessian_error = self.error_from_hesse(sol, solution.fun, h)[1,1]
                    print("Error based on Hessian matrix (in eV):")
                    print(hessian_error)
                    self.CVI_fitfunc_plot(sol)
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
                    plt.fill_between(np.arange(N), sol[1] - Terr, sol[1] + Terr,color='red', alpha=0.2)
                    plt.plot(np.arange(N),T_iters,marker = "o",linestyle="",color="blue")
                    plt.xlim(0,N-1)
                    plt.subplot(212)
                    plt.plot(np.arange(N),chi_iters,marker = "o",linestyle="",color="green")
                    plt.ylabel("$\chi^2$",fontsize = 15)
                    plt.xlabel("number of iterations",fontsize = 15)
                    plt.xlim(0,N-1)
                    plt.grid()

    def CVI_Ti_error_sim(self,mu_add,kbt,A,tstart,tstop,iter_num,scalef,plots=False):
        met = "Powell"
        line_param = np.array([mu_add, kbt, A])

        T_i = np.zeros((iter_num))
        chisq = np.zeros((iter_num))
        
        R_plot = round(self.dataobj.coordinate("Device R")[0][0, (self.current_roi-1), 0], 4)

        for i in range(iter_num):
            print("Iteration "+str(i))
            self.simulated, self.simulated_error = self.CVI_line_simulator(mu_add,
                                    kbt, A, tstart,tstop,scalef=scalef, plots=False,sim=True)
            # raise ValueError("stop")
            if(plots == True):
                es_chisq = self.CVI_fitfunc_plot_sim(line_param)
                plt.title("$\chi^2 = $"+str(round(es_chisq, 6)))
            solution = minimize(self.CVI_fitfunc_sim, line_param, method=met,
                                bounds=((None, None), (0.1, None), (None, None)), tol=1e-8,
                                options={"maxiter": 2000})
            if(solution.success == False):
                raise ValueError("Failed T_i fit")
            sol = solution.x
            print(sol[1])
            T_i[i] = sol[1]
            chisq[i] = solution.fun
            if(plots == True):
                self.CVI_fitfunc_plot_sim(sol)
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
        plt.fill_between(np.arange(iter_num), np.mean(T_i) - Terr, np.mean(T_i) + Terr,color='red', alpha=0.2)
        plt.plot(np.arange(iter_num),T_i,marker = "o",linestyle="",color="blue")
        plt.xlim(0,iter_num-1)
        plt.subplot(212)
        plt.plot(np.arange(iter_num),chisq,marker = "o",linestyle="",color="green")
        plt.ylabel("$\chi^2$",fontsize = 15)
        plt.xlabel("number of iterations",fontsize = 15)
        plt.xlim(0,iter_num-1)
        plt.grid()
        return Terr
    
    def Ti_error_simulation(self,fittype,roi,wstart,wstop,mu_add,kbt,A,dslit,
                            t_start,t_stop,bcg,N,simd,simgrid,scalef,plots=False):
        if(self.campaign == "OP2.1"):
            centre_of_lens = np.array([1.305624, 6.094843, -3.013095])
            roi_pos = np.loadtxt("OP21/ROI_pos.txt")
            self.observation_directions = np.zeros(
                (roi_pos.shape[0], roi_pos.shape[1]))
            for i in range(roi_pos.shape[1]):
                self.observation_directions[:,i] = roi_pos[:, i]-centre_of_lens
            self.B = np.loadtxt("OP21/"+self.expe_id+".txt")
            self.Babs = np.sqrt(np.sum(self.B[:, roi-1]**2))
            divider = np.sqrt(
                np.sum(self.observation_directions[:, roi-1]**2))*self.Babs
            n = self.observation_directions[:, roi-1]
            self.current_theta = np.arccos(
                np.sum(n*self.B[:, roi-1])/divider)*180/np.pi
            self.current_roi = roi
            self.dslit = dslit
            self.wstart = wstart
            self.wstop = wstop
            self.simd = simd
            self.simgrid = simgrid
            # self.errparam = np.array([0.05397319090280941,0.5316958506860017]) 
            # #ROI4 CVI line error between 528.7 and 529.3 for 6s
            
            self.errparam = np.array([0.08126764,0.44721794]) 
            #average ROI1-4 error parameters on Sodium lines for 6s

            self.active = self.active_passive(roi, t_start, t_stop,
                                              background_interval=bcg, error=True, plotting=False)
            self.active = self.active.slice_data(
                slicing={"Wavelength": flap.Intervals(wstart, wstop)})
            if(fittype == "CVI"):
                # location where the web service is hosted
                pc_location = 'http://sv-coda-wsvc-28.ipp-hgw.mpg.de:6055'

                # fetching the fine structure of the predefined line
                fine_structure_query = '/getZeeman.json?name=C-VI-5291&B=' + \
                    str(self.Babs)+'&theta1='+str(self.current_theta)
                fine_structure = requests.get(
                    pc_location + fine_structure_query).json()

                self.zc_locations = np.array(fine_structure['wavelengths'])/10
                self.zc_intensities = np.array(fine_structure['amplitude'])

                Terr = self.CVI_Ti_error_sim(mu_add,kbt,A,t_start,t_stop,N,scalef,plots=plots)
                print("T_i error with the given parameters:")
                print(Terr)
                return Terr