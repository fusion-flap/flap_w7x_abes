#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:11:06 2024

@author: IPP-AD\drefy
"""

import warnings
import matplotlib.pyplot as plt
import flap
import flap_w7x_abes
import requests, json
import re
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.colors import TwoSlopeNorm


flap_w7x_abes.register()
try:
    import flap_w7x_webapi as webapi
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
except ModuleNotFoundError:
    warnings.warn("flap_w7x_webapi not found, unable to download cxrs data")
   
class automated_cxrs_calculation():
    def __init__(self, exp_id, verbose=True,Wdia_threshold=100,temp_resolution=1):       
        if (exp_id is None):
            raise ValueError('exp_id should be set for W7X ABES.')
        self.exp_id=exp_id
        self.temp_resolution=temp_resolution
        if int(self.exp_id[:4]) > 2023:
            self.campaign = "OP2.2"
        else:
            self.campaign = "OP2.1"
        self.verbose=verbose
        self.Wdia_threshold=Wdia_threshold
        self.problem_override=True
        if self.check_if_shot_is_ok_for_evaluation() == False:
           warnings.warn(f"No sufficient data for Experiment ID: {self.exp_id}")

        self.collect_parameters_for_spectral_evaluation()
        self.plot_lightprofile()
        self.calc_temperature_profile()

    def check_if_shot_is_ok_for_evaluation(self):
        # check if there was plasma
        print(f"Start preliminary analysis for shot #"+self.exp_id)
        if flap.get_data('W7X_WEBAPI', name='Wdia', exp_id=self.exp_id, options={'Scale Time': True,'Cache Data': True}).data.max() > self.Wdia_threshold:
            self.plasma_ok=True
            print(f"Plasma discharge looks ok based on Wdia (>{self.Wdia_threshold}kJ)")
        else:
            self.plasma_ok=False
            warnings.warn(f"Wdia is below the threshold: {self.Wdia_threshold} kJ (Experiment ID: {self.exp_id})")
            return False   
        
        if self.get_log() == True:
            # check if there was good ABES data and CXRS 
            if self.check_verdicts(self.logbook) == True:
                print("Data quality is ok for CXRS evaluation according to the logbook.")
                self.diagnostics_ok=True
            else:
                self.diagnostics_ok=False
                warnings.warn("Problems with the diagnostics")
                return False
            
            # check if shot was slow modulation
            if self.check_chopper_mode(self.logbook) == True:
                print("Chopper mode 'Camera'. This is ok for CXRS evaluation.")
                self.chopper_ok=True
            else:
                self.chopper_ok=False
                warnings.warn("Problems with the chopping mode")
                return False
        else:
            print("No logbook entry for the shot, lets assume that it was okay...")
            self.diagnostics_ok=None
            self.chopper_ok=None
        
        # self.check_line_modulation()
        return True
    def transform_id(self,input_id):
        """
        Transforms an input ID from '20241204.049' format to 'XP_20241204.49'.
        """
        parts = input_id.split('.')
        date_part = parts[0]
        suffix_part = str(int(parts[1]))  
    
        transformed_id = f"XP_{date_part}.{suffix_part}"
        return transformed_id
    
    def get_log(self):
        id = self.transform_id(self.exp_id)
        logbook_path = 'https://w7x-logbook.ipp-hgw.mpg.de/api/log/QSI/'
        res = requests.get(logbook_path+id)# check if there was good ABES data and CXRS 
        
        if res.status_code == 200:
            self.logbook = res.json()
            return True
        else:
            return False
    
    def check_verdicts(self,log):
        """
        Checks the 'Verdict' values of the APDCAM and CXRS logs.
        Returns False if either is not 'OK'; otherwise, returns True.
        """
        description = log['_source']['description']
        subsystems = description.split('\n---')
    
        targets = ["APDCAM", "CXRS"]
    
        for subsystem in subsystems:
            lines = subsystem.strip().split('\n')
            title = lines[0].strip()
            verdict_line = next((line for line in lines if 'Verdict:' in line), None)
    
            if verdict_line:
                verdict = verdict_line.split(':')[-1].strip()
                if any(target in title for target in targets) and verdict != "OK":
                    print(f"Error: {title} -> {verdict}")  # Debugging output
                    return False
    
        return True
    
    def check_chopper_mode(self,log):
        """
        Checks the 'Chopper' mode in the APDCAM log.
        Returns True if the mode is 'Camera'; False if it is 'Timed' or anything else.
        """
        description = log['_source']['description']
        subsystems = description.split('\n---')
    
        for subsystem in subsystems:
            lines = subsystem.strip().split('\n')
            title = lines[0].strip()
    
            if "APDCAM" in title:
                chopper_line = next((line for line in lines if 'Chopper:' in line), None)
    
                if chopper_line:
                    try:
                        chopper_mode = chopper_line.split(':')[-1].strip()
                        if "Camera" in chopper_mode:
                            # print("Chopper mode 'Camera'. This is ok for CXRS evaluation.")
                            return True
                        elif "Timed" in chopper_mode:
                            print("Chopper mode contains 'Timed'. This is not ok for CXRS evaluation.")  # Debugging output
                            return False
                        else:
                            print(f"Chopper mode is unknown: {chopper_mode}")  # Debugging output
                            return False
                    except IndexError:
                        print("Invalid chopper_line format. Unable to parse chopper mode.")  # Debugging output
                        return False
                else:
                    print("Chopper line is empty or None.")  # Debugging output
                    return False
                
        print(f"Chopper mode is unknown: {chopper_mode}")
        return False

    def collect_parameters_for_spectral_evaluation(self):
        
        # start time of beam & plasma
        # end time of beam & plasma
        def find_plasma_and_beam_time_overlap():
            Wdia=flap.get_data('W7X_WEBAPI', name='Wdia', exp_id=self.exp_id,options={'Scale Time': True,'Cache Data': False})
            crossings=find_threshold_crossings(Wdia.data, self.Wdia_threshold)
            if len(crossings) == 2:
                self.eval_tstart=Wdia.get_coordinate_object('Time').values[crossings[0][0]]
                tstop=Wdia.get_coordinate_object('Time').values[crossings[1][0]]
                max_tstop=32 # maximum length of the beam injection
                self.eval_tstop=tstop if tstop < max_tstop else max_tstop
                print(f"Calculation for the time range: {self.eval_tstart:.1f}-{self.eval_tstop:.1f} s.")  
            else:
                print(f"There is a problem with this shot, Wdia has crossed {self.Wdia_threshold}kA several times...")
                self.eval_tstart=0
                self.eval_tstop=10
                print(f"Using default {self.eval_tstart}-{self.eval_tstop} s for evaluation.")
        
        def calculate_evaluation_time_ranges():
           
            # Generate the intervals
            tstart_values = np.arange(self.eval_tstart, self.eval_tstop, self.temp_resolution)
            tstop_values = tstart_values + self.temp_resolution
            
            # Filter intervals to ensure they do not exceed eval_tstop
            valid_indices = tstop_values <= self.eval_tstop
            tstart_values = tstart_values[valid_indices]
            tstop_values = tstop_values[valid_indices]
            
            # Combine into a vector of (tstart, tstop) pairs
            intervals = list(zip(tstart_values, tstop_values))
            
            # Print the intervals
            self.intervals=np.array(intervals)
        
        def find_threshold_crossings(signal, threshold):
            crossings = []
            # Calculate the difference between the signal and the threshold
            above_threshold = signal > threshold
            # Find where the signal changes from below to above or above to below the threshold
            for i in range(1, len(above_threshold)):
                if above_threshold[i] != above_threshold[i - 1]:  # Check if crossing happened
                    direction = "up" if above_threshold[i] else "down"
                    crossings.append((i - 1, direction))
            return crossings
        
        # def extract_cxrs_parameters():
        #     """
        #     Extract CXRS parameters from a given log data dictionary.
          
        #     Returns:
        #         dict: A dictionary with extracted CXRS parameters.
        #     """
        #     try:
        #     # Extract the log description from the dictionary
        #         log_text = self.logbook.get('_source', {}).get('description', '')
        #         if not log_text:
        #             return {"Error": "Log description not found."}
            
        #         # Define a pattern to identify the CXRS section
        #         cxrs_pattern = r"---CXRS AUTOLOG---(.*?)---"
            
        #         # Search for the CXRS section in the log
        #         cxrs_section = re.search(cxrs_pattern, log_text, re.DOTALL)
        #         if not cxrs_section:
        #             return {"Error": "CXRS section not found in the log."}
            
        #         # Extracted section text
        #         cxrs_text = cxrs_section.group(1)
            
        #         # Define patterns to extract specific parameters
        #         patterns = {
        #             "Number of channels": r"Number of channels:\s*(\d+)",
        #             "Number of images": r"Number of images:\s*(\d+)",
        #             "Central wavelength": r"Central wavelength:\s*([\d.]+)",
        #             "Grating": r"Grating:\s*(\d+g_per_mm)",
        #             "Slit width": r"Slit width:\s*([\d.-]+)",
        #             "Major radii": r"Major radii\s*([\d.-]+-[\d.-]+)m",
        #             "Remark": r"Remark:\s*(.*)"
        #         }
            
        #         # Extract parameters using the defined patterns
        #         cxrs_parameters = {}
        #         for key, pattern in patterns.items():
        #             match = re.search(pattern, cxrs_text)
        #             if match:
        #                 cxrs_parameters[key] = match.group(1)
        #             else:
        #                 cxrs_parameters[key] = "Not found"
            
        #             self.cxrs_parameters=cxrs_parameters
        #     except:
        #         print("Error: Log description not found.")
        #         return {"Error": "Log description not found."}
        
        def get_cxrs_parameters():
            flap.config.read()
            d = flap.get_data('W7X_ABES_CXRS', name='CXRS_TEST',exp_id=self.exp_id,object_name=f"CXRS_{self.exp_id}")
            self.cxrs_parameters=d.config
            
        
        def parameter_definition():
            try:
                if float(self.cxrs_parameters['Central wavelength']) == 529:
                    self.grating=self.cxrs_parameters['Grating']
                    self.wavelength_setting=float(self.cxrs_parameters['Central wavelength'])
                    self.wstart = 528.25 #CVI line
                    self.wstop = 529.75
                    self.backg = [520,521]
                    self.fittype = "CVI"
                    self.mu = 528.93 # doppler + calibration shift, kozponti hullamhossz is jo
                    self.kbt = 50 # 100-200 kozott jol megy
                    self.A = 1.5e-03 # plotfit fv. - kb mekkora az elso illesztett es mert intenzitas aranya
                    self.itern = 10 # hibaszamolas statisztika
                    if float(self.cxrs_parameters['Slit width']) != -1:
                        self.dslit = float(self.cxrs_parameters['Slit width'])
                    else:
                        self.dslit = 150 # file nevben benne van
                    
                elif float(self.cxrs_parameters['Central wavelength']) == 495:
                    self.grating=self.cxrs_parameters['Grating']
                    self.wavelength_setting=float(self.cxrs_parameters['Central wavelength'])
                    self.wstart = 494
                    self.wstop = 495.25
                    self.backg = [488.5,491.5]
                    self.fittype = "CV"
                    self.mu = 494.72
                    self.kbt = 100
                    self.A = 7e-04
                    self.itern = 5
                    if float(self.cxrs_parameters['Slit width']) != -1:
                        self.dslit = float(self.cxrs_parameters['Slit width'])
                    else:
                        self.dslit = 150
                
                elif float(self.cxrs_parameters['Central wavelength']) == 568:
                    self.grating=self.cxrs_parameters['Grating']
                    self.wavelength_setting=float(self.cxrs_parameters['Central wavelength'])
                    self.wstart = 568
                    self.wstop = 569
                    self.backg = [560,561]
                    self.fittype = "??"
                    self.mu = 568.5
                    self.kbt = 100
                    self.A = 1e-03
                    self.itern = 5
                    if float(self.cxrs_parameters['Slit width']) != -1:
                        self.dslit = float(self.cxrs_parameters['Slit width'])
                    else:
                        self.dslit = 150
                        
                elif float(self.cxrs_parameters['Central wavelength']) == 589:
                    self.grating=self.cxrs_parameters['Grating']
                    self.wavelength_setting=float(self.cxrs_parameters['Central wavelength'])
                    self.wstart = 588
                    self.wstop = 590
                    self.backg = [580,581]
                    self.fittype = "??"
                    self.mu = 589.5
                    self.kbt = 100
                    self.A = 1e-03
                    self.itern = 5
                    if float(self.cxrs_parameters['Slit width']) != -1:
                        self.dslit = float(self.cxrs_parameters['Slit width'])
                    else:
                        self.dslit = 150
                else: 
                    print(f'Central wavelength in the log ({self.cxrs_parameters['Central wavelength']}) does not match any of the predefined ones.')
                    print("Let's assume that it is CVI...")
                    self.grating=self.cxrs_parameters['Grating']
                    self.wavelength_setting=float(self.cxrs_parameters['Central wavelength'])
                    self.wstart = 528.25
                    self.wstop = 529.75
                    self.backg = [520,521]
                    self.fittype = "CVI"
                    self.mu = 528.93
                    self.kbt = 50
                    self.A = 1.5e-03
                    self.itern = 10
                    if float(self.cxrs_parameters['Slit width']) != -1:
                        self.dslit = float(self.cxrs_parameters['Slit width'])
                    else:
                        self.dslit = 150
            except:
                print("No logbook we use some defaults...")
                print("Let's assume that it is CVI...")
                self.grating='1800g_per_mm'
                self.wavelength_setting=528
                self.wstart = 528.25
                self.wstop = 529.75
                self.backg = [520,521]
                self.fittype = "CVI"
                self.mu = 528.93
                self.kbt = 50
                self.A = 1.5e-03
                self.itern = 10
                self.dslit = 150
        find_plasma_and_beam_time_overlap()
        calculate_evaluation_time_ranges()
        # extract_cxrs_parameters()
        get_cxrs_parameters()
        parameter_definition()

    def extract_signal_with_time(self,time, data, roi=None, plot=False):
        """
        Extracts the signal by interpolating the background and subtracting it from the original data.
        
        Parameters:
            time (numpy array): Time points corresponding to the data values.
            data (numpy array): Data values where even indices are background, and odd indices are background + signal.
            roi (str): Region of interest to include in plot titles. Default is None.
            plot (bool): Whether to plot the results. Default is False.
        
        Returns:
            signal (numpy array): The extracted signal for valid frames.
            interpolated_background (numpy array): The interpolated background.
            valid_time (numpy array): Time points for frames with valid signal.
        """
        if len(time) != len(data):
            raise ValueError("The time and data arrays must have the same length.")
        
        # Indices for background and signal frames
        background_indices = np.arange(1, len(data), 2)  # Even indices (background only)
        signal_indices = np.arange(0, len(data), 2)      # Odd indices (signal + background)
        
        # Extract background values and times
        background_time = time[background_indices]
        background = data[background_indices]
        
        # Interpolate the background across the entire time range
        background_interp_func = interp1d(background_time, background, kind='linear', fill_value="extrapolate")
        interpolated_background = background_interp_func(time)
        
        # Subtract interpolated background from the original data to get the signal
        signal = data - interpolated_background
    
        # Retain only valid signal points (odd frames)
        valid_signal = signal[signal_indices]
        valid_time = time[signal_indices]
    
        # Plot the results if requested
        if plot:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            
            # Define ROI title addition
            roi_title = f" (ROI: {roi})" if roi else ""
            
            # Plot the original data and interpolated background
            axs[0].plot(time, data, 'o-', label="Original Data (Background + Signal)")
            axs[0].plot(time, interpolated_background, '--', label="Interpolated Background")
            axs[0].set_title(f"Original Data and Interpolated Background{roi_title}")
            axs[0].set_ylabel("Value")
            axs[0].legend()
            axs[0].grid()
    
            # Plot the extracted signal for valid points only
            axs[1].plot(valid_time, valid_signal, 'o-', label="Extracted Signal", color='orange')
            axs[1].set_title(f"Extracted Signal (Valid Frames Only){roi_title}")
            axs[1].set_xlabel("Time")
            axs[1].set_ylabel("Signal Value")
            axs[1].legend()
            axs[1].grid()
            
            plt.tight_layout()
            plt.show()
        
        return valid_signal, interpolated_background, valid_time

    def plot_lightprofile(self, integration_time=1.0):        
        """
        Perform spectral evaluation with boxcar integration, downsampling using midpoint time for intervals, and return time ranges.
        
        Parameters:
            integration_time (float): Time interval (in seconds) for integration to reduce noise.
        
        Returns:
            time_ranges (list): A list of arrays representing the time ranges used for integration for each ROI.
        """
        # Generate spectral data
        spec = flap_w7x_abes.spectra(self.exp_id, campaign=self.campaign, spatcal=True)
        
        # Get active ROI coordinates
        self.active_roi = spec.dataobj.coordinate('ROI')[0][0, :, 0]
        
        # Perform wavelength calibration
        spec.wavelength_calibration(man=True, grid=self.grating, wavelength_setting=self.wavelength_setting)
        
        # Initialize storage for pcolormesh plotting
        roi_values = []  # Stores ROI values
        time_values = None  # To store midpoint time points
        signal_matrix = []  # Stores signal for each ROI
        time_ranges_all = []  # Stores integration time ranges for all ROIs
        
        # Loop over each ROI
        for roi in self.active_roi:
            # Extract data for the current ROI
            d = spec.slice_by_wl(roi, [self.wstart, self.wstop], plotting=False)
            
            # Extract signal and background
            valid_signal, interpolated_background, valid_time = self.extract_signal_with_time(
                time=d.coordinate('Time')[0],
                data=d.data,
                roi=roi,
                plot=False  # Disable individual ROI plotting
            )
            
            # Exclude the first and last seconds from valid_time and valid_signal
            time_mask = (valid_time > valid_time[0] + 0) & (valid_time < valid_time[-1] - 0)
            valid_signal = valid_signal[time_mask]
            valid_time = valid_time[time_mask]
    
            # Calculate rolling window size in points
            dt = np.mean(np.diff(valid_time))  # Average time step
            window_size = int(np.round(integration_time / dt))  # Number of points for the integration time
    
            # Apply boxcar integration (rolling average)
            boxcar_signal = np.convolve(valid_signal, np.ones(window_size)/window_size, mode='valid')
    
            # Compute the integration time ranges
            start_times = valid_time[:-window_size+1]
            end_times = valid_time[window_size-1:]
            mid_times = (start_times + end_times) / 2  # Midpoints for plotting
    
            # Downsample the signal and collect time ranges
            downsample_indices = np.arange(0, len(mid_times), window_size)
            downsampled_signal = boxcar_signal[downsample_indices]
            downsampled_mid_times = mid_times[downsample_indices]
            integration_ranges = np.column_stack((start_times[downsample_indices], end_times[downsample_indices]))
            
            # Append results for pcolormesh plot and time ranges
            roi_values.append(roi)
            signal_matrix.append(downsampled_signal)
            time_ranges_all.append(integration_ranges)
            
            # Store time values (only once, as they should be consistent across ROIs)
            if time_values is None:
                time_values = downsampled_mid_times
    
        # Convert signal_matrix to a NumPy array for plotting
        signal_matrix = np.array(signal_matrix)
    
        # Determine the color normalization
        vmin = np.min(signal_matrix)
        vmax = np.max(signal_matrix)
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)  # Center the colormap at 0
    
        # Create the pcolormesh plot
        plt.figure(figsize=(10, 6))
        mesh = plt.pcolormesh(
            time_values,
            roi_values,
            signal_matrix,
            cmap='seismic',
            norm=norm,
            shading='auto'  # Ensures no gaps between grid cells
        )
        plt.colorbar(mesh, label='Signal Intensity')
        plt.title('Signal Intensity Map (Boxcar Integration)')
        plt.xlabel('Time (s)')
        plt.ylabel('ROI')
        plt.show()
    
        # Optionally store the results in the class for further analysis
        self.signal_matrix = np.array(signal_matrix)
        self.time_values = np.array(time_values)
        self.time_ranges_all = np.array(time_ranges_all)
        
    def calc_temperature_profile(self, threshold=3.0):
        """
        Calculates the temperature profile by fitting the data for each ROI and each time interval 
        where the signal exceeds the specified threshold, and plots a contour plot of the results.
        
        Parameters:
            threshold (float): Minimum signal intensity required to run the analysis.
        """
        # Ensure the required properties are available
        if not hasattr(self, 'signal_matrix') or not hasattr(self, 'time_ranges_all'):
            raise AttributeError("Required properties (signal_matrix or time_ranges_all) are missing. "
                                 "Run plot_lightprofile first.")
    
        spec = flap_w7x_abes.spectra(self.exp_id, campaign=self.campaign, spatcal=True)
        spec.wavelength_calibration(man=True, grid=self.grating, wavelength_setting=self.wavelength_setting)
        # Initialize storage for temperature and error matrices, and full solutions
        temp_matrix = np.full(self.signal_matrix.shape, np.nan)  # NaN for skipped channels
        error_matrix = np.full(self.signal_matrix.shape, np.nan)
        full_solutions = [[None for _ in range(self.signal_matrix.shape[1])] for _ in range(self.signal_matrix.shape[0])]
    
        # Loop over each ROI and its corresponding time ranges
        for roi_idx, (roi, intervals) in enumerate(zip(self.active_roi, self.time_ranges_all)):
            for interval_idx, (tstart, tstop) in enumerate(intervals):
                # Check if the signal exceeds the threshold
                if self.signal_matrix[roi_idx, interval_idx] > threshold:
                    print(f"Processing ROI: {roi}, Interval: {tstart} - {tstop}, Signal: {self.signal_matrix[roi_idx, interval_idx]}")
    
                    # Perform temperature fitting
                    try:
                        solution, Terr = spec.tempfit(
                            self.fittype, roi, self.wstart, self.wstop,
                            self.mu, self.kbt, self.A, self.dslit,
                            tstart, tstop, self.backg, self.itern,
                            plots=False, calc_err=False
                        )
                        # Store full solution
                        full_solutions[roi_idx][interval_idx] = solution
    
                        # Extract temperature and error
                        temp_matrix[roi_idx, interval_idx] = solution.x[1]
                        error_matrix[roi_idx, interval_idx] = Terr
                    except Exception as e:
                        print(f"Error fitting ROI: {roi}, Interval: {tstart} - {tstop}: {e}")
                        temp_matrix[roi_idx, interval_idx] = np.nan
                        error_matrix[roi_idx, interval_idx] = np.nan
                        full_solutions[roi_idx][interval_idx] = None
                else:
                    print(f"Skipping ROI: {roi}, Interval: {tstart} - {tstop}, Signal below threshold.")
    
        # Store the results for further analysis
        self.temp_matrix = temp_matrix
        self.error_matrix = error_matrix
        self.full_solutions = full_solutions
    
        # Generate a contour plot of the temperature profile
        plt.figure(figsize=(10, 6))
        mesh = plt.pcolormesh(
            self.time_values,
            np.arange(len(self.active_roi)),
            temp_matrix,
            cmap='hot',
            shading='auto'
        )
        plt.colorbar(mesh, label='Temperature (arbitrary units)')
        plt.title('Temperature Profile')
        plt.xlabel('Time (s)')
        plt.ylabel('ROI Index')
        plt.show()
    
        # Optionally generate a contour plot for the error matrix
        plt.figure(figsize=(10, 6))
        mesh = plt.pcolormesh(
            self.time_values,
            np.arange(len(self.active_roi)),
            error_matrix,
            cmap='cool',
            shading='auto'
        )
        plt.colorbar(mesh, label='Temperature Error')
        plt.title('Temperature Error Profile')
        plt.xlabel('Time (s)')
        plt.ylabel('ROI Index')
        plt.show()
    # def calc_temperature_profile(self):
        # for roi in self.active_roi:
            # for tstart, tstop in self.intervals:
                # print(tstart,tstop)
                # spec.slice_by_wl(roi,[self.wstart, self.wstop])
                # spec.tempfit(self.fittype,roi,self.wstart,self.wstop,self.mu,self.kbt,self.A,self.dslit,self.tstart,self.tstop,self.backg,self.itern,plots=True)
                
# def check_line_modulation(self):
#     tstart = 2
#     tstop = 3
#     backg = [520,521]
#     spec = flap_w7x_abes.spectra(self.exp_id,campaign="OP2.2",spatcal=True)
#     self.active_roi=spec.dataobj.coordinate('ROI')[0][0,:,0]
#     spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
#     for roi in self.active_roi:extract_signal_with_time(time, data, plot=False)
#         spec.autocorr(roi,tstart,tstop,527,530)
#         # spec.sliceconvert_nanosecs_by_wl(roi,[528.5, 529.5])
#         # spec.active_passive(roi, tstart, tstop,background_interval=backg,error=True)

# def test_tempfit_op22():
#     exp_id = '20240926.031'
#     roi = 34
#     tstart = 2
#     tstop = 3
#     wstart = 528.25 #CVI line
#     wstop = 529.75
#     backg = [520,521]
#     fittype = "CVI"
#     mu = 528.93
#     kbt = 50
#     A = 1.5e-03
#     itern = 10
#     dslit = 150
#     spec = flap_w7x_abes.spectra(exp_id,campaign="OP2.2", spatcal=True)
#     spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    
#     spec.tempfit(fittype,roi,wstart,wstop,mu,kbt,A,dslit,
#                   tstart,tstop,backg,itern,plots=False)
    
plt.close("all")    
a=automated_cxrs_calculation("20240926.031")
# test_tempfit_op22()