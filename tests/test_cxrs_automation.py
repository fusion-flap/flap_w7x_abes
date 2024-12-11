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


flap_w7x_abes.register()
try:
    import flap_w7x_webapi as webapi
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
except ModuleNotFoundError:
    warnings.warn("flap_w7x_webapi not found, unable to download cxrs data")
   
class automated_cxrs_calculation():
    def __init__(self, exp_id, verbose=True,Wdia_threshold=200):       
        if (exp_id is None):
            raise ValueError('exp_id should be set for W7X ABES.')
        self.exp_id=exp_id
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
        self.spectral_evaluation()

    def check_if_shot_is_ok_for_evaluation(self):
        # check if there was plasma
        print(f"Start preliminary analysis for shot #"+self.exp_id)
        if flap.get_data('W7X_WEBAPI', name='Wdia', exp_id=self.exp_id, options={'Scale Time': True,'Cache Data': True}).data.max() > self.Wdia_threshold:
            self.plasma_ok=True
            print(f"Plasma discharge looks ok based on Wdia (>{self.Wdia_threshold}kA)")
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
                self.tstart=Wdia.get_coordinate_object('Time').values[crossings[0][0]]
                tstop=Wdia.get_coordinate_object('Time').values[crossings[1][0]]
                max_tstop=32 # maximum length of the beam injection
                self.tstop=tstop if tstop < max_tstop else max_tstop
                print(f"Calculation for the time range: {self.tstart:.1f}-{self.tstop:.1f} s.")  
            else:
                print(f"There is a problem with this shot, Wdia has crossed {self.Wdia_threshold}kA several times...")
                self.tstart=0
                self.tstop=10
                print(f"Using default {self.tstart}-{self.tstop} s for evaluation.")
                
        
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
        
        def extract_cxrs_parameters(log_data=self.logbook):
            """
            Extract CXRS parameters from a given log data dictionary.
        
            Args:
                log_data (dict): The input dictionary containing the log data.
        
            Returns:
                dict: A dictionary with extracted CXRS parameters.
            """
            # Extract the log description from the dictionary
            log_text = log_data.get('_source', {}).get('description', '')
            if not log_text:
                return {"Error": "Log description not found."}
        
            # Define a pattern to identify the CXRS section
            cxrs_pattern = r"---CXRS AUTOLOG---(.*?)---"
        
            # Search for the CXRS section in the log
            cxrs_section = re.search(cxrs_pattern, log_text, re.DOTALL)
            if not cxrs_section:
                return {"Error": "CXRS section not found in the log."}
        
            # Extracted section text
            cxrs_text = cxrs_section.group(1)
        
            # Define patterns to extract specific parameters
            patterns = {
                "Number of channels": r"Number of channels:\s*(\d+)",
                "Number of images": r"Number of images:\s*(\d+)",
                "Central wavelength": r"Central wavelength:\s*([\d.]+)",
                "Grating": r"Grating:\s*(\d+g_per_mm)",
                "Slit width": r"Slit width:\s*([\d.-]+)",
                "Major radii": r"Major radii\s*([\d.-]+-[\d.-]+)m",
                "Remark": r"Remark:\s*(.*)"
            }
        
            # Extract parameters using the defined patterns
            cxrs_parameters = {}
            for key, pattern in patterns.items():
                match = re.search(pattern, cxrs_text)
                if match:
                    cxrs_parameters[key] = match.group(1)
                else:
                    cxrs_parameters[key] = "Not found"
        
            self.cxrs_parameters=cxrs_parameters
        
        def parameter_definition():
            self.grating=self.cxrs_parameters['Grating']
            if float(self.cxrs_parameters['Central wavelength']) == 529:
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
            else: 
                print(f'Central wavelength in the log ({self.cxrs_parameters['Central wavelength']}) does not match any of the predefined ones.')
        
        find_plasma_and_beam_time_overlap()
        extract_cxrs_parameters()
        parameter_definition()
        
    def spectral_evaluation(self):        
        spec = flap_w7x_abes.spectra(self.exp_id,campaign=self.campaign,spatcal=True)
        self.active_roi=spec.dataobj.coordinate('ROI')[0][0,:,0]
        spec.wavelength_calibration(man=True,grid=self.grating,wavelength_setting=self.wavelength_setting)
        for roi in self.active_roi:
            spec.tempfit(self.fittype,roi,self.wstart,self.wstop,self.mu,self.kbt,self.A,self.dslit,self.tstart,self.tstop,self.backg,self.itern,plots=True)
                
# def check_line_modulation(self):
#     tstart = 2
#     tstop = 3
#     backg = [520,521]
#     spec = flap_w7x_abes.spectra(self.exp_id,campaign="OP2.2",spatcal=True)
#     self.active_roi=spec.dataobj.coordinate('ROI')[0][0,:,0]
#     spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
#     for roi in self.active_roi:
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
    
# plt.close("all")    
a=automated_cxrs_calculation("20241210.060")
# test_tempfit_op22()