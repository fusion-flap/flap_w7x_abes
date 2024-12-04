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

plt.close("all")
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
        self.verbose=verbose
        self.Wdia_threshold=Wdia_threshold
        self.problem_override=True
        if self.check_if_shot_is_ok_for_evaluation() == False:
           warnings.warn(f"No sufficient data for Experiment ID: {self.exp_id}")

    def check_if_shot_is_ok_for_evaluation(self):
        # check if there was plasma
        if flap.get_data('W7X_WEBAPI', name='Wdia', exp_id=self.exp_id).data.max() > self.Wdia_threshold:
            self.plasma_ok=True
            print("Plasma ok")
        else:
            self.plasma_ok=False
            warnings.warn(f"Wdia is below the threshold: {self.Wdia_threshold} kJ (Experiment ID: {self.exp_id})")
            return False   
        
        if self.get_log() == True:
            # check if there was good ABES data and CXRS 
            if self.check_verdicts(self.logbook) == True:
                self.diagnostics_ok=True
            else:
                self.diagnostics_ok=False
                warnings.warn("Problems with the diagnostics")
                return False
            
            # check if shot was slow modulation
            if self.check_chopper_mode(self.logbook) == True:
                self.chopper_ok=True
            else:
                self.chopper_ok=False
                warnings.warn("Problems with the chopping mode")
                return False
        else:
            print("No logbook entry for the shot, lets assume that it was okay...")
            self.diagnostics_ok=None
            self.chopper_ok=None
        
        self.check_line_modulation()
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
                    chopper_mode = chopper_line.split(':')[-1].strip()
                    if chopper_mode == "Camera":
                        return True
                    elif chopper_mode == "Timed":
                        print(f"Chopper mode is 'Timed'. This is not ok for CXRS evaluation.")  # Debugging output
                        return False
                    else:
                        print(f"Chopper mode is unknown: {chopper_mode}")  # Debugging output
                        return False
                    
        print(f"Chopper mode is unknown: {chopper_mode}")
        return False

    def check_line_modulation(self):
        tstart = 2
        tstop = 3
        backg = [520,521]
        spec = flap_w7x_abes.spectra(self.exp_id,campaign="OP2.2",spatcal=True)
        self.active_roi=spec.dataobj.coordinate('ROI')[0][0,:,0]
        spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
        for roi in self.active_roi:
            spec.autocorr(roi,tstart,tstop,527,530)
            # spec.slice_by_wl(roi,[528.5, 529.5])
            # spec.active_passive(roi, tstart, tstop,background_interval=backg,error=True)

    def collect_parameters_for_spectral_evaluation():
        print("fasz")
        # start time of beam & plasma
        
        # end time of beam & plasma
        
        # wavelength of interest start, stop
        
        # wavelength range with no background
        
        # mu : initial guess for a term that constitutes from the Doppler shift of the line, the wavelength uncertainty of the spectrometer (float), and the line central wavelength
        
        # kbt : initial guess for temperature times Boltzmann constant (float)
        
        # A : initial guess for line intensity factor (float)
        
        # dslit : slit width in the spectrometer in the discharge (int)
        
        # central wavelength setting of the spectormeter
        
        # grid setting of the spectrometer
        
        # campaign name
        
        # fittype : (string) The kind of line that is to be generated. At the moment, it can be "CV" or "CVI".
        
        # N : Number of iterations for the Monte Carlo error calculation process (int)
    

def test_tempfit_op22():
    exp_id = '20240926.031'
    roi = 34
    tstart = 2
    tstop = 3
    wstart = 528.25 #CVI line
    wstop = 529.75
    backg = [520,521]
    fittype = "CVI"
    mu = 528.93
    kbt = 50
    A = 1.5e-03
    itern = 10
    dslit = 150
    spec = flap_w7x_abes.spectra(exp_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    
    spec.tempfit(fittype,roi,wstart,wstop,mu,kbt,A,dslit,
                  tstart,tstop,backg,itern,plots=False)
    
a=automated_cxrs_calculation("20240926.031")
# test_tempfit_op22()