# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:39:51 2024

@author: bcsillag
"""

import warnings
import matplotlib.pyplot as plt
import flap
import flap_w7x_abes

plt.close("all")
flap_w7x_abes.register()
try:
    import flap_w7x_webapi as webapi
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
except ModuleNotFoundError:
    warnings.warn("flap_w7x_webapi not found, unable to download cxrs data")
    
def test_read_webapi():
    expe_id = '20230316.072' 
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    flap.list_data_objects(spec.dataobj)
    
def test_calibration():
    expe_id = '20230316.072' 
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration(man=True,grid="1200g_per_mm",wavelength_setting=530)
    flap.list_data_objects(spec.dataobj)
    
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    flap.list_data_objects(spec.dataobj)
    
def test_slice_wavelength():
    expe_id = '20230316.072' 
    roi = 4
    w = 529
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.slice_by_wl(roi,w)
    
def test_passive():
    expe_id = '20230314.063' 
    roi = 4
    tstart = 1
    tstop = 5
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.passive(roi, tstart, tstop)
    
def test_slice_wavelength_range():
    expe_id = '20230314.063' 
    roi = 4
    wstart = 588.7
    wstop = 589
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.slice_by_wl(roi,[wstart, wstop])
    
def test_active_passive():
    expe_id = '20230330.035' 
    roi = 4
    tstart = 1
    tstop = 4
    backg = [665,666]
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.active_passive(roi, tstart, tstop,background_interval=backg,error=True)
    
def test_autocorr():
    expe_id = '20230330.035'
    roi = 4
    tstart = 1
    tstop = 5
    wstart = 670.1
    wstop = 670.3
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.autocorr(roi,tstart,tstop,wstart,wstop)
    
def test_temporal_shift():
    expe_id = '20230314.063' 
    roi = 4
    tstart = 1
    tstop = 5
    wstart = 588.7
    wstop = 589
    backg = [592,594]
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.tshif(roi,tstart,tstop,wstart,wstop,100,background_interval=backg)
    
def test_get_line_intensity():
    expe_id = '20230330.035'
    roi = 4
    tstart = 1
    tstop = 5
    wstart = 669.8
    wstop = 670.6
    backg = [665,666]
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    
    print(spec.get_line_intensity(roi,tstart,tstop,wstart,wstop,
                            background_interval=backg,plotting=True))
    
def test_error_distribution():
    expe_id = '20230330.035'
    roi = 4
    tstart = 1
    tstop = 5
    wstart = 670
    wstop = 670.4
    backg = [665,666]
    spec = flap_w7x_abes.spectra('W7X_WEBAPI',expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.error_distr(tstart,tstop,wstart,wstop,
                            background_interval=backg,plotting=True,ROI=roi)
    
    
# test_read_webapi()
# test_calibration()
# test_slice_wavelength()
# test_passive()
# test_slice_wavelength_range()
# test_active_passive()
# test_autocorr()
# test_temporal_shift()
# test_get_line_intensity()
test_error_distribution()