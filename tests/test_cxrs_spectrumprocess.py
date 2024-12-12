# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:39:51 2024

@author: bcsillag
"""

import warnings
import matplotlib.pyplot as plt
import flap
import flap_w7x_abes

# plt.close("all")
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
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    flap.list_data_objects(spec.dataobj)
    
def test_calibration():
    expe_id = '20230316.072' 
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration(man=True,grid="1200g_per_mm",wavelength_setting=530)
    flap.list_data_objects(spec.dataobj)
    
    spec = flap_w7x_abes.spectra(expe_id,catest_slice_wavelength_op22_2mpaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    flap.list_data_objects(spec.dataobj)
    
def test_slice_wavelength():
    expe_id = '20230316.072' 
    roi = 4
    w = 529
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.slice_by_wl(roi,w)
    
def test_passive():
    expe_id = '20230314.063' 
    roi = 4
    tstart = 1
    tstop = 5
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.passive(roi, tstart, tstop)
    
def test_slice_wavelength_range():
    expe_id = '20230314.063' 
    roi = 4
    wstart = 588.7
    wstop = 589
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.slice_by_wl(roi,[wstart, wstop])
    
def test_active_passive():
    expe_id = '20230330.035' 
    roi = 4
    tstart = 1
    tstop = 4
    backg = [665,666]
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
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
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
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
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
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
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
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
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.error_distr(tstart,tstop,wstart,wstop,
                            background_interval=backg,plotting=True,ROI=roi)
    

def test_tempfit():
    #test 1 (slow!)
    expe_id = '20230316.047' 
    roi = 2
    tstart = 1
    tstop = 10
    wstart = 494
    wstop = 495.25
    backg = [488.5,491.5]
    fittype = "CV"
    mu = 494.72
    kbt = 100
    A = 7e-04
    itern = 5
    dslit = 100
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                 spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    
    spec.tempfit(fittype,roi,wstart,wstop,mu,kbt,A,dslit,
                  tstart,tstop,backg,itern,plots=True)
    
    #test 2
    expe_id = '20230316.072' 
    roi = 4
    tstart = 6
    tstop = 12
    wstart = 528.25
    wstop = 529.65
    backg = [535,540]
    fittype = "CVI"
    mu = 528.85
    kbt = 160
    A = 6.49104713e-04
    itern = 10
    dslit = 100
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                  spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    
    spec.tempfit(fittype,roi,wstart,wstop,mu,kbt,A,dslit,
                  tstart,tstop,backg,itern,plots=False)
    
def test_error_simulation():
    expe_id = '20230316.072' 
    roi = 4
    tstart = 6
    tstop = 12
    wstart = 528.25
    wstop = 529.65
    backg = [535,540]
    fittype = "CVI"
    mu = 528.85
    kbt = 160
    A = 6.49104713e-04
    itern = 10
    dslit = 100
    
    sf = ((20/7)**2)/6 #future measurement, 1 fiber, 1s, CVI
    simd = 100
    simgrid = "1800g_per_mm"
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.1",
                                  spatcal=True, time_correction=True)
    spec.wavelength_calibration()
    spec.Ti_error_simulation(fittype,roi,wstart,wstop,mu,kbt,A,dslit,
                    tstart,tstop,backg,itern,simd,simgrid,sf,plots=True)
    
def test_read_webapi_op22():
    expe_id = '20240926.031' 
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2",
                                 spatcal=True, time_correction=True)
    flap.list_data_objects(spec.dataobj)
    
def test_calibration_op22():
    expe_id = '20240926.031' 
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    flap.list_data_objects(spec.dataobj)
    
def test_slice_wavelength_op22():
    expe_id = '20240926.031' 
    roi = 34
    w = 529
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    spec.slice_by_wl(roi,w)
    
def test_passive_op22():
    expe_id = '20240926.031'
    roi = 33
    tstart = 2
    tstop = 3
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    spec.passive(roi, tstart, tstop)
    
def test_slice_wavelength_range_op22():
    expe_id = '20240926.027'
    roi = 34
    wstart = 588
    wstop = 590
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal="edge")
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=589)
    spec.slice_by_wl(roi,[wstart, wstop])
    
def test_active_passive_op22():
    expe_id = '20240926.031'
    roi = 34
    tstart = 2
    tstop = 3
    backg = [520,521]
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    spec.active_passive(roi, tstart, tstop,background_interval=backg,error=True)
    
def test_autocorr_op22():
    expe_id = '20240926.027'
    roi = 34
    tstart = 0
    tstop = 100
    wstart = 588
    wstop = 590
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal="edge")
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=589)
    spec.autocorr(roi,tstart,tstop,wstart,wstop)
    
def test_tempfit_op22():
    expe_id = '20240926.031'
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
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=529)
    
    spec.tempfit(fittype,roi,wstart,wstop,mu,kbt,A,dslit,
                  tstart,tstop,backg,itern,plots=False)
    
def test_slice_wavelength_range_op22_2():
    expe_id = '20241212.033'
    roi = 20
    wstart = 588
    wstop = 590
    spec = flap_w7x_abes.spectra(expe_id,campaign="OP2.2", spatcal=True)
    spec.wavelength_calibration(man=True,grid="1800g_per_mm",wavelength_setting=589)
    spec.slice_by_wl(roi,[wstart, wstop])
# test_read_webapi()
# test_calibration()
# test_slice_wavelength()
# test_passive()
# test_slice_wavelength_range()
# test_active_passive()
# test_autocorr()
# test_temporal_shift()
# test_get_line_intensity()
# test_error_distribution()
# test_tempfit()
# test_error_simulation()

# test_read_webapi_op22()
# test_calibration_op22()
# test_slice_wavelength_op22()
# test_passive_op22()
# test_slice_wavelength_range_op22()
# test_active_passive_op22()
# test_autocorr_op22()
# test_tempfit_op22()
test_slice_wavelength_range_op22_2()