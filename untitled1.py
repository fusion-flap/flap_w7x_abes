#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:37:07 2024

@author: apdcam
"""

import flap
import flap_w7x_abes

# flap_w7x_abes.plot_power("20241106.027",signals="ABES-*", plot_type="grid", signal_type='raw',
#                          beam_on=True, timerange=[0.5,1],
#                          crosspower_amplitude=True, resample=None)

# clean_abes = flap_w7x_abes.get_clean_abes("20241106.055",signals="ABES-20", signal_type="beam", timerange=[0.1,1.8])
# clean_abes.plot()

flap_w7x_abes.test_chopper_timing("20241106.055", signal="ABES-20", timerange=[3,3.001])
