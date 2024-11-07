#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 15:43:44 2024

@author: zoletnik
"""
import matplotlib.pyplot as plt

import flap
import flap_w7x_abes
flap_w7x_abes.register()

plt.close('all')

# For testing a fast chopping measurement
flap_w7x_abes.test_chopper_timing(exp_id='20230323.062', 
                                  timerange=[3,3.001],
                                  signal='ABES-15',
                                  resample=None,
                                  x_axis='Time',
                                  )
plt.figure()
# For testing a camera chopping measurement
flap_w7x_abes.test_chopper_timing(exp_id='20181018.003', 
                                  timerange=[0,3.],
                                  signal='ABES-15',
                                  resample=1e5,
                                  x_axis='Time',
                                  beam_on_start_delay=1e3,
                                  beam_on_end_delay=-1e3,
                                  beam_off_start_delay=1e3,
                                  beam_off_end_delay=-1e3,
                                  )