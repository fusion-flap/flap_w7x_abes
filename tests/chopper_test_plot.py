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
                                  start_shift=0,
                                  end_shift=-3
                                  )
plt.figure()
# For testing a camera chopping measurement
flap_w7x_abes.test_chopper_timing(exp_id='20181018.003', 
                                  timerange=[0,3.],
                                  signal='ABES-15',
                                  resample=1e4,
                                  x_axis='Time',
                                  start_shift=1000,
                                  end_shift=-1000
                                  )