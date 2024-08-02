# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 15:20:58 2024

@author: Zoletnik
"""

import flap
default_options = {'Datapath':'data'}
options = flap.config.merge_options(default_options, {}, section='APDCAM control GUI')
print(options)