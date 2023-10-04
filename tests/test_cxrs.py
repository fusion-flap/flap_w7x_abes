#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:10:03 2020

@author: mvecsei
"""

import matplotlib.pyplot as plt
import os

import flap
import flap_w7x_abes


if __name__ == '__main__':

    # shotID = '20230314.032'
    shotID = '20240314.032'

    # flap_w7x_abes.register()
    # d=flap.get_data('W7X_ABES_CXRS',exp_id=shotID, name="QSI")
    
    config = flap_w7x_abes.cxrs_util.read_fibre_config(exp_id=shotID)
    flap_w7x_abes.cxrs_util.plot_fibre_config(config)
    