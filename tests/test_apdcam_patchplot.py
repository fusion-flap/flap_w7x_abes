#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:10:03 2020

@author: mvecsei
"""

import flap
import flap_w7x_abes

if __name__ == '__main__':

    panel = flap_w7x_abes.APDCAMPanel()
    panel.plot_channel_locations(info="ADC")
    panel.plot_channel_locations(info="Signal", exp_id="20250226.038")
    panel.plot_channel_locations(info="Fiber")