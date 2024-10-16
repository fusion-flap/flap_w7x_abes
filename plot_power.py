# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:02:34 2024

@author: Zoletnik
"""

import flap
import flap_w7x_abes

flap_w7x_abes.register()

def plot_power(exp_ID,timerange=None,signals=['ABES-10','ABES-15','ABES-19', 'ABES-23'],
               datapath=None,background=False,fastchop=True,
               frange=None,fres=None,flog=None):
    
    if (fastchop):
        d = flap_w7x_abes.chopped_signals(exp_ID,timerange=timerange,signals=signals,datapath=datapath,background=background)
        options = {}
        if (fres is not None):
            options['Resolution'] = fres
        if (frange is not None):
            options['Range'] = frange
        if (flog is not None):
            options['Logarithmic'] = flog
            
        p = d.apsd(coordinate='Time',options=options)
        options = {'Log x':True,'Log y':True}
        p.plot(axes='Frequency',options=options)
    else:
        raise NotImplementedError("SLow chopping power plot is not implemented yet.")