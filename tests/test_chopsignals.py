# -*- coding: utf-8 -*-
"""
Created on Fri May 10 18:50:28 2019

@author: Sandor Zoletnik  (zoletnik.sandor@wigner.mta.hu)
"""

import os

import flap
import flap_w7x_abes

flap_w7x_abes.register()

print("This program is obsolete!!!")

# Reading configuration file in the test directory
# thisdir = os.path.dirname(os.path.realpath(__file__))
# fn = os.path.join(thisdir,"w7x_config.cfg")
# flap.config.read(file_name=fn)

# flap_w7x_abes.proc_chopsignals(exp_id='20181017.024',
#                                timerange=[2,2.0005],
#                                signals='ABES-20',
#                                on_options={'Start':2.5,'End':0}, 
#                                off_options={'Start':2.5,'End':0},
#                                test=True)
# flap.list_data_objects()
# flap_w7x_abes.proc_chopsignals(exp_id='20181017.024',
#                                timerange=[2,2.1],
#                                on_options={'Start':2.5,'End':0}, 
#                                off_options={'Start':2.5,'End':0},
#                                test=False)
# flap.list_data_objects()
   