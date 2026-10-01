#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 17:22:47 2026

@author: mive
"""

import matplotlib.pyplot as plt
import numpy as np
import os

import flap
from .w7x_abes import abes_get_config

class APDCAMPanel():
    
    def __init__(self):
        self.fiber_centers = {'1':[1294, 1181],
                         '2':[1534, 1147],
                         '3':[1400, 1046],
                         '4':[1137, 843],
                         '5':[285, 541],
                         '6':[123, 488],
                         '7':[586, 839],
                         '8':[415, 711],
                         '9':[201, 334],
                         '10':[366, 389],
                         '11':[477, 523],
                         '12':[590, 663],
                         '13':[1608, 995],
                         '14':[1300, 896],
                         '15':[1641, 827],
                         '16':[1468, 889],
                         '17':[1307, 680],
                         '18':[1135, 664],
                         '19':[1638, 656],
                         '20':[1471, 715],
                         '21':[861, 753],
                         '22':[486, 266],
                         '23':[625, 372],
                         '24':[745, 487],
                         '25':[242, 714],
                         '26':[86, 656],
                         '27':[256, 889],
                         '28':[82, 826],
                         '29':[1433, 1290],
                         '30':[1181, 1035],
                         '31':[1301, 1401],
                         '32':[1131, 1428],
                         '33':[293, 1290],
                         '34':[429, 1186],
                         '35':[525, 1051],
                         '36':[770, 1041],
                         '37':[975, 485],
                         '38':[1279, 92],
                         '39':[862, 362],
                         '40':[1081, 79],
                         '41':[1416, 201],
                         '42':[1241, 266],
                         '43':[1105, 374],
                         '44':[1362, 390],
                         '45':[647, 1178],
                         '46':[424, 1401],
                         '47':[778, 1426],
                         '48':[598, 1428],
                         '49':[863, 1215],
                         '50':[946, 1426],
                         '51':[1079, 1178],
                         '52':[965, 1034],
                         '53':[1524, 332],
                         '54':[1598, 488],
                         '55':[1440, 542],
                         '56':[1232, 498],
                         '57':[444, 90],
                         '58':[311, 198],
                         '59':[668, 76],
                         '60':[884, 81],
                         '61':[326, 1049], 
                         '62':[187, 1154],
                         '63':[424, 891],
                         '64':[119, 993],
                         }
        self.state_2021 = {'1': "APD",
                         '2': "APD",
                         '3': "APD",
                         '4': "MPPC",
                         '5': "APD",
                         '6': "APD",
                         '7': "APD",
                         '8': "APD",
                         '9': "BROKEN",
                         '10': "APD",
                         '11': "NOISY",
                         '12': "APD",
                         '13': "APD",
                         '14': "MPPC",
                         '15': "APD",
                         '16': "MPPC",
                         '17': "APD",
                         '18': "APD",
                         '19': "DISC.",
                         '20': "APD",
                         '21': "APD",
                         '22': "APD",
                         '23': "APD",
                         '24': "APD",
                         '25': "APD",
                         '26': "APD",
                         '27': "APD",
                         '28': "DISC.",
                         '29': "APD",
                         '30': "APD",
                         '31': "APD",
                         '32': "APD",
                         '33': "MPPC",
                         '34': "MPPC",
                         '35': "MPPC",
                         '36': "MPPC",
                         '37': "APD",
                         '38': "APD",
                         '39': "APD",
                         '40': "APD",
                         '41': "APD",
                         '42': "NOISY",
                         '43': "APD",
                         '44': "APD",
                         '45': "MPPC",
                         '46': "MPPC",
                         '47': "APD",
                         '48': "MPPC",
                         '49': "MPPC",
                         '50': "APD",
                         '51': "APD",
                         '52': "APD",
                         '53': "APD",
                         '54': "APD",
                         '55': "APD",
                         '56': "APD",
                         '57': "APD",
                         '58': "APD",
                         '59': "MPPC",
                         '60': "APD",
                         '61': "MPPC",
                         '62': "MPPC",
                         '63': "MPPC",
                         '64': "MPPC",
                         }
        
        self.state_2024 = {'1': "APD",
                         '2': "APD",
                         '3': "APD",
                         '4': "MPPC",
                         '5': "APD",
                         '6': "APD",
                         '7': "Amplifier Error",
                         '8': "APD",
                         '9': "APD",
                         '10': "APD",
                         '11': "APD",
                         '12': "APD",
                         '13': "APD",
                         '14': "MPPC",
                         '15': "APD",
                         '16': "MPPC",
                         '17': "Broken",
                         '18': "APD",
                         '19': "APD",
                         '20': "APD",
                         '21': "APD",
                         '22': "APD",
                         '23': "APD",
                         '24': "APD",
                         '25': "APD",
                         '26': "APD",
                         '27': "APD",
                         '28': "APD",
                         '29': "APD",
                         '30': "APD",
                         '31': "APD",
                         '32': "APD",
                         '33': "MPPC",
                         '34': "MPPC",
                         '35': "MPPC",
                         '36': "MPPC",
                         '37': "Occ. NOISY",
                         '38': "Bad signal 2022",
                         '39': "APD",
                         '40': "APD",
                         '41': "APD",
                         '42': "APD",
                         '43': "APD",
                         '44': "APD",
                         '45': "MPPC",
                         '46': "MPPC",
                         '47': "APD",
                         '48': "MPPC",
                         '49': "MPPC",
                         '50': "APD",
                         '51': "Noisy",
                         '52': "DISC",
                         '53': "APD",
                         '54': "APD",
                         '55': "APD",
                         '56': "APD",
                         '57': "APD",
                         '58': "No Signal",
                         '59': "MPPC",
                         '60': "Noisy",
                         '61': "MPPC",
                         '62': "MPPC",
                         '63': "MPPC",
                         '64': "MPPC",
                         }
        
        self.state = self.state_2024
        self.exp_id = None
        self.patching = None
        
    def get_exp_config(self, exp_id, options={}):        
        options_default = {
            'Datapath': flap.config.get("Module W7X_ABES", "Datapath")}
        options = {**options_default, **options}
        
        self.exp_id = exp_id
        if int(self.exp_id.split(".")[0])<20240101:
            self.state = self.state_2021

        # Get micrometer settings
        datapath_base = options['Datapath']
        datapath = os.path.join(datapath_base, self.exp_id)
        xmlfile = os.path.join(datapath, self.exp_id + '_config.xml')
        xml = flap.FlapXml()
        xml.read_file(xmlfile)
        try:
            if (xml.head.tag != 'ShotSettings'):
                raise ValueError(
                    "No ShotSettings entry found in XML file " + xmlfile + ".")
            if (xml.head.attrib['Experiment'] != "W7-X A-BES"):
                raise ValueError(xmlfile + " is not a W7-X ABES config file.")
        except Exception:
            raise ValueError("File format error in " + xmlfile)
        config = abes_get_config(xml)
        
        self.patching = dict()
        for index, adc in enumerate(config['ADC_list']):
            self.patching[str(adc)]= [config['signal_list'][index],
                                      config['fibre_list'][index]]
        
        return config

    def plot_channel_locations(self, info="ADC", exp_id=None):
        if exp_id != None:
            self.exp_id = exp_id
            self.get_exp_config(self.exp_id)
        if info != "ADC" and self.exp_id == None:
            raise ValueError("Without experiment id can only plot the ADC numbers")
        
        # Create the plot
        plt.figure(figsize=(12, 10))
        
        color_dict = {"APD": "tab:blue",
                      "MPPC": "tab:green",
                      "DEAD": "tab:red"}
        
        # Create circles at each location
        for channelid in self.fiber_centers.keys():
            x_coord = self.fiber_centers[channelid][0]
            y_coord = self.fiber_centers[channelid][1]
            channel_type = self.state[channelid]
            if channel_type in color_dict.keys():
                chcolor = color_dict[channel_type]
            else:
                chcolor = color_dict["DEAD"]
            circle = plt.Circle((x_coord, y_coord), radius=70, 
                               color=chcolor, alpha=0.5, linewidth=1)
            plt.gca().add_patch(circle)
            
            # Add the ChannelID text at the center
            text = channelid
            textbg = "black"
            textc = "white"
            if info == "Signal":
                if channelid in self.patching.keys():
                    text = self.patching[channelid][0]
                    textbg = "xkcd:light lilac"
                    textc = "black"
            elif info == "Fiber":
                if channelid in self.patching.keys():
                    text = self.patching[channelid][1]
                    textbg = "xkcd:light lilac"
                    textc = "black"

            plt.text(x_coord, y_coord, text, 
                    ha='center', va='center', fontsize=12, fontweight='bold',
                    color=textc, bbox=dict(boxstyle="round,pad=0.3", facecolor=textbg, alpha=0.7))
        
        circle = plt.Circle((40, 40), radius=60, 
                           color="tab:blue", alpha=0.5, linewidth=1)
        plt.gca().add_patch(circle)
        plt.text(40, 40, "APD", 
                ha='center', va='center', fontsize=12, fontweight='bold',
                color="white")
        circle = plt.Circle((40, 170), radius=60, 
                           color="tab:green", alpha=0.5, linewidth=1)
        plt.gca().add_patch(circle)
        plt.text(40, 170, "MPPC", 
                ha='center', va='center', fontsize=12, fontweight='bold',
                color="white")
        circle = plt.Circle((40, 300), radius=60, 
                           color="tab:red", alpha=0.5, linewidth=1)
        plt.gca().add_patch(circle)
        plt.text(40, 300, "N.A.", 
                ha='center', va='center', fontsize=12, fontweight='bold',
                color="white")
        
        if info != "ADC":
            plt.text(1600, 40, "ADC num. (free)", 
                    ha='center', va='center', fontsize=10, fontweight='bold',
                    color="white", bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.7))
            plt.text(1600, 80, f"{info} num. (connected)", 
                    ha='center', va='center', fontsize=10, fontweight='bold',
                    color=textc, bbox=dict(boxstyle="round,pad=0.3", facecolor=textbg, alpha=0.7))
        
        # Set axis properties
        plt.xlabel('')
        plt.ylabel('')
        plt.title(f'APDCAM backpanel\n{info} numbers')
        
        # Set equal scaling to ensure circles look circular
        plt.axis('equal')
        
        bgfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils", "panel_bg.png")
        plt.imshow(plt.imread(bgfile))
        plt.gca().invert_yaxis()
        
        plt.xticks([])
        plt.yticks([])
        plt.axis("off")
        plt.ylim([-50, 1550])
        
        # Show the plot
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.show()
