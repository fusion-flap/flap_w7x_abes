#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:53:41 2023

@author: refydani
"""

from tdms_data_reading import read_tdms_data
import numpy as np
import matplotlib.pyplot as plt

class W7X_ABES_diagnostic():
    def __init__(self,shot=None):
        self.shot=shot     
        self.calc_resistor_chain()
        self.print_resistor_chain()
        self.plot_resistor_chain()
        self.calc_emitter_current()
        self.print_emitter_current()
        self.plot_emitter_current()
        
    def print_emitter_current(self):
        print('Ion current: {:3.1f} mA'.format(self.ion_current))
    def calc_emitter_current(self):
        Uem=read_tdms_data('HV Em Meas Voltage',shot=self.shot,group_name='Beam')
        Uex=read_tdms_data('HV Ex Meas Voltage',shot=self.shot,group_name='Beam')
        Iem=read_tdms_data('HV Em Meas Current',shot=self.shot,group_name='Beam')
        Iex=read_tdms_data('HV Ex Meas Current',shot=self.shot,group_name='Beam')
        Iem_res=Uem['data']/self.em_chn_res
        Iex_res=Uex['data']/self.ex_chn_res
        self.t=(Iem['time']-Iem['time'][0])/ np.timedelta64(1, 's') #this converts us to seconds
        self.Iem=Iem['data']-Iem_res
        self.Iex=Iex['data']-Iex_res
        self.eta=1
        self.sigma=-1/((1+self.eta)*((self.Iem/self.Iex)+(self.eta/(1+self.eta))))
        self.Iion=np.divide(self.Iem, (1+self.eta*self.sigma), out=np.zeros_like(self.Iem), where=(1+self.eta*self.sigma)!=0)
        self.ion_current=np.mean(self.Iion)
    def plot_emitter_current(self):
        fig, ax = plt.subplots(3,sharex=True, sharey=False)
        fig.suptitle("W7-X ABES beamgun analysis\n shot#"+self.shot)
        ax[0].plot(self.t,self.Iem,color='orange',label="emitter")
        ax[0].plot(self.t,self.Iex,color='blue',label="extractor")
        ax[0].axhline(y=0,linewidth=1,color='gray',linestyle='--')
        ax[0].set_title('Currents corrected for resistor chain')
        ax[0].set_ylabel('mA')
        ax[0].legend()
        ax[1].set_title('Beam ratio hitting the extractor ($\eta$={:3.1f})'.format(self.eta))
        ax[1].plot(self.t,self.sigma*100,color='black',label="ion current on extractor")
        ax[1].set_ylabel('%')
        ax[1].set_ylim([0,20])
        ax[1].legend() 
        ax[2].set_title('Ion beam current ($\eta$={:3.1f})'.format(self.eta))
        ax[2].plot(self.t,self.Iion,color='green',label="emitter")
        ax[2].set_ylabel('mA')
        ax[2].set_xlabel('time [s]')
        # ax[2].set_ylim([0,20])
        ax[2].legend(title='calculated from') 
        plt.tight_layout()
        fig.show()
    def plot_resistor_chain(self):
        plt.figure(0)
        plt.suptitle("W7-X ABES resistor chain resistance\n HV ramp-up \n shot#"+self.shot)
        plt.plot(self.t_r_em,self.r_em,color='orange',label="emitter chain")
        plt.plot(self.t_r_em,self.r_ex,color='blue',label="extractor chain")
        plt.axhline(y=self.em_chn_res,linewidth=1,color='orange',linestyle='--')
        plt.axhline(y=self.ex_chn_res,linewidth=1,color='blue',linestyle='--')
        plt.annotate('{:3.1f} MOhm'.format(self.em_chn_res), xy=[self.t_r_em[-22],self.em_chn_res+2],color='orange')
        plt.annotate('{:3.1f} MOhm'.format(self.ex_chn_res), xy=[self.t_r_em[-22],self.ex_chn_res+2],color='blue')
        plt.ylabel('MOhm')
        plt.xlabel('time [s]')
        plt.ylim([0,100])
        plt.legend() 
        plt.show()
        
    
    def print_resistor_chain(self):
        print('Emitter resistor chain: {:3.1f} MOhm'.format(self.em_chn_res))
        print('Extractor resistor chain: {:3.1f} MOhm'.format(self.ex_chn_res))

    def calc_resistor_chain(self):
        Iem=read_tdms_data('HV Em Meas Current',shot=self.shot,group_name='Raise')
        Uem=read_tdms_data('HV Em Meas Voltage',shot=self.shot,group_name='Raise')
        Iex=read_tdms_data('HV Ex Meas Current',shot=self.shot,group_name='Raise')
        Uex=read_tdms_data('HV Ex Meas Voltage',shot=self.shot,group_name='Raise')
        self.t_r_em=(Iem['time']-Iem['time'][0])/ np.timedelta64(1, 's')
        self.r_em = np.divide(Uem['data'], Iem['data'], out=np.zeros_like(Uem['data']), where=Iem['data']!=0)
        self.r_ex = np.divide(Uex['data'], Iex['data'], out=np.zeros_like(Uem['data']), where=Iex['data']!=0)
        self.em_chn_res=np.mean(self.r_em[-20:])
        self.ex_chn_res=np.mean(self.r_ex[-20:])


if __name__ == '__main__':  
        shot='T20230210.003'
        plt.close('all')
        W7X_ABES_diagnostic(shot=shot)
        # emitter_current(shot=shot)
        # resistor_chain(shot=shot)
        # plt.plot(time,data)
        # plt.ylabel(unit)
        # plt.xlabel('time')
        # plt.show()
        