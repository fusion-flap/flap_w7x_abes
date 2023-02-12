#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:53:41 2023

@author: refydani
"""

from tdms_data_reading import read_tdms_data
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import datetime

def nearest_ind(items, pivot):
    time_diff = np.abs([date - pivot for date in items])
    return time_diff.argmin(0)

class W7X_ABES_diagnostic():
    def __init__(self,shot=None):
        self.shot=shot     
        self.calc_resistor_chain()
        self.print_resistor_chain()
        # self.plot_resistor_chain()
        self.calc_emitter_current()
        self.calc_space_charge_limit()
        self.print_space_charge_limitt()
        self.print_emitter_current()
        self.plot_emitter_current()
    def print_space_charge_limitt(self):
        print('Ion current: {:3.2f} mA'.format(self.mean_space_charge_limit))
    def calc_space_charge_limit(self):
        ε0=const.epsilon_0
        e=const.e
        m=22.989769*const.atomic_mass
        d=2 # emitter-extractor distance[mm]
        K=(4/9)*ε0*(2*e/m)**(0.5)
        d_em=40 #emitter diameter [mm] - this is the scaling parameter, this should be measured with a series
        Uem=read_tdms_data('HV Em Meas Voltage',shot=self.shot,group_name='Beam')
        Uex=read_tdms_data('HV Ex Meas Voltage',shot=self.shot,group_name='Beam')
        du=Uem['data']-Uex['data']
        du[du<0] = 0 # suppress negative values for power calculation to avoid warnings
        self.space_charge_limit=K*(np.sqrt((du)*1e3)**3)/(d**2)*(d_em/2)**2*np.pi
        self.mean_space_charge_limit=np.mean(self.space_charge_limit)
    def print_emitter_current(self):
        print('Ion current: {:3.2f} mA'.format(self.ion_current))
    def calc_emitter_current(self):
        Uem=read_tdms_data('HV Em Meas Voltage',shot=self.shot,group_name='Beam')
        Uex=read_tdms_data('HV Ex Meas Voltage',shot=self.shot,group_name='Beam')
        Iem=read_tdms_data('HV Em Meas Current',shot=self.shot,group_name='Beam')
        Iex=read_tdms_data('HV Ex Meas Current',shot=self.shot,group_name='Beam')
        self.Uem=Uem['data']
        self.Uex=Uex['data']
        Iem_res=Uem['data']/self.em_chn_res
        Iex_res=Uex['data']/self.ex_chn_res
        self.t=(Iem['time']-Iem['time'][0])/ np.timedelta64(1, 's') #this converts us to seconds
        self.Iem=Iem['data']-Iem_res
        self.Iex=Iex['data']-Iex_res
        self.eta=1
        Iexs=self.Iex
        Iexs[Iexs==0] = 1e-6 # suppress zero values for power calculation to avoid warnings
        self.sigma=-1/((1+self.eta)*((self.Iem/Iexs)+(self.eta/(1+self.eta))))
        self.Iion=np.divide(self.Iem, (1+self.eta*self.sigma), out=np.zeros_like(self.Iem), where=(1+self.eta*self.sigma)!=0)
        ind=nearest_ind(Iem['time'], Iem['time'][0]+np.timedelta64(1,'s')) #cacluate only after beam stablizes
        self.ion_current=np.mean(self.Iion[ind:])
    def plot_emitter_current(self):
        fig, ax = plt.subplots(4,sharex=True, sharey=False)
        fig.set_figheight(15)
        fig.suptitle("W7-X ABES beamgun analysis\n shot#"+self.shot)
        ax[0].plot(self.t,self.Uem,color='orange',label="emitter")
        ax[0].plot(self.t,self.Uex,color='blue',label="extractor")
        ax[0].set_title('Voltages')
        ax[0].set_ylabel('kV')
        ax[0].legend()
        
        ax[1].plot(self.t,self.Iem,color='orange',label="emitter")
        ax[1].plot(self.t,self.Iex,color='blue',label="extractor")
        ax[1].axhline(y=0,linewidth=1,color='gray',linestyle='--')
        ax[1].set_title('Currents corrected for resistor chain')
        ax[1].set_ylabel('mA')
        ax[1].legend()
        
        ax[2].set_title('Beam ratio hitting the extractor ($\eta$={:3.1f})'.format(self.eta))
        ax[2].plot(self.t,self.sigma*100,color='black',label="extractor ion current ratio")
        ax[2].set_ylabel('%')
        ax[2].set_ylim([0,50])
        ax[2].legend() 
        
        ax[3].set_title('Ion beam current ($\eta$={:3.1f})'.format(self.eta))
        ax[3].plot(self.t,self.Iion,color='green',label="extracted ion current: {:3.2f} mA".format(self.ion_current))
        ax[3].axhline(y=self.ion_current,linewidth=1,color='green',linestyle='--')
        # ax[3].annotate('{:3.2f} mA'.format(self.ion_current), xy=[self.t[0]-0.1,self.ion_current],color='green')
        ax[3].plot(self.t,self.space_charge_limit,color='purple',label="space charge limit: {:3.2f} mA".format(self.mean_space_charge_limit),linestyle='--')
        ax[3].axhline(y=self.mean_space_charge_limit,linewidth=1,color='purple',linestyle='--')
        ax[3].set_ylabel('mA')
        ax[3].set_xlabel('time [s]')
        # ax[2].set_ylim([0,20])
        ax[3].legend() 
        # plt.tight_layout()
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
        Iemb=read_tdms_data('HV Em Meas Current',shot=self.shot,group_name='Beam')
        t0=Iemb['time'][0]
        Iem=read_tdms_data('HV Em Meas Current',shot=self.shot,group_name='Raise')
        Uem=read_tdms_data('HV Em Meas Voltage',shot=self.shot,group_name='Raise')
        Iex=read_tdms_data('HV Ex Meas Current',shot=self.shot,group_name='Raise')
        Uex=read_tdms_data('HV Ex Meas Voltage',shot=self.shot,group_name='Raise')
        ind=[nearest_ind(Iem['time'], t0-np.timedelta64(1,'s')),nearest_ind(Iem['time'], t0)]
        self.t_r_em=(Iem['time']-Iem['time'][0])/ np.timedelta64(1, 's')
        self.r_em = np.divide(Uem['data'], Iem['data'], out=np.zeros_like(Uem['data']), where=Iem['data']!=0)
        self.r_ex = np.divide(Uex['data'], Iex['data'], out=np.zeros_like(Uem['data']), where=Iex['data']!=0)
        self.em_chn_res=np.mean(self.r_em[ind[0]:ind[1]])
        self.ex_chn_res=np.mean(self.r_ex[ind[0]:ind[1]])

if __name__ == '__main__':  
        # shot='T20230207.002'
        shot='T20230210.003'
        plt.close('all')
        a=W7X_ABES_diagnostic(shot=shot)
        # emitter_current(shot=shot)
        # resistor_chain(shot=shot)
        # plt.plot(time,data)
        # plt.ylabel(unit)
        # plt.xlabel('time')
        # plt.show()
        