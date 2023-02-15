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

def smooth(y, box_pts):
    # box = np.ones(box_pts)/box_pts
    # y_smooth = np.convolve(y, box, mode='same')
    y_smooth=[]
    for i in range(0, len(y)):
        li = int(i - (box_pts-1)/2)
        if li < 0:
            li = 0
        ui = li + box_pts
        if ui > len(y)-1:
            ui = len(y)-1
            li = ui - box_pts
        y_smooth.append(np.mean(y[li:ui]))
    return np.array(y_smooth)

class W7X_ABES_diagnostic():
    def __init__(self,shot=None):
        self.shot=shot     
        self.plotpath='/data/W7-X/Beam/tdms_figures/'
        self.calc_resistor_chain()
        self.calc_emitter_current()
        self.calc_space_charge_limit()
        self.calc_toroidal_aiming_test()
        self.calc_fc2()
        
        self.print_shotnumber()
        self.print_HV()
        self.print_resistor_chain()       
        self.print_space_charge_limitt()
        self.print_emitter_current()
    
        # self.plot_resistor_chain()
        # self.plot_main_parameter()
        self.plot_emitter_current()
        # self.plot_focusing_test()
        # self.plot_toroidal_aiming_test()
        
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
        self.Uextractor=np.mean(Uex['data'][ind:])
        self.Uemitter=np.mean(Uem['data'][ind:])
        
    def calc_space_charge_limit(self):
        ε0=const.epsilon_0
        e=const.e
        m=22.989769*const.atomic_mass
        d=2 # emitter-extractor distance[mm]
        K=(4/9)*ε0*(2*e/m)**(0.5)
        d_em=60 #emitter diameter [mm] - this is the scaling parameter, this should be measured with a series
        Uem=read_tdms_data('HV Em Meas Voltage',shot=self.shot,group_name='Beam')
        Uex=read_tdms_data('HV Ex Meas Voltage',shot=self.shot,group_name='Beam')
        du=Uem['data']-Uex['data']
        du[du<0] = 0 # suppress negative values for power calculation to avoid warnings
        self.dU=du
        self.space_charge_limit=K*(np.sqrt((du)*1e3)**3)/(d**2)*(d_em/2)**2*np.pi
        ind=nearest_ind(Uem['time'], Uem['time'][0]+np.timedelta64(1,'s')) #cacluate only after beam stablizes
        self.mean_space_charge_limit=np.mean(self.space_charge_limit[ind:])
    
    def calc_fc2(self):
        fc=read_tdms_data('FC2 Resistor Current mA',shot=self.shot,group_name='Beam')
        I_fc_ps=read_tdms_data('- FC Current mA',shot=self.shot,group_name='Beam')
        ind=nearest_ind(fc['time'], fc['time'][0]+np.timedelta64(1,'s'))
        self.I_fc=fc['data']
        self.fc_curr=np.mean(fc['data'][ind:])
        self.I_fc_ps=I_fc_ps['data']
        
    def calc_toroidal_aiming_test(self):
        U_aim_tor=read_tdms_data('- Aiming Control (tor) V',shot=self.shot,group_name='Beam')
        self.U_aim_tor=U_aim_tor['data']
    
    def print_HV(self):
        print('Emitter voltage: {:3.1f} kV'.format(self.Uemitter))
        print('Extractor voltage: {:3.1f} kV'.format(self.Uextractor))
        print('Voltage difference: {:3.1f} kV'.format(abs(self.Uextractor-self.Uemitter) ))
    
    def print_shotnumber(self):
        print("W7-X ABES beamgun parameters, shot#"+self.shot)
        
    def print_resistor_chain(self):
        print('Emitter resistor chain: {:3.1f} MOhm'.format(self.em_chn_res))
        print('Extractor resistor chain: {:3.1f} MOhm'.format(self.ex_chn_res))
        
    def print_space_charge_limitt(self):
        print('Space charge limit: {:3.2f} mA'.format(self.mean_space_charge_limit))
        
    def print_emitter_current(self):
        print('Ion current: {:3.2f} mA'.format(self.ion_current))
        
    def plot_emitter_current(self):
        fig, ax = plt.subplots(5,sharex=True, sharey=False)
        fig.set_figheight(8)
        fig.suptitle("W7-X ABES beamgun analysis\n shot#"+self.shot)
        ax[0].plot(self.t,self.Uem,color='orange',label="emitter: {:3.1f} kV".format(self.Uemitter))
        ax[0].plot(self.t,self.Uex,color='blue',label="extractor: {:3.1f} kV".format(self.Uextractor))
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
        ax[2].set_ylim([0,20])
        ax[2].legend() 
        
        ax[3].set_title('Ion beam current ($\eta$={:3.1f})'.format(self.eta))
        ax[3].plot(self.t,self.Iion,color='green',label="extracted ion current: {:3.2f} mA".format(self.ion_current))
        ax[3].axhline(y=self.ion_current,linewidth=1,color='green',linestyle='--')
        # ax[3].annotate('{:3.2f} mA'.format(self.ion_current), xy=[self.t[0]-0.1,self.ion_current],color='green')
        ax[3].plot(self.t,self.space_charge_limit,color='purple',label="space charge limit: {:3.2f} mA".format(self.mean_space_charge_limit),linestyle='--')
        # ax[3].plot(self.t,smooth(self.fc,1),color='black',label="FC2 current: {:3.2f} mA".format(self.fc_curr),linestyle='-')
        
        ax[3].axhline(y=self.mean_space_charge_limit,linewidth=1,color='purple',linestyle='--')
        ax[3].set_ylabel('mA')
        ax[3].set_xlabel('time [s]')
        ax[3].legend() 
        
        ax[4].set_title('Ion beam current ($\eta$={:3.1f})'.format(self.eta))
        # ax[4].axhline(y=self.ion_current,linewidth=1,color='green',linestyle='--')
        ax[4].plot(self.t,smooth(self.I_fc,1),color='black',label="FC2 current",linestyle='-')
        # ax[4].plot(self.t,smooth(self.I_fc_ps,1),color='blue',label="FC2 house PS",linestyle='--')
        ax[4].plot(self.t,smooth(self.I_fc,5)-smooth(self.I_fc_ps,5),color='pink',label="FC2-FC2 house PS",linestyle='--')
        ax[4].plot(self.t,self.Iion,color='green',label="extracted ion current")
        ax[4].set_ylabel('mA')
        ax[4].set_xlabel('time [s]')
        ax[4].legend() 
        
        # plt.tight_layout()
        plt.savefig(self.plotpath+self.shot+'_processed_beam_properties.png')
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
        
    def plot_main_parameter(self):
        channels=['HV Em Meas Voltage','HV Em Meas Current','HV Ex Meas Voltage','HV Ex Meas Current','FC2 Resistor Current mA','- Aiming Control (tor) V']
        fig6, ax = plt.subplots(len(channels),sharex=True, sharey=False)
        fig6.suptitle("W7-X ABES beam parameters, shot#"+shot)
        for i in range(len(channels)):
            d=read_tdms_data(channels[i],shot=self.shot,group_name=None)
            ax[i].plot(d['time'],d['data'])
            ax[i].set_title(channels[i])
            ax[i].set_ylabel(d['unit'])
            if i == len(channels)-1:
                ax[i].set_xlabel("time")
        plt.tight_layout()
        plt.show()
        plt.savefig(self.plotpath+self.shot+'_main_parameter_plot.png')
        
    def plot_focusing_test(self):
        plt.figure(10)
        plt.suptitle("W7-X ABES focusing test shot#"+self.shot)
        plt.scatter(smooth(((self.Uem-self.Uex)[40:]),5),smooth(((self.fc)[40:]),5),color='black',label="FC current")
        plt.ylabel('mA')
        plt.xlabel('kV')
        plt.legend() 
        plt.show()
        
    def plot_toroidal_aiming_test(self):
        plt.figure(10)
        plt.suptitle("W7-X ABES toroidal aiming test shot#"+self.shot)
        plt.scatter(smooth((self.U_aim_tor[40:]),5),smooth(((self.fc)[40:]),5),color='black',label="FC current")
        plt.ylabel('mA')
        plt.xlabel('V')
        plt.legend() 
        plt.show()

    
if __name__ == '__main__':  
        shot='20230215.027'
        # shot='pina'
        # plt.close('all')
        a=W7X_ABES_diagnostic(shot=shot)
        # emitter_current(shot=shot)
        # resistor_chain(shot=shot)
        # plt.plot(time,data)
        # plt.ylabel(unit)
        # plt.xlabel('time')
        # plt.show()
        