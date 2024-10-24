# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 21:15:27 2021

@author: Zoletnik
"""

import matplotlib.pyplot as plt
import numpy as np
try:
    try:
        from .nptdms_mod import TdmsFile
    except ImportError:
        from nptdms_mod import TdmsFile
except ModuleNotFoundError:
    from flap_w7x_abes.nptdms_mod import TdmsFile
    
import datetime
import os
import copy
from scipy.interpolate import interp1d
import flap
import matplotlib.gridspec as gridspec


try:
    try:
        from .utc_offset import UTC_offset
    except ImportError:
        from utc_offset import UTC_offset
except ModuleNotFoundError:
    from flap_w7x_abes.utc_offset import UTC_offset

class BORIMonitor():
    def __init__(self, date=None, exp_id=None, datapath='/data',
                 material="Na"):
        if date is None and exp_id is None:
            raise ValueError("Either date or experiment id should be given to read a beam log")
        self.date = date
        self.exp_id = exp_id
        self.datapath = datapath
        
        self.read_tdms()
        self.get_vap_pressure(material=material)
        

    def read_tdms(self):
        data_names = {'Emit Current A':"Current:A",
                      'HV Em Meas Voltage':"Voltage:kV",
                      'HV Ex Meas Voltage':"Voltage:kV",
                      'HV Em Meas Current':"Current:mA",
                      'HV Ex Meas Current':"Current:mA",
                      'TC Oven Top':"Temperature:C",
                      'TC Oven Bottom':"Temperature:C",
                      'TC Torus Side Cone':"Temperature:C",
                      'TC Emit Side Cone':"Temperature:C",
                      'FC1 in':"n.a.:n.a.",
                      'FC2 in':"n.a.:n.a.",
                      'FC Polarity':"n.a.:n.a",
                      'FC1 Resistor Current mA':"Current:mA",
                      'FC2 Resistor Current mA':"Current:mA",
                      'VG HighVac1':"Pressure:mbar",
                      'VG HighVac2':"Pressure:mbar",
                      'VG ForeVac':"Pressure:mbar",
                      'Neut Shut Closed':"n.a.:n.a",
                      'PB Feedback-Open':"n.a.:n.a"}
        if self.date is not None:
            t,d,u = read_date_tdms(data_names=list(data_names.keys()),startdate=self.date,datapath=self.datapath)
            reftime = np.datetime64(self.date)
        elif self.exp_id is not None:
            t,d,u = read_exp_tdms(data_names=list(data_names.keys()),exp_id=self.exp_id,datapath=self.datapath)
            reftime = t[0]
        
        #normalizing the time and creating a flap dataobject
        t = t-reftime
        t=t/np.timedelta64(1,"s")
        t= t.astype("float")
        time_coord = flap.Coordinate(name='Time', unit='Second',values=t,
                                     shape=t.shape,
                                     mode=flap.CoordinateMode(equidistant=False),
                                     dimension_list=[0])
        
        self.data=dict()
        for index, data_name in enumerate(data_names):
            self.data[data_name] = flap.DataObject(data_array=d[index],
                                                   data_unit=flap.coordinate.Unit(
                                                       name=data_names[data_name].split(":")[0],
                                                       unit=data_names[data_name].split(":")[1]),
                                                   coordinates=[time_coord], exp_id=self.date,
                                                   data_title=data_name, data_shape=d[index].shape,
                                                   info="BORI log data")
        

    def slice_time(self, last_minutes=20):
        return_dataobject = copy.deepcopy(self)
        timevec = self.data['Emit Current A'].get_coordinate_object("Time").values
        xlim = [timevec[np.where(timevec > timevec[-1] - last_minutes*60)][0], timevec[-1]]
        for key in self.data.keys():
            return_dataobject.data[key] = return_dataobject.data[key].slice_data(slicing={"Time":flap.Intervals(xlim[0], xlim[1])})
        return return_dataobject

    def get_load_resistance(self):
        """
        Checks the locations where the Emitter voltage and the extractor voltage are stable, 
        where the extractor voltage is larger than the emitter voltage and where both of them are above 100V.
        It calculates the resistance from these data
        -------
        Returns the emitter and extractor resistance in MOhm

        """
        em_voltage_der = np.abs(self.data['HV Em Meas Voltage'].data[1:]-self.data['HV Em Meas Voltage'].data[:-1])
        ex_voltage_der = np.abs(self.data['HV Ex Meas Voltage'].data[1:]-self.data['HV Ex Meas Voltage'].data[:-1])
        
        rel_time = np.where((em_voltage_der<0.01)*\
                            (ex_voltage_der<0.01)*\
                            (self.data['HV Em Meas Voltage'].data[1:]>0.1)*\
                            (self.data['HV Ex Meas Voltage'].data[1:]>self.data['HV Em Meas Voltage'].data[1:]))
        
        self.em_resistance = np.nanmean(self.data['HV Em Meas Voltage'].data[rel_time]/self.data['HV Em Meas Current'].data[rel_time])
        self.ex_resistance = np.nanmean(self.data['HV Ex Meas Voltage'].data[rel_time]/self.data['HV Ex Meas Current'].data[rel_time])

        return self.em_resistance, self.ex_resistance
    
    def get_overcurrent(self):
        """
        Calculates the overcurrents on the extractor and the emitter and plots them
        """
        
        if (hasattr(self, "em_resistance") is False) or (hasattr(self, "ex_resistance") is False):
            self.get_load_resistance()

        extractor_overcurrent = copy.deepcopy(self.data['HV Ex Meas Current'])
        extractor_overcurrent.data = self.data['HV Ex Meas Current'].data-self.data['HV Ex Meas Voltage'].data/self.ex_resistance
        extractor_overcurrent.name = "Extractor overcurrent"
        self.data["Extractor overcurrent"] = extractor_overcurrent
        
        emitter_overcurrent = copy.deepcopy(self.data['HV Em Meas Current'])
        emitter_overcurrent.data = self.data['HV Em Meas Current'].data-self.data['HV Em Meas Voltage'].data/self.em_resistance
        emitter_overcurrent.name = "Emitter overcurrent"
        self.data["Emitter overcurrent"] = emitter_overcurrent

    def get_extracted_charge(self):
        if ("Emitter overcurrent" not in list(self.data.keys())) or ("Extractor overcurrent" not in list(self.data.keys())):
            self.get_overcurrent()
        timevect = self.data["Emitter overcurrent"].get_coordinate_object("Time").values
        rel_data = self.data["Emitter overcurrent"].data.clip(min=0)
        failed_data = np.where(np.isnan(self.data["Emitter overcurrent"].data))
        rel_data[failed_data] = 0
        
        extracted_charge = copy.deepcopy(self.data["Emitter overcurrent"])
        extracted_charge.data_title = "Extracted charge"
        extracted_charge.data_unit.name = "Charge"
        extracted_charge.data_unit.unit = "mC"
        extracted_charge.data[0] = 0
        for index in range(len(extracted_charge.data)-1):
            extracted_charge.data[index+1] = extracted_charge.data[index] + rel_data[index]*(timevect[index+1]-timevect[index])
        
        self.data["Extracted charge"] = extracted_charge
        

    def get_child_langmuir(self, neutralizer_shutter=None):
        """
        Checks the locations where the Emitter voltage and the extractor voltage are stable, 
        where the emitter voltage is larger than the extractor voltage and where both of them are above 100V.
        It calculates the Child-Langmuir for the data
        -------
        Populates the self.child_langmuir dictionary: It is a dictionary of lists. The keys are the heating currents.
        The first element of the list is a numpy array containing the voltage difference in kV, the second is the corresponding beam current in mA
        """
        
        if ("Emitter overcurrent" not in list(self.data.keys())) or ("Extractor overcurrent" not in list(self.data.keys())):
            self.get_overcurrent()
        
        em_voltage_der = np.abs(self.data['HV Em Meas Voltage'].data[1:]-self.data['HV Em Meas Voltage'].data[:-1])
        ex_voltage_der = np.abs(self.data['HV Ex Meas Voltage'].data[1:]-self.data['HV Ex Meas Voltage'].data[:-1])
        
        if neutralizer_shutter is None:
            rel_time = np.where((em_voltage_der<0.1)*\
                                (ex_voltage_der<0.1)*\
                                (self.data['HV Em Meas Voltage'].data[1:]>0.1)*\
                                (self.data['HV Ex Meas Voltage'].data[1:]<self.data['HV Em Meas Voltage'].data[1:]))
        else:
            rel_time = np.where((em_voltage_der<0.1)*\
                                (ex_voltage_der<0.1)*\
                                (abs(self.data['Neut Shut Closed'].data[1:] - neutralizer_shutter)<0.5)*\
                                (self.data['HV Em Meas Voltage'].data[1:]>0.1)*\
                                (self.data['HV Ex Meas Voltage'].data[1:]<self.data['HV Em Meas Voltage'].data[1:]))
        voltage_difference = self.data['HV Em Meas Voltage'].data[1:][rel_time] - self.data['HV Ex Meas Voltage'].data[1:][rel_time]
        heating_current = self.data['Emit Current A'].data[1:][rel_time]
        beam_current = self.data["Emitter overcurrent"].data[1:][rel_time]

        #rounds the heating_current to 1A accuracy
        heating_current_approx = (heating_current+0.5).astype(int)
        heating_current_values = np.unique(heating_current_approx)
        
        self.child_langmuir = dict()

        for index, heating_curr in enumerate(heating_current_values):
            rel_points = np.where(heating_current_approx == heating_curr)
            self.child_langmuir[str(heating_curr)] = [voltage_difference[rel_points],
                                                      beam_current[rel_points]]
        
    def plot_child_langmuir(self, newfigure=True, plotcolor=None, neutralizer_shutter=None,
                            label=None, alpha=None):
        if (hasattr(self, "child_langmuir") is False) or neutralizer_shutter is not None:
            self.get_child_langmuir(neutralizer_shutter=neutralizer_shutter)

        if newfigure is True:
            plt.figure()
            plt.title(f"Child-Langmuir plot based on BORI logs")
        
        color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
                      "tab:purple", "tab:brown", "tab:pink", "tab:gray",
                      "tab:olive", "tab:cyan"]
        marker_list = ["o", "v", "*", "D"]
        
        for index, heating_curr in enumerate(self.child_langmuir.keys()):
            if plotcolor is None:
                current_color = color_list[index%len(color_list)]
                current_marker = marker_list[(index//len(color_list))%len(marker_list)]
            else:
                current_color = plotcolor
                current_marker = marker_list[index%len(marker_list)]
            label=f"{heating_curr}A ({self.date})"
            if alpha is None:
                alpha = 0.2
            plt.scatter(self.child_langmuir[heating_curr][0],
                        self.child_langmuir[heating_curr][1],
                        c=current_color,
                        marker=current_marker,
                        label=label,
                        alpha=alpha)
        plt.xlabel("Extraction voltage [kV]")
        plt.ylabel("Emitter extra current [mA]")
        plt.legend()
        
    def plot_neutralizer(self, newfigure=True, plotcolor=None):

        if ("Emitter overcurrent" not in list(self.data.keys())) or ("Extractor overcurrent" not in list(self.data.keys())):
            self.get_overcurrent()
        
        em_voltage_der = np.abs(self.data['HV Em Meas Voltage'].data[1:]-self.data['HV Em Meas Voltage'].data[:-1])
        ex_voltage_der = np.abs(self.data['HV Ex Meas Voltage'].data[1:]-self.data['HV Ex Meas Voltage'].data[:-1])
        

        rel_time = np.where((em_voltage_der<0.1)*\
                            (ex_voltage_der<0.1)*\
                            (self.data['HV Em Meas Voltage'].data[1:]>0.1)*\
                            (self.data['HV Ex Meas Voltage'].data[1:]<self.data['HV Em Meas Voltage'].data[1:])*\
                            (self.data['Neut Shut Closed'].data[1:]==0)*\
                            (self.data['FC2 in'].data[1:]==1)*\
                            (self.data['FC1 in'].data[1:]==0) )
                
            
        oven_temp = self.data['TC Oven Bottom'].data[1:][rel_time]
        beam_current = self.data["Emitter overcurrent"].data[1:][rel_time]+self.data["Extractor overcurrent"].data[1:][rel_time]
        fc2_current = self.data["FC2 Resistor Current mA"].data[1:][rel_time]
        
        #rounds the heating_current to 1A accuracy
        oven_temp_approx = (oven_temp+0.5).astype(int)
        oven_temp_values = np.unique(oven_temp_approx)
        
        if newfigure is True:
            plt.figure()
            plt.title(f"Neutralizer efficiency based on BORI logs")
        
        if plotcolor is None:
            plotcolor = "tab:blue"
        
        plt.scatter(oven_temp_approx,
                    fc2_current/beam_current,
                    c=plotcolor,
                    alpha=0.2,
                    label=f"{self.date}")
        plt.xlabel("Oven bottom temperature [C]")
        plt.ylabel("FC2 current/Beam current")
        plt.legend()
    
    def get_vap_pressure(self, material="Na"):
        datatemp = copy.deepcopy(self.data['TC Oven Top'])
        datatemp.data_unit.name = "Pressure"
        datatemp.data_unit.unit = "mbar"
        datatemp.data = vap_pressure(self.data['TC Oven Top'].data, material=material)
        datatemp.name = "Vap Press Oven Top"
        self.data["Vap Press Oven Top"] = copy.deepcopy(datatemp)
        datatemp.data = vap_pressure(self.data['TC Oven Bottom'].data, material=material)
        datatemp.name = "Vap Press Oven Bottom"
        self.data["Vap Press Oven Bottom"] = copy.deepcopy(datatemp)
        datatemp.data = vap_pressure(self.data['TC Torus Side Cone'].data, material=material)
        datatemp.name = "Vap Press Torus Side Cone"
        self.data["Vap Press Torus Side Cone"] = copy.deepcopy(datatemp)
        datatemp.data = vap_pressure(self.data['TC Emit Side Cone'].data, material=material)
        datatemp.name = "Vap Press Emit Side Cone"
        self.data["Vap Press Emit Side Cone"] = copy.deepcopy(datatemp)
        

def find_files(startdate=None,starttime='0000',start_datetime=None,
               enddate=None,endtime='2359',end_datetime=None,
               datapath='',UTC_offset_minutes=None,verbose=True
               ):
    """
    Finds the files in which data for the indicated time interval (local time) can be found.

    Parameters
    ----------
    startdate : string, optional
        The start date YYYYMMDD. If not set the actual daye is assumed.
    starttime : string, optional
        The start time in local time, format: HHMM. The default is '0000'. The time is local time
    start_datetime: numpy datetime64 object
        Alternative to specify the start time. If this is set startdate and stattime is not used.
    enddate: string, optional
        The end date YYYYMMDD. If not set the same day is assumed as the start date.
    endtime: string, optional
        The end time in local time, format: HHMM. The default is '2459'. The time is local time.
    end_datetime: numpy datetime64 object
        Alternative to specify the end time. If this is set enddate and endtime is not used.
    datapath : string
        The data path.
    UTC_offset_minutes : int, optional
        The time offset of the local time realtive to UTC. 
        If None it will be determined from the date.
    verbose: boolean
        Print messages about progress of processing.

    Raises
    ------
    ValueError
        No suitable file found.

    Returns
    -------
    Returns the filenames, start and end times sorted by increasing starttime.
    file_list: list of strings
        The list of suitable filenames (with full path)
    file_starttime_list: list of np.datetime64 objects
        The start times of the files in local times.
    file_endtime_list: list of np.datetime64 objects
        The end times of the files in local times.
        

    """
    
    if (start_datetime is None):
        if (startdate is None):
            _startdate = str(datetime.date.today().year) + str(datetime.date.today().month) + str(datetime.date.today().day)
        else:
            _startdate = startdate
        _start_datetime = np.datetime64(_startdate[:4]+'-'+_startdate[4:6]+'-'+_startdate[6:8]+'T'+
                                        starttime[:2]+':'+starttime[2:4]
                                        )
    else:
        _start_datetime = start_datetime

    if (end_datetime is None):        
        if (enddate is None):
            _enddate = _startdate
        else:
            _enddate = enddate
        _end_datetime = np.datetime64(_enddate[:4]+'-'+_enddate[4:6]+'-'+_enddate[6:8]+'T'+
                                      endtime[:2]+':'+endtime[2:4]
                                      )
    else:
        _end_datetime = end_datetime
        
    date_str =str(_start_datetime)
    _date = date_str[:4] + date_str[5:7] + date_str[8:10]
    _UTC_offset_minutes = UTC_offset(UTC_offset_minutes=UTC_offset_minutes,date=_date)   
    file_starttime_list = []
    file_endtime_list = []
    fname_list = []
    # Starting search 13 h before the desired time as one log file is about 11 hour long
    act_datetime = _start_datetime - np.timedelta64(13,'h')
    dirname = datapath
    try:
        l = os.listdir(dirname)
    except FileNotFoundError:
        raise ValueError("No logfiles found, bad datapath?")
            
    for f in l:
        ind = f.find('.tdms')
        if (ind < 0):
            continue
        if (f[ind:] != '.tdms'):
            continue
        if (f[:5] == 'bori_'):
            try:
               file_start_datetime = np.datetime64(f[5:9]+'-'+f[9:11]+'-'+f[11:13]+'T'+f[14:16]+':'+f[16:18]+':'+f[18:20]) 
            except ValueError:
                continue
        else:
            try:
               file_start_datetime = np.datetime64(f[0:4]+'-'+f[4:6]+'-'+f[6:8]+'T'+f[9:11]+':'+f[11:13]+':'+f[13:15]) 
            except ValueError:
                continue            
        if (file_start_datetime > _end_datetime):
            continue
        if (file_start_datetime < _end_datetime - np.timedelta64(2,'D')):
            continue
        fname = os.path.join(dirname,f)
        if (verbose):
            print('Checking {:s}'.format(fname),flush=True)
        try:
            with TdmsFile.open(fname) as tdms_file:                    
#                tdms_version =  tdms_file.properties['Version']
                try:
                    ch_t = tdms_file['MonitorData']['TimeStamp']
                    file_start_datetime_from_file = ch_t[0] + np.timedelta64(_UTC_offset_minutes, 'm')
                    file_end_datetime_from_file = ch_t[-1] + np.timedelta64(_UTC_offset_minutes, 'm')
                    MAXDIFF = 500 # maximum time difference between time in file name and contents [s]
                    for i in range(11):
                        if (abs(file_start_datetime_from_file - file_start_datetime) > np.timedelta64(MAXDIFF,'s')):  
                            file_start_datetime_from_file = ch_t[i + 1] + np.timedelta64(_UTC_offset_minutes, 'm')
                            print("Warning: timestamp [{:d}] and time from file name are significantly different. (file: {:s})".format(i,fname),flush=True)
                        else:
                            break
                    if (abs(file_start_datetime_from_file - file_start_datetime) > np.timedelta64(MAXDIFF,'s')):  
                        print("Warning: Cannot find reasonable start time in file (file: {:s})".format(fname),flush=True)
                        print("Omitting this file.",flush=True)
                        continue                 
                    if (file_start_datetime_from_file > _end_datetime):
                        continue
                    for i in range(10):
                        if (file_end_datetime_from_file < file_start_datetime_from_file):
                            file_end_datetime_from_file = ch_t[-(2 + i)] + np.timedelta64(_UTC_offset_minutes, 'm')
                            print("Warning: timestamp [{:d}] in file is bad. (file: {:s})".format(-(i + 1),fname),flush=True)
                        else:
                            break
                except (IndexError, IOError, ValueError):
                    print("Broken file? {:s}".format(fname))
                    continue
        except ValueError:
            print("Broken file: {:s}".format(fname))
            continue
        if (file_end_datetime_from_file < file_start_datetime_from_file):  
            print("Warning: Cannot find reasonable end time in file (file: {:s})".format(fname),flush=True)
            print("Omitting this file.",flush=True)
            continue  
        if (file_end_datetime_from_file < _start_datetime):
            continue
        fname_list.append(fname)
        file_starttime_list.append(file_start_datetime)
        file_endtime_list.append(file_end_datetime_from_file)
    if (len(fname_list) == 0):
        raise ValueError('No data found.')
    rel_starttime = (np.array(file_starttime_list) - np.amin(np.array(file_starttime_list))).astype(float)
    ind = np.argsort(rel_starttime)
    fname_sorted = []
    starttime_sorted = []
    endtime_sorted = []
    for i in ind.tolist():
        fname_sorted.append(fname_list[i])
        starttime_sorted.append(file_starttime_list[i])
        endtime_sorted.append(file_endtime_list[i])
    return fname_sorted,starttime_sorted,endtime_sorted
            
def read_date_tdms(data_names=None,startdate=None,starttime='0000',start_datetime=None,
              enddate=None,endtime='2359',end_datetime=None,
               datapath='',UTC_offset_minutes=None,verbose=True):
    """
    Reads multiple data from the TDMS files specified by the start and end dates/times.

    Parameters
    ----------
    data_names : string or list of strings
        The channel names to read. The default is None.
    startdate : string, optional
        The start date YYYYMMDD. If not set the actual daye is assumed.
    starttime : string, optional
        The start time HHMM. The default is '0000'. The time is local time
    start_datetime: numpy datetime64 object
        Alternative to specify the start time. If this is set startdate and stattime is not used.
    enddate: string, optional
        The end date YYYYMMDD. If not set the same day is assumed as the start date.
    endtime: string, optional
        The end time HHMM. The default is '2459'. The time is local time.
    end_datetime: numpy datetime64 object
        Alternative to specify the end time. If this is set enddate and endtime is not used.
    datapath : string, optional
        The data path. The default is for the beam computer.
    UTC_offset_minutes : int, optional
        The time offset of the local time realtive to UTC. 
        If None it will be determined from the date.
    verbose: boolean
        Print messages about progress of processing.

    Raises
    ------
    ValueError
        No suitable file found.
    IOError
        Error deadeing data.

    Returns
    -------
    time : numpy datetime64 array
        The timstamps for the measurements. 
        Use (time-time[0])/np.timedelta64(1,'s') to convert to time relative to first time in seconds.
    data : List of numpy float arrays
        The data.
    data_unit : List of string(s)
        List of unit(s) corresponding to data

    """
    if (starttime is None):
        _starttime = '0000'
    else:
        _starttime= starttime
    if (endtime is None):
        _endtime = '2359'
    else:
        _endtime= endtime
    if (len(_starttime) != 4):
        raise ValueError('starttime should be given as HHMM')
    if (len(_endtime) != 4):
        raise ValueError('endtime should be given as HHMM')
    if (data_names is None):
        raise ValueError('Data names should be set.')
    if (type(data_names) is not list):
        _data_names = [data_names]
    else:
        _data_names = data_names
    if (start_datetime is None):
        if (startdate is None):
            _startdate = str(datetime.date.today().year) + str(datetime.date.today().month) + str(datetime.date.today().day)
        else:
            if (len(startdate) != 8):
                raise ValueError("date should be given as YYYMMDD")
            _startdate = startdate
        _start_datetime = np.datetime64(_startdate[:4]+'-'+_startdate[4:6]+'-'+_startdate[6:8]+'T'+
                                        _starttime[:2]+':'+_starttime[2:4]
                                        )
    else:
        _start_datetime = start_datetime

    if (end_datetime is None):        
        if (enddate is None):
            _enddate = _startdate
        else:
            if (len(enddate) != 8):
                raise ValueError("date should be given as YYYMMDD")
            _enddate = enddate
        _end_datetime = np.datetime64(_enddate[:4]+'-'+_enddate[4:6]+'-'+_enddate[6:8]+'T'+
                                      _endtime[:2]+':'+_endtime[2:4]
                                      )
    else:
        _end_datetime = end_datetime

    date_str =str(_start_datetime)
    _date = date_str[:4] + date_str[5:7] + date_str[8:10]
    _UTC_offset_minutes = UTC_offset(UTC_offset_minutes=UTC_offset_minutes,date=_date)   

    fnames, starttimes, endtimes = find_files(start_datetime=_start_datetime,end_datetime=_end_datetime,
                                              datapath=datapath,UTC_offset_minutes=_UTC_offset_minutes,verbose=verbose
                                              )

    time = np.ndarray(0,dtype=np.datetime64)
    data = [np.ndarray(0,dtype=float)]*len(_data_names)
    for fn in fnames:
        if (verbose):
            print('Processing {:s}'.format(fn),flush=True)
        with TdmsFile.open(fn) as tdms_file:
            # for channelname in list(sorted(tdms_file['MonitorData']._channels.keys())):
            #     print(channelname)
            try:
                tdms_version =  tdms_file.properties['Version']
            except KeyError:
                tdms_version = 1

            ch_t = tdms_file['MonitorData']['TimeStamp']
            try:
                t = ch_t[:] + np.timedelta64(_UTC_offset_minutes, 'm')
                ind = np.nonzero(np.logical_and(t >= _start_datetime,
                                                t <= _end_datetime 
                                                )
                                 )[0]
            except (IndexError, IOError, ValueError):
                print("Error reading from file: {:s}".format(fn))
                continue      
            if (len(ind) == 0):
                print("Cannot find data in time interval in file: {:s}".format(fn))
                continue
            data_unit=[]
            read_error = False
            data_save = copy.deepcopy(data)
            for i,signame in enumerate(_data_names):
                if (verbose):
                    print("Reading '{:s}'".format(signame),flush=True)
                ch = tdms_file['MonitorData'][signame] 
                d = ch[:]
                try:
                    data[i] = np.concatenate((data[i],d[ind].astype(float)))
                except IndexError:
                    read_error = True
                    break
            if (read_error):
                data = data_save
                if (verbose):
                    print("Broken file '{:s}'".format(fn),flush=True)           
            else:
                time = np.concatenate((time,t[ind]))
    if (len(time) < 2):
        raise ValueError("NO data found.")
    dt = np.diff((time - time[0]) / np.timedelta64(1,'s'))
    ind_bad = np.nonzero(dt <= 0)[0]
    if (len(ind_bad) != 0):
        print('Removed {:d} bad time points.'.format(len(ind_bad)),flush=True)
        ind_good = np.nonzero(dt > 0)[0]
        if (len(ind_good) == 0):
            raise ValueError("No good time points found.")
        time = np.concatenate((np.array([time[0]]),time[ind_good + 1]))
        for i in range(len(data)):
            data[i] = np.concatenate((np.array([data[i][0]]),data[i][ind_good + 1]))
    
    return time, data, data_unit

def read_channels(startdate,datapath,page='MonitorData'):
    channels = []
    file = find_files(startdate,datapath=datapath)
    with TdmsFile.open(file[0][0]) as tdms_file:
            ch_t = tdms_file[page]
    channel_list=ch_t.channels()
    for c in channel_list:
        channels.append(c.name)
    return channels

def channel_list(file,page='MonitorData'):
    """
    Return the list of channels in the TDMS file.

    Parameters
    ----------
    file : string
        The full file name.

    Returns
    -------
    channels : list of strings
        The channel list.

    """
    channels = []
    with TdmsFile.open(file) as tdms_file:
        ch_t = tdms_file[page]
    channel_list=ch_t.channels()
    for c in channel_list:
        channels.append(c.name)
    return channels
    
def page_list(file):
   with TdmsFile.open(file) as tdms_file:
       groups = tdms_file.groups()
       group_names = []
       for g in groups:
           group_names.append(g.name)
   return group_names
    
    
def read_exp_tdms(data_names, exp_id,datapath='/data/W7X/APDCAM'):
    dirname = os.path.join(datapath,exp_id)
    tdms_files = [filename for filename in os.listdir(dirname) if ("tdms" in filename and "tdms_index" not in filename)]
    all_data = dict()
    for file in tdms_files:
        with TdmsFile.open(os.path.join(dirname,file)) as tdms_file:
            for group in tdms_file.groups():
                currtime =  group['TimeStamp'][:]
                #getting the times
                if 'time_vect' not in locals():
                    time_vect = currtime
                else:
                    time_vect = np.concatenate([time_vect, currtime])
                # print([groupname.path.split("/")[2] for groupname in group.channels()])
                #getting the data
                for data in data_names:
                    currdata =  group[data][:]
                    #getting the times
                    if data not in all_data.keys():
                        all_data[data] = currdata
                    else:
                        all_data[data] = np.concatenate([all_data[data], currdata])
    data_unit = []
    #Sorting the time vector
    for key in all_data.keys():
        all_data[key] = np.array([x for _, x in sorted(zip(time_vect, all_data[key]))])
    time_vect = np.array(sorted(time_vect))
    
    return_data = [all_data[key] for key in data_names]
    
    return time_vect, return_data, data_unit

def plot_beamdata(startdate=None,starttime=None,endtime=None,enddate=None,datapath='',start_datetime=None,end_datetime=None,figure=None,
                  R_emit=81,R_ext=73, last_minutes=None):
    """
    Plot beam data.

    Parameters
    ----------
    startdate : string
        The start date of processing, YYYYMMDD
    starttime : string, optional
        The start time of processing, HHMM The default is None, whiche means 0000
    endtime : string, optional
        The end time of processing, HHMM. The default is None, which means 2400
    enddate : string, optional
        The end date of processing, YYYYMMDD. The default is None which means use the same date as startdate.
    datapath : string, optional
        The data path. The default is ''.
    start_datetime: numpy datetime64 object
        Alternative to specify the start time. If this is set startdate and stattime is not used.
   end_datetime: numpy datetime64 object
       Alternative to specify the end time. If this is set enddate and endtime is not used.
    figure : Matplotlib figure, optional
        The figure to use. The default is None, meaning make new figure.
    R_emit : float, optional
        The resistor on the emitter power supply [MOhm]. The default is 81.
        If None, will calculate from currents.
    R_ext : float, optional
        The resistor on the emitter power supply [MOhm]. The default is 73.
        If None, will calculate from currents.
    last_minutes : float
        Plot only the last minutes indicated by this argument.

    Returns
    -------
    None.

    """
    
    labelsize = 10
    linewidth = 1
    plt.rcParams['lines.linewidth'] = linewidth
    plt.rcParams['axes.linewidth'] = linewidth
    plt.rcParams['axes.labelsize'] = labelsize 
    plt.rcParams['axes.titlesize'] = labelsize 
    plt.rcParams['xtick.labelsize'] =  labelsize 
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = linewidth
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.minor.width'] = linewidth
    plt.rcParams['xtick.minor.size'] = 2
    plt.rcParams['ytick.labelsize'] = labelsize 
    plt.rcParams['ytick.major.width'] = linewidth
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.minor.width'] = linewidth
    plt.rcParams['ytick.minor.size'] = 2
    plt.rcParams['legend.fontsize'] = labelsize 
#    plt.rcParams['suptitle.fontsize'] = 10
    
    data_names = ['Emit Current A','HV Em Meas Voltage','HV Ex Meas Voltage','HV Em Meas Current','HV Ex Meas Current',
                  'TC Oven Top','TC Oven Bottom','TC Torus Side Cone','TC Emit Side Cone','FC1 in','FC2 in','FC Polarity',
                  'FC1 Resistor Current mA','FC2 Resistor Current mA','VG HighVac1','VG HighVac2',
                  'Neut Shut Closed']
    t,d,u = read_data(data_names=data_names,startdate=startdate,starttime=starttime,endtime=endtime,datapath=datapath)
    d_dict = {}
    for key,data in zip(data_names,d):
        d_dict[key] = data 
    
    if ((t[-1] -t[0]) / np.timedelta64(1,'s') > 300):
        time = (t - t[0]) / np.timedelta64(1,'m')
        time_unit = 'min'
    else:
        time = (t - t[0]) / np.timedelta64(1,'s')
        time_unit = 's'
        
    if (figure is None):
        plt.close('all')
        figure= plt.figure(figsize=(25,17))
    else:
        plt.figure(figure)
    gs = gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=plt.gca(),hspace=0.5,wspace=0.5)
    plt.suptitle('Start time: {:s}, R_emit={:3.0f}MOhm  R_ext={:3.0f}MOhm'.format(str(t[0]),R_emit,R_ext),fontsize=16)


    if (last_minutes is not None):
        if (time_unit == 'min'):
            xlim = [time[np.nonzero(time > time[-1] - last_minutes)[0][0]], time[-1] ]
        else:
            xlim = [time[np.nonzero(time > time[-1] - last_minutes)[0][0] * 60], time[-1]]
    else:
        xlim = None
    
    ax = plt.subplot(gs[0,0:2])
    plt.plot(time,d_dict['Emit Current A'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[A]')
    plt.title('Emitter heating current')
    if (xlim is not None):
        plt.xlim(*xlim)
    
    plt.subplot(gs[0,2:4],sharex=ax)
    plt.plot(time,d_dict['HV Em Meas Voltage'])
    plt.plot(time,d_dict['HV Ex Meas Voltage'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[kV]')
    plt.title('Accelerator voltages')
    plt.legend(['$HV_{emit}$','$HV_{ext}$'])
    
    plt.subplot(gs[1,2:4],sharex=ax)
    plt.plot(time,d[3])
    plt.plot(time,d[4])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[mA]')
    plt.title('HV PS currents')
    plt.legend(['$HV_{emit}$','$HV_{ext}$'])
    trange = plt.xlim()
    
    nonzero_current = np.logical_and(d[3] > 0.05,d[4] > 0.05)
    ind = np.nonzero(np.logical_and(d_dict['HV Em Meas Voltage'] < d_dict['HV Ex Meas Voltage'],nonzero_current))[0]
    if (len(ind) != 0):
        if (R_emit is None):
            R_emit = np.mean(d_dict['HV Em Meas Voltage'][ind] / d_dict['HV Em Meas Current'][ind])
        if (R_ext is None):
            R_ext = np.mean(d_dict['HV Ex Meas Voltage'][ind] / d_dict['HV Ex Meas Current'][ind])
        print("R_emit={:3.0f}MOhm  R_ext={:3.0f}MOhm".format(R_emit,R_ext))
        I_emit_ext = d_dict['HV Ex Meas Voltage'] / R_ext -  d[4] 
        I_beam = d_dict['HV Em Meas Current'] - d_dict['HV Ex Meas Voltage'] / R_emit - I_emit_ext
        plt.subplot(gs[2,2:4],sharex=ax)
        plt.plot(time,I_beam)
        plt.plot(time,I_emit_ext)
        plt.xlabel('Time [{:s}]'.format(time_unit))
        plt.ylabel('[mA]')
        plt.title('Beam current')
        plt.legend(['Beam current','E-ext current'])
        
    plt.subplot(gs[1,0:2],sharex=ax)
    legend = []
    ind = np.nonzero(np.logical_and(d_dict['TC Oven Bottom'] > 20, d_dict['TC Oven Bottom'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Oven Bottom'],20,300))
        legend.append('Oven')
    ind = np.nonzero(np.logical_and(d_dict['TC Oven Top'] > 20, d_dict['TC Oven Top'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Oven Top'],20,300))
        legend.append('Top')
    ind = np.nonzero(np.logical_and(d_dict['TC Torus Side Cone'] > 20, d_dict['TC Torus Side Cone'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Torus Side Conedata_unit'],20,300))
        legend.append('Cone torus side')
    ind = np.nonzero(np.logical_and(d_dict['TC Emit Side Cone'] > 20, d_dict['TC Emit Side Cone'] < 300))[0]
    if (len(ind) > 0):
        plt.plot(time,np.clip(d_dict['TC Emit Side Cone'],20,300))
        legend.append('Cone emit side')
    plt.legend(legend)    
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[C]')
    plt.title('Neutralizer temperatures')
    plt.xlim(*trange)

    ind_shut = np.nonzero(d_dict['Neut Shut Closed'] != 0)[0]
    fc2_act = d_dict['FC2 in']
    fc2_act[ind_shut] = 0                         
    ind1 = np.nonzero(d_dict['FC1 in'] == 1)[0]
    ind2 = np.nonzero(fc2_act == 1)[0]   
    if ((len(ind1) != 0) or (len(ind2 != 0))):
        legend = []
        plt.subplot(gs[3,2:4],sharex=ax)
        if (len(ind1) != 0):
            plt.plot(time[ind1],d_dict['FC1 Resistor Current mA'][ind1])
            legend.append('FC1 Resistor Current mA')
        if (len(ind2) != 0):
            plt.plot(time[ind2],d_dict['FC2 Resistor Current mA'][ind2])
            legend.append('FC2 Resistor Current mA')
        plt.legend(legend)    
        plt.xlabel('Time [{:s}]'.format(time_unit))
        plt.ylabel('[mA]')
        plt.title('Faraday cup currents')
        plt.xlim(*trange)
        

    plt.subplot(gs[2,0:2],sharex=ax)
    plt.plot(time,d_dict['VG HighVac1'])
    plt.plot(time,d_dict['VG HighVac2'])
    plt.legend(['High vac 1','High vac 2'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[mbar]')
    plt.title('Vacuum pressure')
    plt.xlim(*trange)
    plt.yscale('log')
    
    plt.subplot(gs[3,0:2],sharex=ax)
    I_extra_emit = d_dict['HV Em Meas Voltage'] / R_emit - d_dict['HV Em Meas Current']
    I_extra_ext = d_dict['HV Ex Meas Voltage'] / R_ext - d_dict['HV Ex Meas Current']
    plt.plot(time,I_extra_emit)
    plt.plot(time,-I_extra_ext)
    plt.legend(['Into Emit. PS','From Ext. PS'])
    plt.xlabel('Time [{:s}]'.format(time_unit))
    plt.ylabel('[mA]')
    plt.title('Extra currents')
    plt.xlim(*trange)

def vap_pressure(T, material="Na"):
    # Returns the equilibirum vapor pressure of the given alkali metal
    # T in celsius
    # Based on Honig (Courtesy RCA laboratories)
    T = T + 273.15
    
    temp_list_kelvin = np.array([300, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450,
                 460, 470, 480, 490, 500, 510, 520, 530])
    if material == "Na":
        if str(type(T)) != "<class 'numpy.ndarray'>":
            if T < 97.79:
                return np.nan
        press_list = [1e-8, 3e-8, 7e-8, 2e-7, 4e-7, 8e-7, 2e-6, 4e-6, 8e-6,
                      2e-5, 4e-5, 6e-5, 1e-4, 2e-4, 3.5e-4, 6.5e-4, 1.5e-3,
                      3.5e-3, 8e-3]
        vap_press = interp1d(temp_list_kelvin[1:], np.log(press_list), fill_value="extrapolate")
        result = np.exp(vap_press(T))
        if str(type(T)) == "<class 'numpy.ndarray'>":
            result[np.where(T < 97.79)] = np.nan
    if material == "K":
        if str(type(T)) != "<class 'numpy.ndarray'>":
            if T < 63.5:
                return np.nan
        press_list = [2e-8, 3e-6, 7e-6, 1.2e-5, 3e-5, 6e-5, 1e-4, 2e-4, 3.2e-4, 6e-4, 1e-3,
                      2e-3, 3e-3, 5e-3, 8e-3, 1.2e-2, 2.2e-2, 5e-2, 1e-1, 2e-1]
        vap_press = interp1d(temp_list_kelvin, np.log(press_list), fill_value="extrapolate")
        result = np.exp(vap_press(T))
        if str(type(T)) == "<class 'numpy.ndarray'>":
            result[np.where(T < 63.5)] = np.nan
    
#    plt.plot(temp_list_kelvin-273.15, press_list)
#    plt.yscale("log")
#    plt.ylim([1e-8, 1e-2])
    return result

if __name__ == "__main__":
    # plot_beamdata(startdate="20240924",datapath='/data',last_minutes=20)
    beam_log = BORIMonitor(exp_id="20240919.035")

            
