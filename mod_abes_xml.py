# -*- coding: utf-8 -*-
"""
Created on Mon May 27 22:10:49 2019

@author: Zoletnik
"""

import sys
import os

import flap

class ABES_Xml:
    def __init__(self,filename):
        self.filename = filename
        self.sections = []
        self.sectionElements = []
        self.top = None
        
    def createHead(self,shotID,startTime):
        self.top = ET.Element('ShotSettings',attrib={"Version":"1.1", "ShotID": shotID, "Experiment":"W7-X A-BES",
                                                     "Date":str(startTime.year)+"-"+str(startTime.month)+"-"+str(startTime.day),
                                                     "Time":str(startTime.hour)+":"+str(startTime.minute)+":"+str(startTime.second)})

    def addElement(self,section=None,element=None,value=None,unit=None,comment=None, value_type=None):
        if (section == None) or (element == None) or (value == None):
            return "ABES_Xml.addElement: Missing input data"
        

        if (type(value) == int):
            value_str = str(value)  
            type_str = 'long'
        elif (type(value) == float):
            value_str = str(value)
            type_str = "float"
        elif (type(value) == str):
            value_str = value
        else:
            return " ABES_Xml.addElement: unsuitable input data type"

        if (value_type != None):
            type_str = value_type
            
        if (unit == None):
            unit_str = "none"
        else:
            unit_str = unit
                
        try:
            section_index = self.sections.index(section)
            s = self.sectionElements[section_index]
        except:
            s = ET.SubElement(self.top, section)
            self.sections.append(section)
            self.sectionElements.append(s)
            section_index = len(self.sectionElements)-1
        
        if (comment == None):
            child = ET.SubElement(s, element, attrib={"Value":value_str, "Unit":unit_str,"Type":type_str,})
        else:    
            child = ET.SubElement(s, element, attrib={"Value":value_str, "Unit":unit_str,"Type":type_str,\
                                                  "Comment":comment})

    def writeFile(self):
        ET.ElementTree(self.top).write(self.filename)
 

def abes_get_config_2017(xml):
    retval = {}
    try:
        retval['version'] = xml.head.attrib['Version']
    except:
        retval['version'] = '1.0'
    try:
        retval['ShotID'] = xml.head.attrib['ShotID']
    except KeyError:
        raise ValueError("Invalid config file head format.")
    if (retval['version'] != '1.0'):    
        try:
            retval['Time'] = xml.head.attrib['Time']
            retval['Date'] = xml.head.attrib['Date']
        except KeyError:
            raise ValueError("Invalid config file head format.")
    try:
        retval['TriggerTime'] = Decimal(xml.get_element('System','TriggerTime')['Value'])
        retval['APDCAM_state'] = int(xml.get_element('APDCAM','State')['Value'])
        if (retval['APDCAM_state'] is 1):
            ADCDiv = Decimal(xml.get_element('APDCAM', 'ADCDiv')['Value'])
            ADCMult = Decimal(xml.get_element('APDCAM', 'ADCMult')['Value'])
            retval['APDCAM_f_ADC'] = Decimal(20e6) * ADCMult / ADCDiv
            samplediv = Decimal(xml.get_element('APDCAM', 'Samplediv')['Value'])
            retval['APDCAM_f_sample'] = retval['APDCAM_f_ADC'] / samplediv
            retval['APDCAM_sampletime'] = Decimal(1.) / retval['APDCAM_f_sample']
            retval['APDCAM_samplenumber'] = int(xml.get_element('APDCAM', 'SampleNumber')['Value'])
            retval['APDCAM_bits'] = int(xml.get_element('APDCAM', 'Bits')['Value'])
            trigger = Decimal(xml.get_element('APDCAM', 'Trigger')['Value'])
            if (trigger < 0):
                trigger = Decimal(0)
            retval['APDCAM_starttime'] = trigger + retval['TriggerTime']
            mask1 = int(xml.get_element('APDCAM', 'ChannelMask1')['Value'],16)
            mask2 = int(xml.get_element('APDCAM', 'ChannelMask2')['Value'],16)
            mask3 = int(xml.get_element('APDCAM', 'ChannelMask3')['Value'],16)
            mask4 = int(xml.get_element('APDCAM', 'ChannelMask4')['Value'],16)
            chmask = mask1 + (mask2 << 32) + (mask3 << 64) + (mask4 << 96)
            retval['APDCAM_chmask'] = chmask
            retval['APDCAM_bias1'] = float(xml.get_element('APDCAM','DetectorBias1')['Value'])
            retval['APDCAM_bias2'] = float(xml.get_element('APDCAM','DetectorBias2')['Value'])
            retval['APDCAM_det_temp'] = float(xml.get_element('APDCAM','DetectorTemp')['Value'])
            s = int(xml.get_element('APDCAM','ClockSource')['Value'])
            if (s == 1):
                retval['APDCAM_clock_source'] = 'External'
            else:
                retval['APDCAM_clock_source'] = 'Internal'
        chopmode = int(xml.get_element('Chopper','Mode')['Value'])
        clk = xml.get_element('Chopper','BaseClockFrequency')
        clk_freq = Decimal(clk['Value'])
        if (clk['Unit'] == 'Hz'):
            pass
        elif (clk['Unit'] == 'kHz'):
            clk_freq = clk_freq * Decimal(1000)
        elif (clk['Unit'] == 'MHz'):
            clk_freq = clk_freq * Decimal(1E6)
        else:
            raise ValueError("Unknown chopper clock frequency unit.")
        retval['Chopper clock'] = clk_freq
        sch = xml.get_element('Chopper','SchemeFileContents')['Value']
        sch = sch.replace("\n","")
        sch = sch.split('<NL>')
        retval['Chopper scheme'] = sch
        if (chopmode == 1):
            retval['Chopper mode'] = 'Camera'
            retval['CMOS exptime'] = Decimal(xml.get_element('CMOS','Exptime')['Value'])
            retval['CMOS frametime'] = Decimal(xml.get_element('CMOS','Frametime')['Value'])
            retval['CMOS frame number'] = int(xml.get_element('CMOS','FrameNumber')['Value'])
        else:
            retval['Chopper mode'] = 'Timed'
            retval['Chopper period'] = Decimal(xml.get_element('Chopper','PeriodTime')['Value']) \
                                         /Decimal(1000000)
    except ValueError as e:
        raise e
    ch_list = []
    fibre_list = []
    det_type_list = []
    adc_list = []
    for i in range(64):
        try:
            ch = xml.get_element('Optics', 'ADC' + str(i + 1))['Value']
            adc_list.append(i + 1)
            ch_arr = ch.split('-')
            if ((ch_arr[0][0] >= '1') and (ch_arr[0][0] <= '9')):
                ch_arr[0] = 'ABES-' + ch_arr[0]
            ch_list.append(ch_arr[0])
            fibre_list.append(ch_arr[1])
            if (len(ch_arr) > 2):
                det_type_list.append(ch_arr[2])
            else:
                det_type_list.append(' ')
        except Exception:
            pass

    def sortkey(in_string):
        signal_name = in_string[0]
        if (signal_name[0:5] == 'ABES-'):
            if (len(signal_name[5:]) == 1):
                signal_name = 'ABES-0' + signal_name[5:]
        return signal_name


    # Ordering the channel list
    ch_list, adc_list, fibre_list,det_type_list = zip(*sorted(zip(ch_list,
                                                                  adc_list,
                                                                  fibre_list,
                                                                  det_type_list),
                                                      key=sortkey))
    retval['ADC_list'] = list(adc_list)
    retval['signal_list'] = list(ch_list)
    retval['fibre_list'] = list(fibre_list)
    retval['det_type_list'] = list(det_type_list)

    return retval

class ShotTime():
    __init__(self,year=0,month=0,day=0,hour=0,minute=0,second=0)
    self.year = year
    self.month = month
    self.day = day
    self.hour = hour
    self.minute = minute
    self.second = second
    

def mod_abes_xml(exp_id):
    try:
        datapath_base = _options['Datapath']
    except (KeyError, TypeError):
        datapath_base = 'data'
    if (type(exp_id) is not str):
        raise ValueError("exp_id should be a string of format yyyymmdd.xxx")
    datapath = os.path.join(datapath_base,exp_id)
    xmlfile = os.path.join(datapath, exp_id + '_config.xml')
    old_xmlfile = os.path.join(datapath, exp_id + '_config_2017.xml')
    if (sys.platform == 'win32'):
        cmd = 'rename ' + xmlfile + ' ' + old_xmlfile 
    else:
        cmd = 'mv ' + xmlfile + ' '+ old_xmlfile
    d=subprocess.run([cmd],check=False,shell=True,stdout=subprocess.PIPE)
    old_xml = flap.FlapXml()
    try:
        old_xml.read_file(old_xmlfile)
    except Exception:
        raise IOError("Error reading XML file:" + old_xmlfile)
    m = ABES_Xml()
    
    shot_start_time=ShotTime(year=int(exp_id[0:4]),month=exp_id[4:6],day=exp_id[6:8],
                             hour=0, minute=0,second=0)
    m.createHead(exp_id,shot_start_time)
   
    m.addElement(section="System",element="TriggerTime", \
                 value=self.GUI_status.config.triggerTime,\
                 unit = 's',\
                 value_type='float',\
                 comment="Time of the trigger to APDCAM relative to heating start.")
    m.addElement(section="System",element="GateValve", \
                 value=with_valve,\
                 value_type='int',\
                 comment="1: Gate valve controlled, 0: Gate valve not controlled")

    m.addElement(section="Beam",element="Status", \
                 value=with_beam,\
                 value_type='int',\
                 comment="1: Beam operated, 0: beam not operated")
    if (with_beam == 1) :
        m.addElement(section="Beam",element="Energy", \
                 value=float(beam_settings.HV_emit_beam)/1000,\
                 value_type='float',\
                 unit = 'kV',\
                 comment="Beam energy")    
        m.addElement(section="Beam",element="ExtractorVoltage", \
                 value=float(beam_settings.HV_extrac_beam)/1000,\
                 value_type='float',\
                 unit = 'kV',\
                 comment="Beam extraction voltage")    

    if (self.GUI_status.APDCAM_connected):
        st = 1
    else:
        st = 0
    m.addElement(section="APDCAM",element="State", \
                 value= st,\
                 value_type='int',\
                 comment="APDCAM state; 0: off, 1: on")
    m.addElement(section="APDCAM",element="Type", \
                 value= "FC64",\
                 value_type='string',\
                 comment="APDCAM type: FC64 normally")        
    if (measurePara.externalTriggerPolarity == None):
        trig = -1
    else:
        trig = float(measurePara.triggerDelay)/1e6
    m.addElement(section="APDCAM",element="Trigger", \
                 value= float(trig),\
                 unit = 's',\
                 comment="Trigger: <0: manual,otherwise external with this delay")
    err,f,source = apd.getAdcClock()
    if (err != ""):
        return err
    m.addElement(section="APDCAM",element="ClockSource", \
                 value= int(source),\
                 value_type='int',\
                 comment="Clock source. 0: internal, 1: external")        
    m.addElement(section="APDCAM",element="ADCMult", \
                 value= int(status.CC_settings[codes_CC.CC_REGISTER_BASE_PLL_MULT]),\
                 value_type='int')
    m.addElement(section="APDCAM",element="ADCDiv", \
                 value= int(status.CC_settings[codes_CC.CC_REGISTER_BASE_PLL_DIV_ADC]),\
                 value_type='int')
    m.addElement(section="APDCAM",element="Samplediv", \
                 value= measurePara.sampleDiv,\
                 value_type='int')
    m.addElement(section="APDCAM",element="SampleNumber", \
                 value= measurePara.numberOfSamples,\
                 value_type='long')
    for i in range(4):
        m.addElement(section="APDCAM",element="ChannelMask"+str(i+1), \
                 value= "{:X}".format(measurePara.channelMasks[i]),\
                 value_type='long',\
                 comment="The channel mask for ADC block "+str(i+1)+"(hex). Each bit is one channel.")
    m.addElement(section="APDCAM",element="Bits", \
                 value= measurePara.bits,\
                 value_type='long',\
                 comment="The bit resolution of the ADC")
    m.addElement(section="APDCAM",element="DetectorBias1", \
                 value = float(status.HV_act[0]),\
                 value_type='float',\
                 unit = "V",\
                 comment="Detector 1 bias voltage. (APD)")
    m.addElement(section="APDCAM",element="DetectorBias2", \
                 value= float(status.HV_act[1]),\
                 value_type='float',\
                 unit = "V",\
                 comment="Detector 2 bias voltage. (MPPC)")
    m.addElement(section="APDCAM",element="DetectorTemp", \
                 value= float(status.temps[ABES_GUI_APDCAM_class.DETECTOR_TEMP_SENSOR-1]),\
                 value_type='float',\
                 unit = "C",\
                 comment="Detector 1 measured temperature.")

    if (self.GUI_status.CMOS_connected):
        st = 1
    else:
        st = 0
    m.addElement(section="CMOS",element="State", \
                 value= st,\
                 value_type='int',\
                 comment="CMOS state; 0: off, 1: on")
    m.addElement(section="CMOS",element="Exptime", \
                 value= float(CMOS_settings.exptime),\
                 value_type='float',\
                 unit = "ms",\
                 comment="CMOS camera exposure time")
    m.addElement(section="CMOS",element="Frametime", \
                 value= float(CMOS_settings.frametime),\
                 value_type='float',\
                 unit = "ms",\
                 comment="CMOS camera frame time")
    m.addElement(section="CMOS",element="FrameNumber", \
                 value= CMOS_settings.number_of_frames,\
                 value_type='int',\
                 comment="CMOS camera number of requested frames.")
    
    if (chop_settings != None):
        if (chop_settings.mode == chopper_settings_class.CAMERA):
            mode = 1
        elif (chop_settings.mode == chopper_settings_class.TIMED) :
            mode = 0
        else:
            return "Invalid chopper mode."
        m.addElement(section="Chopper",element="Mode", \
                 value=mode,\
                 value_type='int',\
                 comment="1: camera sync, 0: timed")
        m.addElement(section="Chopper",element="SchemeName", \
                 value=chop_settings.scheme_name,\
                 value_type='string',\
                 comment="Name of the timer scheme.")
        m.addElement(section="Chopper",element="SchemeFile", \
                 value=chop_settings.filename,\
                 value_type='string',\
                 comment="Name of the timer scheme file.")
        try:
            f = open(chop_settings.filename,"r")
        except:
            return "Error reading file: "+chop_settings.filename
        try:
            lines = f.readlines()
        except:
            f.close()
            return "Error reading file: "+chop_settings.filename
        f.close
        all_lines = ""
        for i in range(len(lines)):
            all_lines = all_lines+lines[i]+"<NL>"
        m.addElement(section="Chopper",element="SchemeFileContents", \
                 value=all_lines,\
                 value_type='string',\
                 comment="Contents of the timer scheme file.")
        if (chop_settings.mode == chopper_settings_class.TIMED):
            m.addElement(section="Chopper",element="PeriodTime", \
                     value=chop_settings.period_time,\
                     value_type='float',\
                     unit = 'microsec',\
                     comment="Period time of the chopping scheme for timed mode.")        
        m.addElement(section="Chopper",element="BaseClockFrequency", \
                 value=20e3/apd.CAMTIMER.clockDiv,\
                 value_type='float',\
                 unit = 'kHz',\
                 comment="The base clock frequency of the timer")
        if (beam_settings != None):
            m.addElement(section="Chopper",element="VoltBottom0", \
                     value=beam_settings.chop_bottom_0,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on bottom deflection plate when chopper control is 0")    
            m.addElement(section="Chopper",element="VoltBottom1", \
                     value=beam_settings.chop_bottom_1,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on bottom deflection plate when chopper control is 1")    
            m.addElement(section="Chopper",element="VoltTop0", \
                     value=beam_settings.chop_top_0,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on top deflection plate when chopper control is 0")                
            m.addElement(section="Chopper",element="VoltTop1", \
                     value=beam_settings.chop_top_1,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on top deflection plate when chopper control is 1")    
            m.addElement(section="Chopper",element="VoltLeft0", \
                     value=beam_settings.chop_left_0,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on left deflection plate when chopper control is 0")    
            m.addElement(section="Chopper",element="VoltLeft1", \
                     value=beam_settings.chop_left_1,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on left deflection plate when chopper control is 1")    
            m.addElement(section="Chopper",element="VoltRight0", \
                     value=beam_settings.chop_right_0,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on right deflection plate when chopper control is 0")    
            m.addElement(section="Chopper",element="VoltRight1", \
                     value=beam_settings.chop_right_1,\
                     value_type='float',\
                     unit = 'V',\
                     comment="Voltage on right deflection plate when chopper control is 1")    


        m.addElement(section="System",element="APD_H-Micrometer", \
                 value=self.GUI_status.config.micrometer_H,\
                 value_type='float',\
                 unit = 'mm',\
                 comment="Position of the H (beam-perp.) micrometer of the APD fibre plate.")   
        m.addElement(section="System",element="APD_V-Micrometer", \
                 value=self.GUI_status.config.micrometer_V,\
                 value_type='float',\
                 unit = 'mm',\
                 comment="Position of the V (beam-parallel) micrometer of the APD fibre plate.")   
 
   
    
    for i in range(64):
        try:
            ind = self.GUI_status.config.ADCs.index(i+1)
            m.addElement(section="Optics",element="ADC"+str(self.GUI_status.config.ADCs[ind]), \
                         value= self.GUI_status.config.signals[ind],\
                         value_type='string',\
                         comment="Optical channel - fibre number")
        except:
            pass
    
    try:
        m.writeFile()
    except: return "Error writing file "+fn
    
    return ""
"""

mod_abes_xml('20171207.027')