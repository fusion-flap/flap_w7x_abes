# -*- coding: utf-8 -*-
"""
Created on Mon May 27 22:10:49 2019

@author: Zoletnik
"""
def mod_abes_xml():
    

#def writeConfigFile(self,with_beam,with_valve,shot_start_time,datapath="data",shotID=None,CMOS_settings=None,chop_settings=None,\
#                    beam_settings=None):
#    """ Writes an XML configuration file.
#    Call this routine only after a succesful measure() !!!
#    """
    
    apd = self.GUI_status.APDCAM_reg
    apd.readStatus()
    measurePara = apd.measurePara
    status = apd.status
    codes_CC = apd.codes_CC
    
    fn = datapath+"/"+shotID+"_config.xml"
    m = ABES_Xml(fn)
    m.createHead(shotID,shot_start_time)
   
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
