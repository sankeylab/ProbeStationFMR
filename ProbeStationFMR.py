# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:29:09 2015

@author: CHILDRESSLAB
"""

###################
# To do
###################

# Add functions / buttons for moving the magnet arbitrarily
# Add support for 2d frequency-power sweeps
# Add averaging functions for multiple iterations
# Add support for power sweeps
# Add support for arbitrary power / frequency lists
# Add diagnostic plots


# Import the basics
import numpy as _n
import scipy as _s
import scipy.optimize as spop
import time

# DAQ stuff
import spinmob.daqmx as _daqmx
import pyvisa
import visa as _visa

# Plotting/gui stuff
import matplotlib.pyplot as plt
import pyqtgraph as pg
#import pyqtgraph.widgets.MatplotlibWidget as mpw
import spinmob.egg as egg
from gui.ColorWFPlot import *

# File import stuff
import sys
import csv
sys.path.append("M:\\Lab Files\\Magnetometry\\Simulations, Theory\\Spin Hall Effect Macrospin Simulations")
sys.path.append("X:\\Lab Files\\Magnetometry\\Simulations, Theory\\Spin Hall Effect Macrospin Simulations")
import os
from Tkinter import *
import tkFileDialog as filedialog

# Import analysis routines
import SHEmacrospin_fitters as FMRfit

# Import instruments
from instruments import signalgenerators
from instruments import camgui

# Import camera libraries
try:
    import cv2 as _cv2
    cam = True
except ImportError:
    cam = False
    print "No camera interface"
    
# Python ActiveX Client for interfacing with Labview    
import win32com.client  




###############################################
### HARDWARE INITIALIZATION
###############################################   
daqmx = _daqmx.daqmx_system()


# Generate list of connected RF sources     
resource_manager = _visa.ResourceManager()
device_names     = list(resource_manager.list_resources()); 
sg_strings   = ["121-227450000-0385" , "176142"       , "178228"      , "113055"]
sg_names     = ["Anapico"            , "SMB100A_CMC"  , "SMB100A_RnS" , "SMB100A_6GHz"]
sg_connected = []
i = 0
for j in sg_strings:
    for k in range(len(device_names)):         
        if device_names[k].find(j) >= 0:
            sg_connected.append(sg_names[i])
    i += 1
if len(sg_connected) == 0:
    print "WARNING: No signal generator connected"
sg_connected.append("None")


def init_sources():
    """
    Initializes the main and calibration signal generators.
    """
    global a, b
    
    source = cfg_hardware["RFSource"]
    if   source == "Anapico":    
        a = signalgenerators.Anapico()
    elif source == "SMB100A_CMC":
        a = signalgenerators.SMB100A_CMC()
    elif source == "SMB100A_RnS":
        a = signalgenerators.SMB100A_RnS()
    elif source == "SMB100A_6GHz":
        a = signalgenerators.SMB100A_6GHz()
    elif source == "None":		
        a = signalgenerators.NoneGen()
        
    if source != "None":
        a.write("PULM:SOUR EXT")
        a.write("PULM:STAT OFF")
        a.write("OUTP OFF")

    # calibration source    
    source = cfg_hardware["CalibSource"]
    if source != cfg_hardware["RFSource"]:
        if   source == "Anapico":    
            b = Anapico()
        elif source == "SMB100A_CMC":
            b = SMB100A_CMC()
        elif source == "SMB100A_RnS":
            b = SMB100A_RnS()
        elif source == "SMB100A_6GHz":
            b = SMB100A_6GHz()
    
        if source != "None":
            b.write("PULM:SOUR EXT")
            b.write("PULM:STAT OFF")
            b.write("OUTP OFF")
    else:
        cfg_hardware["CalibSource"] = "None"


def init_crystaldet():
    """
    Initializes the calibration file for the crystal detector.
    """
    global cd_calib
    
    # Default calibration (voltage proportional to power)
    cd_calib = _n.vstack( ( _n.array(range(-30,31)) , 10**(_n.array(range(-30,31))/10) ) )   
    
    # Import calibration file
    cd = cfg_hardware["CrystalDet"]
    
    try:
        f = open("instruments/" + cd + ".cfg")
        lines = f.readlines()
        cd_calib = _n.zeros([len(lines), 2])
        k = 0
        for line in lines:
            cd_calib[k, 0] = float(line.split()[0])
            cd_calib[k, 1] = float(line.split()[1])
            k += 1
        # Put values in increasing order
        if cd_calib[-1, 1] < cd_calib[1,1]: cd_calib = _n.flipud(cd_calib)
        print "Crystal detector: calibration loaded"
    except IOError:
        if cd != "None": print "Crystal detector: no calibration file found"

        

# Conex actuators interface
from instruments.Conex import CONEX_controller as CONEX_c

conex_status = ["None"]
try:
    CONEX_c.CONEX("COM5")
    CONEX_c.CONEX("COM3")
    CONEX_c.CONEX("COM4")
    conex_status.append("Available")
    from Conex import ConexXYZMove as conex_move
    
except pyvisa.errors.VisaIOError:
    print "No magnet actuator interface available"        


###############################################
##### GUI LAYOUT
###############################################      

# Window properties
window = egg.gui.Window("Probe Station FMR Measurements")
window.set_size([900,680])
window.set_position([0,0])
window.set_column_stretch(1,100)



### Upper left global controls
g_controls = window.place_object(egg.gui.GridLayout(False))

# First Row: Buttons and settings for the RF-Sweep
b_sweep      = g_controls.place_object(egg.gui.Button("Sweep")).set_checkable(True).set_width(60)
b_fourpt     = g_controls.place_object(egg.gui.Button("IV Curve")).set_checkable(True).set_width(60)
b_isweep     = g_controls.place_object(egg.gui.Button("I Sweep")).set_checkable(True).set_width(60)
b_bsweep     = g_controls.place_object(egg.gui.Button("B Sweep")).set_checkable(True).set_width(60)
g_controls.new_autorow()
g_controls.place_object(egg.gui.Label('Iteration:'), alignment=2)
n_fire_count = g_controls.place_object(egg.gui.NumberBox(0, 1, [0,None],True)).set_width(70)


############################
### Tabs - Experiment Config
############################
window.new_autorow()
tabs_controls = window.place_object(egg.gui.TabArea(False))


t_sweep       = tabs_controls.add_tab("Sweep")
cfg_sweep     = t_sweep.place_object(egg.gui.TreeDictionary("cfg_sweep.cfg")).set_width(320)
t_fourpt      = tabs_controls.add_tab("IV Curve")
cfg_fourpt    = t_fourpt.place_object(egg.gui.TreeDictionary("cfg_fourpt.cfg")).set_width(320)
t_pcalib      = tabs_controls.add_tab("Power Calib.")
cfg_pcalib    = t_pcalib.place_object(egg.gui.TreeDictionary("cfg_pcalib.cfg")).set_width(320)
t_isweep      = tabs_controls.add_tab("I Sweep")
cfg_isweep    = t_isweep.place_object(egg.gui.TreeDictionary("cfg_isweep.cfg")).set_width(320)
t_bsweep      = tabs_controls.add_tab("B Sweep")
cfg_bsweep    = t_bsweep.place_object(egg.gui.TreeDictionary("cfg_bsweep.cfg")).set_width(320)


# Hardware settings
t_hardware   = tabs_controls.add_tab("Hardware")
cfg_hardware = t_hardware.place_object(egg.gui.TreeDictionary("cfg_hardware.cfg")).set_width(320)
daqmx_device_names = daqmx.get_device_names() # get the names of all devices attached to the system

# calibration source
cfg_hardware.add_parameter("RFSource",    "Anapico",    type='list',    values=sg_connected)
cfg_hardware.add_parameter("CalibSource", "Anapico",    type='list',    values=sg_connected)
cfg_hardware.add_parameter("SaturationPow", 25,         type='float', limits = (0,None), siPrefix = True, suffix = 'dBm') #Adrian addition
cfg_hardware.add_parameter('SettlingTime', 0.001, type = 'float', limits = (0, None), siPrefix = True, suffix = 's', minStep = 0.001)
cfg_hardware.add_parameter("DCSourceR", 1000, type = 'float', limits = (0,None), siPrefix = True, suffix = 'Ohm')
crystal_detectors = ["None", "Narda4503A"]
cfg_hardware.add_parameter("CrystalDet" , crystal_detectors[0], type='list', values=crystal_detectors)
cfg_hardware.add_parameter("SA/LowPass",     100e6,     type='float', limits = (0,None), siPrefix = True, suffix = 'Hz')

# analog in
cfg_hardware.add_parameter("AI/Device",             0,      type='list',  values=daqmx_device_names)
cfg_hardware.add_parameter("AI/Sweep",            '0',      type='str')
cfg_hardware.add_parameter("AI/FourPT",           '1',      type='str')
cfg_hardware.add_parameter("AI/dVdI",             '2',      type='str')
cfg_hardware.add_parameter('AI/Crystal_Detector', '3',      type='str')
cfg_hardware.add_parameter("AI/Range",    '10.0',  type='str')
cfg_hardware.add_parameter("AI/Rate",      1e3,    type='float', limits=(10,None), siPrefix=True, suffix='Hz', dec=True, minStep=1)

# analog out
cfg_hardware.add_parameter("AO/Device",         0,      type='list',    values=daqmx_device_names)
cfg_hardware.add_parameter("AO/PULM",           0,      type='int',     limits=(0,None))
cfg_hardware.add_parameter("AO/DC",             1,      type='int',     limits=(0,None))
cfg_hardware.add_parameter("AO/DCRampTime",    10,    type='float',   limits=(0,None), siPrefix=True, suffix='s/mA')
cfg_hardware.add_parameter("AO/Trigger",        2,      type='int',     limits=(0,None))
cfg_hardware.add_parameter("AO/Range",          10.0,   type='float',   limits=(0,10.0), siPrefix=True, suffix='V',  dec=True, minStep=0.001)
cfg_hardware.add_parameter("AO/Rate",           1e3,    type='float',   limits=(10,None), siPrefix=True, suffix='Hz', dec=True, minStep=1)

# amplifier
cfg_hardware.add_parameter("Amp/Gain", 2500, type='int', limits=(0,10000))
cfg_hardware.add_parameter("Amp/Attenuation", 0, type='float', limits=(None,0) , siPrefix=False, suffix='dB', dec = True, minStep = 0.01)

# DC source
cfg_hardware.add_parameter("DC/V2A",           0.02, type='float', limits=(None,None) , siPrefix=True, suffix='A/V', dec=True, minStep=0.01)
cfg_hardware.add_parameter("DC/MaxV",           5.0, type='float', limits=(0,10), siPrefix=True, suffix = 'V', dec=True, minstep=0.1)

# Magnet actuators
cfg_hardware.add_parameter("Conex/Status", conex_status[-1], type = 'list', values = conex_status)
cfg_hardware["Conex/Status"] = conex_status[-1]

# Labview communication
cfg_hardware.add_parameter("Labview/VI_Path", "", type="str")
b_labviewvi = cfg_hardware.add_button("Select Labview VI Path")


# Filters

# Sweep parameters
experiments = ["Driven Dynamics", "Spectrum Analyzer"]
cfg_sweep.add_parameter("Experiment", "Driven Dynamics", type='list', values = experiments)

pow_modes = ['Constant', 'Calibrated', 'Mixed']; sweep_modes = ['LIN', 'LOG']
cfg_sweep.add_parameter("Power",             0, type='float', limits=(-30,25)         , siPrefix=True, suffix='dBm', dec=False)
cfg_sweep.add_parameter("PowerMode",'Constant', type='list', values = pow_modes)
cfg_sweep.add_parameter("DC_Current",      0.0, type='float', limits=(-cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"],cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"]) , siPrefix=True, suffix='A' , step=0.0001)
cfg_sweep.add_parameter("Start",         1.0e9, type='float', limits=(100.0e6,20.0e9) , siPrefix=True, suffix='Hz' , dec=True, Step=1e8)
cfg_sweep.add_parameter("Stop",         10.0e9, type='float', limits=(100.0e6,20.0e9) , siPrefix=True, suffix='Hz' , dec=True, Step=1e8)
cfg_sweep.add_parameter("Steps",           100, type='int'  , limits=(1,None)         , siPrefix=True)
cfg_sweep.add_parameter("Mode",          'LIN', type='list' , values = sweep_modes)
cfg_sweep.add_parameter("Iterations",        1, type='int'  , limits=(0,None))

# Pulse modulation parameters
cfg_sweep.add_parameter("PULM/Period",         0.001, type='float', limits=(1e-6,1e-2)      , siPrefix=True, suffix='s' , dec=True)
cfg_sweep.add_parameter("PULM/Repeats",            1, type='int'  , limits=(1,1e6)) #Number of times that the pulse train is repeated for a single frequency
cfg_sweep.add_parameter("PULM/Amplitude",        3.0, type='float', limits=(0,5), suffix=' V' )
cfg_sweep.add_parameter("PULM/DC",                 0, type='float', limits=(-2,2), suffix=' V' )   

# DC current and 4-point measurement parameters
cfg_fourpt.add_parameter("Start",        -1e-3, type='float', limits=(-cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"],cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"]) , siPrefix=True, suffix='A' , step=0.0001)
cfg_fourpt.add_parameter("Stop",          1e-3, type='float', limits=(-cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"],cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"]) , siPrefix=True, suffix='A' , step=0.0001)
cfg_fourpt.add_parameter("Steps",          100, type='int'  , limits=(1,None)         , siPrefix=True)
cfg_fourpt.add_parameter("Mode",         'LIN', type='list' , values = sweep_modes)
cfg_fourpt.add_parameter("ModPeriod",     1e-3, type='float', limits=(1e-6,1e0)      , siPrefix=True, suffix='s' , dec=True, minStep=1)
cfg_fourpt.add_parameter("ModRepeats",     100, type='int'  , limits=(1,1e6)) # Number of modulation periods for a single current step
cfg_fourpt.add_parameter("ModAmp",     10.0e-6, type='float', limits=(0,cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"]/10)   , siPrefix=True, suffix='A' , step=1e-6)
cfg_fourpt.add_parameter("PeriodCutoff",   0.0, type='int'  , limits=(0,None)         , siPrefix=True, dec=False)

# Spectrum analyzer calibration parameters
dd_calib_modes = ['High Damping', 'Background Subtraction']
cfg_pcalib.add_parameter("DD_Calib/Mode", 'High Damping', type='list' , values = dd_calib_modes)
cfg_pcalib.add_parameter("DD_Calib/Iterations",        1, type='int'  , limits=(0,None)         , siPrefix=True)
cfg_pcalib.add_parameter("DD_Calib/DC_Current"  ,     1e-3, type='float', limits=(-cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"],cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"]) , siPrefix=True, suffix='A' , step=0.0001)
cfg_pcalib.add_parameter("DD_Calib/dVdI_Current",     1e-3, type='float', limits=(-cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"],cfg_hardware["DC/MaxV"]*cfg_hardware["DC/V2A"]) , siPrefix=True, suffix='A' , step=0.0001)
cfg_pcalib.add_parameter("DD_Calib/dVdI_Steps",         11, type='int'  , limits=(1,None)         , siPrefix=True)

# Spectrum analyzer calibration parameters
cfg_pcalib.add_parameter("SA_Calib/Range",      100e6, type='float', limits=(0,None)         , siPrefix=True, suffix='Hz', dec=False)
cfg_pcalib.add_parameter("SA_Calib/Steps",        100, type='int'  , limits=(0,None)         , siPrefix=True)
cfg_pcalib.add_parameter("SA_Calib/PowerStart",  -100, type='float', limits=(None,None)      , siPrefix=True, suffix='dBm')
cfg_pcalib.add_parameter("SA_Calib/PowerStop",   -100, type='float', limits=(None,None)      , siPrefix=True, suffix='dBm')
cfg_pcalib.add_parameter("SA_Calib/PowerSteps",     1, type='int'  , limits=(0,None)         , siPrefix=True, dec=False)

# I sweep parameters
cfg_isweep.add_parameter("Start", -0.008, type='float', limits=(-.01,.01) , siPrefix=True, suffix='A' , step=0.0001)
cfg_isweep.add_parameter("Stop", 0.008, type='float', limits=(-.01,.01) , siPrefix=True, suffix='A' , step=0.0001)
cfg_isweep.add_parameter("Steps", 17, type='int')

# B-Sweep parameters
cfg_bsweep.add_parameter("TrajectoryFile", None, type="str")
b_bsweepfile = cfg_bsweep.add_button("Select Trajectory File")
cfg_bsweep.add_parameter("Start",          1e-3, type='float', siPrefix=True, suffix='T' , dec=True)
cfg_bsweep.add_parameter("Stop",         100e-3, type='float', siPrefix=True, suffix='T' , dec=True)
cfg_bsweep.add_parameter("Steps",            10, type='int'  , limits=(1,None)         , siPrefix=True)
cfg_bsweep.add_parameter("Mode",          'LIN', type='list' , values = sweep_modes)
cfg_bsweep.add_parameter("Offsets/Alpha",     0, type='float', siPrefix=True, suffix='rad')
cfg_bsweep.add_parameter("Offsets/x0",        0, type='float', siPrefix=True, suffix='mm')
cfg_bsweep.add_parameter("Offsets/y0",        0, type='float', siPrefix=True, suffix='mm')
cfg_bsweep.add_parameter("Offsets/z0",        0, type='float', siPrefix=True, suffix='mm')




### Initialize hardware

# Initialize signal generator(s)
cfg_hardware.connect_signal_changed("RFSource", init_sources)
cfg_hardware.connect_signal_changed("CalibSource", init_sources)
init_sources()

# Initialize crystal detector   
cfg_hardware.connect_signal_changed("CrystalDet", init_crystaldet)
init_crystaldet()



### Signal generator communication
window.new_autorow()
aq_controls = window.place_object(egg.gui.GridLayout(False))
b_outp = aq_controls.place_object(egg.gui.Button("RF Off").set_checkable(True).set_checked(False)).set_width(70) 
b_pulm = aq_controls.place_object(egg.gui.Button("PULM On").set_checkable(True).set_checked(True)).set_width(70) ; pulm_state = 1   
b_sweepref   = aq_controls.place_object(egg.gui.Button("Reference")).set_checkable(True).set_width(70) 
b_resis   = aq_controls.place_object(egg.gui.Button("Measure R")).set_checkable(True).set_checked(True).set_width(70) 
#b_multiplot= aq_controls.place_object(egg.gui.Button("Plot Multiple"))

# Functions for switches
def outp_switch():
    """
    RF On / RF Off switch
    """
# FIRMWARE
    b = int(a.query("OUTP?"))
    if b == 1:
        a.write("OUTP 0")
        b_outp.set_checked(False)
        b_outp.set_text("RF Off")
    else:
        a.write("OUTP 1")
        b_outp.set_checked(True)
        b_outp.set_text("RF On")      
# connect this function 
window.connect(b_outp.signal_clicked,outp_switch) 

def outp_off():
    """
    Turn output off
    """    
    a.write("OUTP 0")
    b_outp.set_checked(False)
    b_outp.set_text("RF Off") 
# connect this function to a window close event
# window.connect(window._window.closeEvent,outp_off)    not sure how to implement...
    

def outp_state():
    """
    Queries the Anapico for its output state and sets the RF On/Off button accordingly
    """
    if cfg_hardware["RFSource"] != "None": 
        print 
        b = int(a.query("OUTP?"))
        if b == 1:
            b_outp.set_checked(True)
            b_outp.set_text("RF On")          
        else:    
            b_outp.set_checked(False)
            b_outp.set_text("RF Off")
outp_state()


def pulm_switch():
    """
    PULM On / PULM Off switch
    """ 
    global pulm_state    
    pulm_state = 1 - pulm_state
    if pulm_state == 1:
        b_pulm.set_text("PULM On")
    else:
        b_pulm.set_text("PULM Off")      
# connect this function 
window.connect(b_pulm.signal_clicked,pulm_switch)      


# Send arbitrary commands to the signal generator
window.new_autorow()
mw_gen_control = window.place_object(egg.gui.GridLayout(False),0,3)

mw_gen_combo = mw_gen_control.place_object(egg.gui.ComboBox())
mw_commands = ["Command to signal generator",
               "OUTP ","FREQ:FIX ","FREQ:MODE ","POW ","LIST:COUN ","LIST:DEL ","LIST:FREQ ",
               "LIST:MAN ","LIST:MODE ","LIST:POWER ","PULM:STAT "]
for n in mw_commands: mw_gen_combo.add_item(n);
mw_gen_combo.insert_separator(1)

mw_gen_command = mw_gen_control.place_object(egg.gui.TextBox('',False))

# Corresponding functions
def mw_gen_combotext():
    if mw_gen_combo.get_current_index() != 0:
        string = mw_gen_combo.get_text()
        mw_gen_command.set_text(string)
        mw_gen_combo.set_current_index(0)
# connect this function      
window.connect(mw_gen_combo.signal_activated,mw_gen_combotext)

def mw_gen_comm():
    """
    Sends arbitrary commands to the signal generator
    """
    mw_gen_str = mw_gen_command.get_text()
    
    if "?" not in mw_gen_str:
        a.write(mw_gen_str)
# FIRMWARE
    else:
        a.query(mw_gen_str)
    outp_state() # Update RF ON/OFF switch
# connect this function    
window.connect(mw_gen_command.signal_return_pressed,mw_gen_comm)




#########################
### Tabs - Plots
#########################
tabs_plots = window.place_object(egg.gui.TabArea(False), 1,0, row_span=4, alignment=0)

# tabs
t_cam       = tabs_plots.add_tab("Camera")
#t_ao        = tabs_plots.add_tab("Output")
t_calib    = tabs_plots.add_tab("Power Calibration")
t_raw       = tabs_plots.add_tab("Raw (Sweep)")
t_sweep    = tabs_plots.add_tab("Sweep")
t_fmr       = tabs_plots.add_tab("FMR")
t_rawiv     = tabs_plots.add_tab("Raw (IV)")
t_iv        = tabs_plots.add_tab("IV Curve")
t_isweep    = tabs_plots.add_tab("I Sweep")
t_bsweep    = tabs_plots.add_tab("B Sweep")
t_test      = tabs_plots.add_tab("Test")

# Camera interface
t_cam_row1      = t_cam.place_object(egg.gui.GridLayout())
t_cam_row1.set_column_stretch(2)
video_input     = t_cam_row1.place_object(egg.gui.NumberBox(0, 1).set_width(40))
video_label     = t_cam_row1.place_object(egg.gui.Label('Input channel'))
button_stream   = t_cam_row1.place_object(egg.gui.Button('Stream').set_checkable(True))
# second row: image
t_cam_row2 = t_cam.place_object(egg.gui.GridLayout(), 0,1, alignment=0)

if cam: 
    image_raw = t_cam_row2.place_object(camgui.ImageWithButtons(window), alignment=0)


# databox plotters and goodies for the power calibration tab
tabs_calib   = t_calib  .place_object(egg.gui.TabArea(False), 1,0, row_span=4, alignment=0)
t_calib_dd   = tabs_calib.add_tab("Driven Dynamics")
plot_calib   = t_calib_dd.place_object(egg.gui.DataboxPlot("*.calib" , "plotsw_calib.cfg"), column_span=7)
plot_calib.button_autoscript.set_checked(False)
t_calib_dd.new_autorow()
b_calib       = t_calib_dd  .place_object(egg.gui.Button("Power Calibration"), column = 0, alignment = 0).set_checkable(True)
b_calibthird  = t_calib_dd  .place_object(egg.gui.Button("Step 3"), column = 5, alignment = 0).set_checkable(True).set_checked(True)
b_calibreset  = t_calib_dd  .place_object(egg.gui.Button("Reset"), column = 6, alignment = 0)

t_calib_sa    = tabs_calib.add_tab("Spectrum Analyzer")
plot_calib_sa = t_calib_sa.place_object(ColorWFPlot(size = [cfg_sweep['Steps'], cfg_pcalib['SA_Calib/PowerSteps']], \
    bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_pcalib['SA_Calib/PowerStart'], cfg_pcalib['SA_Calib/PowerStop']]))
plot_calib_sa.add_databox('P_out')
plot_calib_sa.add_databox('V_base')
plot_calib_sa.add_databox('Gain')
plot_calib_sa.add_databox('Gain_dB')
plot_calib_sa.set_databox(3)
plot_calib_sa.set_labels(['Frequency (Hz)', r'$P_{in}$ (dBm)', r''])



# databox plotters for the raw / four point measurements
plot_raw      = t_raw     .place_object(egg.gui.DataboxPlot("*.raw"    , "plotsw_raw.cfg"), column_span=2)
raw_script = "x = [ d[0] ]\ny = [ d[1] ]\n\nxlabels=\'t\'\nylabels=\'V\'"; plot_raw.script.set_text(raw_script)
t_raw.new_autorow()



# databox plotters and goodies for the lock-in tab
plot_sweep   = t_sweep  .place_object(egg.gui.DataboxPlot("*.swp" , "plotsw_sweep.cfg"), column_span=2)
plot_sweep.button_autoscript.set_checked(False)
sweep_script = "x = [ d[0], d[0] ]\ny = [ d[5], d[6] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=\'V\'"; plot_sweep.script.set_text(sweep_script)
t_sweep.new_autorow()
b_send_sweep = t_sweep  .place_object(egg.gui.Button("Send to Plot"), column = 0)
b_useascalib = t_sweep  .place_object(egg.gui.Button("Use as Calibration"), column = 1, alignment = 1)
b_quicksweep = t_sweep  .place_object(egg.gui.Button("Quick Sweep"), column = 1, alignment = 2).set_checkable(True)



# databox plotters for the Ohm tab
plot_fmr      = t_fmr .place_object(egg.gui.DataboxPlot("*.fmr" , "plotsw_fmr.cfg"), column_span=2)
plot_fmr.button_autoscript.set_checked(False)
fmr_script = "x = [ d[0], d[0] ]\ny = [ d[3], d[4] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=\'deltaR (Ohms)\'"; plot_fmr.script.set_text(fmr_script)
t_fmr.new_autorow()
b_send_fmr = t_fmr  .place_object(egg.gui.Button("Send to Plot"), column = 0, alignment = 1)
b_useasref = t_fmr  .place_object(egg.gui.Button("Use as Reference"), column = 1, alignment = 1)
b_resetref = t_fmr  .place_object(egg.gui.Button("Reset Reference"), column = 1, alignment = 2)



# databox plotters and goodies for the IV curve tab
plot_fourpt   = t_rawiv   .place_object(egg.gui.DataboxPlot("*.fourpt" , "plotsw_raw.cfg"), column_span=2)
fourpt_script = "x = [ d[0], d[0] ]\ny = [ d[1], d[2] ]"
plot_fourpt.script.set_text(fourpt_script)

plot_iv       = t_iv  .place_object(egg.gui.DataboxPlot("*.iv" , "plotsw_iv.cfg"), column_span=2)
plot_iv.button_autoscript.set_checked(False)
iv_script = "x = [ d[0], d[0], d[0] ]\ny = [ d[1], d[4], d[5] ]\n\nxlabels=\'Current (A)\'\nylabels=['V', \'dV_amp/dI (Ohms)\', \'dV_nonamp/dI (Ohms)\']"
plot_iv.script.set_text(iv_script)
t_iv.new_autorow()
b_send_iv     = t_iv  .place_object(egg.gui.Button("Send to Plot"), column = 0)
l_imp_iv      = t_iv  .place_object(egg.gui.Label("Offset Voltage = ---mV | Mean Resistance = ---Ohm"), column = 1, alignment=2)



# FMR vs I sweep waterfall plots
# Plot tabs
t_isweep.new_autorow()
plot_isweep = t_isweep.place_object(ColorWFPlot(size = [cfg_sweep['Steps'], cfg_isweep['Steps']], \
    bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_isweep['Start'], cfg_isweep['Stop']]))
plot_isweep.set_labels(['Frequency (Hz)', r'$I_{DC}$ (mA)', r'$\delta R (\Omega)$'])
plot_isweep.add_databox('dR')
plot_isweep.add_databox('dR_c')
plot_isweep.add_databox('std_dR')
plot_isweep.add_databox('std_dR_c')
plot_isweep.set_databox(1)



# Special plots and parameters for the B-field sweep
tabs_bsweep   = t_bsweep  .place_object(egg.gui.TabArea(False), 1,0, row_span=4, alignment=0)
t_bsweep_traj = tabs_bsweep.add_tab("Trajectory")
plot_traj     = t_bsweep_traj.place_object(egg.gui.DataboxPlot("*.traj" , "plotsw_traj.cfg"))

# Results tab triggering the sweep and showing the histogram and waterfall
t_bsweep_res   = tabs_bsweep.add_tab("Results")

plot_bsweep = t_bsweep_res.place_object(ColorWFPlot(size = [cfg_sweep['Steps'], cfg_bsweep['Steps']], \
    bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_bsweep['Start'], cfg_bsweep['Stop']]))
plot_bsweep.set_labels(['Frequency (Hz)', 'Magnetic Field (T)', r'$\delta R (\Omega)$'])
plot_bsweep.add_databox('dR')
plot_bsweep.add_databox('dR_c')
plot_bsweep.add_databox('std_dR')
plot_bsweep.add_databox('std_dR_c')
plot_bsweep.set_databox(1)
tabs_bsweep.set_current_tab(1)



# Test tab  
plot_test = t_test.place_object(ColorWFPlot(size = [cfg_sweep['Steps'], cfg_bsweep['Steps']], \
    bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_bsweep['Start'], cfg_bsweep['Stop']]))
plot_test.set_labels(['Frequency (Hz)', 'Magnetic Field (T)', r'$\delta R (\Omega)$'])




# Performance labels
window.new_autorow()
g_perf = window.place_object(egg.gui.GridLayout((0,4,2,1)), column_span=2)
l_perf_acqtime   = g_perf.place_object(egg.gui.Label("Acquisition time: --s"), column = 0)
l_perf_dutycyc   = g_perf.place_object(egg.gui.Label("Duty cycle: --%"), column = 1)
l_perf_breakdown = g_perf.place_object(egg.gui.Label("ai Start: --s | ao Start: --s | ao Wait: --s | ai Read: --s | Analysis: --s"), column = 2, alignment = 2)
g_perf.place_object(egg.gui.Label("|||||"), column = 3, alignment=2)
l_rliac_fmr      = g_perf  .place_object(egg.gui.Label("Device resistance = -----Ohm, AC current = -----mA, eta = -----Ohm/mA^2"), column = 4, alignment=2)



# set the default gui stuff
tabs_plots.set_current_tab(0)
plot_raw.button_enabled.set_checked(False)





#####################
### GUI FUNCTIONALITY
#####################

def button_stream_pressed(*a):
    '''
    Called whenever the stream button associated to the camera interface is pressed.
    '''
    # let the loop shut itself down
    if not button_stream.is_checked(): 
        return
    
    channel = int(video_input.get_value())
    
    # connect to the camera
    cam = _cv2.VideoCapture(channel)
    
    # for the frames per second calculation    
    #t0 = time.time()
    n  = 0
    
    # global variables for script execution
    g = _n.__dict__
    

    # loop until we're told not to
    while button_stream.is_checked():
        
        # get an image
        success, image = cam.read()
        if success:
            
            # process the image
            image = camgui.image_to_rgb_data(image)
            g.update(dict(r=image[:,:,0],
                          g=image[:,:,1],
                          b=image[:,:,2]))

            # the try/except thing here prevents a bad script from
            # pooping out the program. 
            try:
                # get the plot data based on the script
                #data = eval(text_script.get_text(), g)                
                data = eval("(r+g+b)/3.0", g)

                image_raw.set_data(data)                          

            # If the script pooped. Quietly do nothing. The 
            # script box should be pink already
            except: 
                pass
                
            # update the frames per second
            n = n+1
            if n%10 == 0:
                #label_fps.set_text("FPS = " + str(int(1.0*n/(time.time()-t0))))
                n = 0
        #        t0 = time.time()
        
        # let the gui update every frame. Otherwise it freezes!    
        window.process_events()

    # release the camera
    cam.release()   
    
# connect the button to the function
window.connect(button_stream.signal_clicked, button_stream_pressed)


def b_send_sweep_clicked(*a):
    """
    Sends the current sweep to a pylab plot.
    """  
    plt.figure(1)
    if cfg_sweep["Experiment"] == "Driven Dynamics":
        FMRfit.fit_spectrum(plot_sweep['f'], plot_sweep['V_sinc'], gcf=True)
        
    #    # In-phase plot
    #    plt.title("Rectified Voltage vs. Frequency (In-Phase Component)")
    #    plt.xlabel('Frequency (GHz)'); plt.ylabel(r'$\Delta V$')
    #    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)); plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #    plt.plot(plot_sweep['f'], plot_sweep['V_sinc'])
    #    plt.tight_layout()
        
        # Out-of-phase plot
        plt.figure(2)
        plt.title("Rectified Voltage vs. Frequency (Out-of-Phase Component)")
        plt.xlabel('Frequency (GHz)'); plt.ylabel(r'$\Delta V$')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)); plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.plot(plot_sweep['f']/1e9, plot_sweep['V_cosc'])
        plt.tight_layout()

    elif cfg_sweep["Experiment"] == "Spectrum Analyzer":
        plt.title("Mean Voltage vs. Frequency")        
        plt.xlabel('Frequency (GHz)'); plt.ylabel(r'$\Delta V$')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)); plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.plot(plot_sweep['f']/1e9, plot_sweep['V'])
        plt.tight_layout()        
        
        # RMS value plot        
        plt.figure(2)
        plt.title("Detected Power vs. Frequency")
        plt.xlabel('Frequency (GHz)'); plt.ylabel('Power (mW)')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)); plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.plot(plot_sweep['f']/1e9, plot_sweep['P_out'])
        plt.tight_layout()        

    plt.show(1)
       
b_send_sweep.signal_clicked.connect(b_send_sweep_clicked)


def b_send_fmr_clicked(*a):
    """
    Sends the current FMR curve to a pylab plot.
    """
    # FMR Plot
    plt.figure(3); plt.hold(True)
    if cfg_sweep['Experiment'] == "Driven Dynamics":
        FMRfit.fit_spectrum(plot_fmr['f'], plot_fmr['dR_c'], sigma = plot_fmr['std_dR_c'], subplot=False, gcf=True)
    elif cfg_sweep['Experiment'] == "Spectrum Analyzer":
        plt.xlabel('Frequency (GHz)'); plt.ylabel('Power (mW)')
        plt.errorbar(plot_fmr['f']/1e9, plot_fmr['P_out'], plot_fmr['std_dR'])
    plt.tight_layout()

b_send_fmr.signal_clicked.connect(b_send_fmr_clicked)


def b_send_iv_clicked(*a):
    """
    Sends the current IV / dV/dI curve to a pylab plot.
    """
    # V vs I plot
    plt.figure(5)
    plt.xlabel(r'$I_{dc}$' + ' (A)'); plt.ylabel(r'$V$')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)); plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(plot_iv[0], plot_iv[1])
    plt.tight_layout()
    
    # dV/dI vs I plot
    plt.figure(6)
    plt.xlabel(r'$I_{dc}$' + ' (A)'); plt.ylabel(r'$dV/dI$' + ' ' + r'$(\Omega)$')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0));
    plt.scatter(plot_iv[0], plot_iv[4])
    # Fit    
    plt.hold(True)
    popt, perr = fit_dVdI(plot_iv[0], plot_iv[4])
    xvals = _n.linspace(1.1*_n.min(plot_iv[0]), 1.1*_n.max(plot_iv[0]), 100)
    yvals = map(lambda x: popt[0] + popt[1]*x**2, xvals)
    plt.plot(xvals,yvals)
    plt.xlim(1.1*_n.min(plot_iv[0]), 1.1*_n.max(plot_iv[0]))
    plt.tight_layout()
    
    # dV/dIn vs I plot
    plt.figure(7)
    plt.xlabel(r'$I_{dc}$' + ' (A)'); plt.ylabel(r'$dV/dI$n' + ' ' + r'$(\Omega)$')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0));
    plt.scatter(plot_iv[0], plot_iv[5])
    # Fit
    plt.hold(True)
    popt, perr = fit_dVdI(plot_iv[0], plot_iv[5])
    xvals = _n.linspace(1.1*_n.min(plot_iv[0]), 1.1*_n.max(plot_iv[0]), 100)
    yvals = map(lambda x: popt[0] + popt[1]*x**2, xvals)
    plt.plot(xvals,yvals)
    plt.xlim(1.1*_n.min(plot_iv[0]), 1.1*_n.max(plot_iv[0]))
    plt.tight_layout()    
    
b_send_iv.signal_clicked.connect(b_send_iv_clicked)


#def b_multiplot_clicked(*a):
#    """
#    Plots a bunch of databoxes.
#    """
#    s.plot.xy.files(xlabel = 'Frequency (Hz)', ylabel = r'$\Delta V$', label = '', shell_history=0, yshift = 0.1)
#    
#b_multiplot.signal_clicked.connect(b_multiplot_clicked)

def b_calibload_clicked(*a):
    """
    Loads the power calibration in the global variables.
    """
    if cfg_sweep['Experiment'] == "Driven Dynamics":
        global f_calib, p_calib, p_factor, p_factor_mixed
        f_calib = plot_calib['f']
        p_calib = plot_calib['P']
        p_factor = plot_calib['P_f']
        p_factor_mixed = plot_calib['P_f_m']
        
        global R_L, eta
        R_L = plot_calib.headers['R_L']
        eta = plot_calib.headers['eta']
    
#    elif cfg_sweep['Experiment'] == "Spectrum Analyzer":
#       global f_calib_sa      
#        f_calib_sa = _s.interpolate.interp2d(plot_calib['f'], plot_calib['P_out'], 10**(plot_calib['P_in']/10))
#        global sa_p_out, sa_gain
        
    print "Calibration loaded"
    
plot_calib.after_load_file = b_calibload_clicked

def b_calibreset_clicked(*a):
    """
    Resets the power calibration.
    """
    global f_calib, p_calib, p_factor, p_factor_mixed
    f_calib, p_calib, p_factor = power_calibration()
    p_factor_mixed = [1]*len(p_factor)
    
    plot_calib.plot()
    reinit_calib_databoxes()
    print "Calibration reset."    
    
b_calibreset.signal_clicked.connect(b_calibreset_clicked)


def b_useascalib_clicked(*a):
    """
    Uses the current sweep data to rescale the power calibration
    """
    global f_calib, p_calib, p_factor 
    f_calib, p_c2, p_f2 = power_calibration(freqs,plot_sweep['V_sinc'])    
    p_factor = _n.multiply(p_factor, p_f2)
    p_calib  = -10*_n.log10(p_factor)
    
b_useascalib.signal_clicked.connect(b_useascalib_clicked)


def b_resetref_clicked(*a):
    """
    Resets the FMR reference spectrum to zero
    """
    global fmr_ref
    fmr_ref = [0]*cfg_sweep['Steps']
    
b_resetref.signal_clicked.connect(b_resetref_clicked)    


def b_useasref_clicked(*a):
    """
    Uses the current FMR data as the new background spectrum
    """
    global fmr_ref
    fmr_ref += -plot_fmr['dR_c']
    
b_useasref.signal_clicked.connect(b_useasref_clicked)


def b_bsweepfile_clicked(*a):
    """
    Loads a file dialog to select the magnet trajectory file
    """
    root = Tk()
    root.withdraw()
    cfg_bsweep["TrajectoryFile"] = filedialog.askopenfilename(filetypes = [("Comma-Separated Values", ".csv")])
    
    
b_bsweepfile.signal_clicked.connect(b_bsweepfile_clicked)


def b_labviewvi_clicked(*a):
    """
    Selects the path of the Labview Host Control VI to interface with.
    """
    root = Tk()
    root.withdraw()    
    cfg_hardware["Labview/VI_Path"] = filedialog.askopenfilename(filetypes = [("Labview VIs", ".vi")]) 
    
b_labviewvi.signal_clicked.connect(b_labviewvi_clicked)



# Initialize frequency vector
if cfg_sweep['Mode'] in ['LIN'] :
    freqs  = _n.linspace(cfg_sweep['Start'], cfg_sweep['Stop'] , cfg_sweep['Steps'])
else:
    logbounds = [_n.log10(cfg_sweep['Start']), _n.log10(cfg_sweep['Stop'])]
    freqs  = _n.logspace(logbounds[0], logbounds[1], cfg_sweep['Steps']) 

# Initialize offset mixdown voltage
V_off = 0.0


def power_calibration(freqs = freqs, V_calib = [1]*cfg_sweep['Steps']):
    """
    Uses sweep data to perform a power calibration.
    """
    # Get factor by which we have to multiply the power
    p_factor = _n.array(V_calib)/V_calib[0]
    p_factor[ _n.where(p_factor <= 0) ] = 1 # Just to prevent invalid values

    # Get power offset on a decibel scale    
    p_calib = -10*_n.log10(p_factor)

    # Also return the calibration frequencies        
    return freqs, p_calib, p_factor

# Run the function once to get some initial data (but first, get the frequency vector)
f_calib, p_calib, p_factor = power_calibration()
p_factor_mixed = [1]*len(p_factor)


def get_Iac(P_ac = cfg_sweep['Power']):
    """
    Calculates the AC current going through the device
    """
    global R_L, I_ac
    I_ac = _n.sqrt(4 * 1e-3*10**(pows[0]/10) * 50) / (R_L + 50)
    
# Load some initial values
R_L = 0.0; I_ac = 1.0; eta = 1.0




################################################
### OUTPUT WAVEFORM AND DATABOXES INITIALIZATION
################################################

def init_pulm(ao_t_sw):
    """
    Builds the analog output square wave vector for the pulse modulation
    """
    tsteps = int(cfg_sweep['PULM/Period'] * cfg_sweep['PULM/Repeats'] * cfg_hardware['AO/Rate'])
    
    ao_sw = _n.zeros(tsteps)  
    
    if cfg_sweep["Experiment"] == "Driven Dynamics":
        for n in range(0,int(tsteps)) :
            if _n.ceil(2*ao_t_sw[n]/(cfg_sweep['PULM/Period'])) % 2 == 0 \
            and n != 0: # first and last points are on the low state
                ao_sw[n] = cfg_sweep['PULM/Amplitude'] + cfg_sweep['PULM/DC']
            else:
                ao_sw[n] = cfg_sweep['PULM/DC']
            
    return ao_sw


def init_dc(cur):
    """
    Builds the analog output vector for an applied DC current during the sweep.
    The arguments are the applied current in amperes and the device resistance.
    """ 
    global R_L
    
    # Build time steps vector
    tsteps = int(cfg_sweep['PULM/Period'] * cfg_sweep['PULM/Repeats'] * cfg_hardware['AO/Rate'])  
    
    # Build DC currents vector. Finish off with zero current.
    dc_volt = cur * (R_L + cfg_hardware["DCSourceR"])
    ao_dc = _n.zeros(tsteps) + dc_volt

    return ao_dc


def init_fourpt(cur, mod=0):
    """
    Builds the DC output square wave vector for four-point measurement.
    The arguments are the nominal applied current and the modulation (both in amperes)
    """   
    global R_L
    
    # Compute time steps vector and number of time steps
    t_vec   = _n.arange(0.0,cfg_fourpt["ModPeriod"] * cfg_fourpt["ModRepeats"], 1.0/cfg_hardware['AO/Rate'])
    t_steps = len(t_vec)
        
    # Calculate corresponding voltage and build array. Finish off with a current of 0.
    V_dc = cur * (R_L + cfg_hardware["DCSourceR"])
    V_ac = mod * (R_L + cfg_hardware["DCSourceR"])
    
    # Calculate cosine and sine arrays (used in lock-in analysis of dV/dI)
    cos_iv = _n.cos(2*_n.pi*t_vec / cfg_fourpt["ModPeriod"])
    sin_iv = _n.sin(2*_n.pi*t_vec / cfg_fourpt["ModPeriod"])
    
    ao_fourpt = _n.zeros(t_steps) + V_dc + V_ac*sin_iv
    
    return ao_fourpt, cos_iv, sin_iv


def init_sweep_arrays(dc_cur = cfg_sweep['DC_Current']):
    """
    Initializes the various waveforms and arrays required for the frequency sweep acquisition
    and analysis
    """  
    global ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr, R_L
    
    # Build the time steps vector
    ao_t = _n.linspace(1/cfg_hardware['AO/Rate'], cfg_sweep['PULM/Period']*cfg_sweep['PULM/Repeats'], 1.0*cfg_sweep['PULM/Period']*cfg_sweep['PULM/Repeats']*cfg_hardware['AO/Rate'])

    # Build square wave for PULM
    ao_sw = init_pulm(ao_t)
    
    # Build DC current and four-point measurement arrays
    ao_dc = init_dc(dc_cur)
          
    # Plot new analog output curve
    #ao_plot.clear()    
    #ao_curve = egg.pyqtgraph.PlotCurveItem(ao_t, ao_sw, stepMode=False, fillLevel=0)
    #ao_plot.addItem(ao_curve)
    #ao_plot.autoRange()
        
    # Prepare sine and cosine array for lock-in analysis
    cos_li = _n.cos(2*_n.pi*_n.array(ao_t)/(cfg_sweep['PULM/Period']))
    sin_li = _n.sin(2*_n.pi*_n.array(ao_t)/(cfg_sweep['PULM/Period']))        
    
    # Build frequency vector
    if cfg_sweep['Mode'] in ['LIN'] :
        freqs  = _n.linspace(cfg_sweep['Start'], cfg_sweep['Stop'] , cfg_sweep['Steps'])
    else:
        logbounds = [_n.log10(cfg_sweep['Start']), _n.log10(cfg_sweep['Stop'])]
        freqs  = _n.logspace(logbounds[0], logbounds[1], cfg_sweep['Steps']) 
        
    # Build power vector (accounting for calibration and saturation)
    if cfg_sweep['PowerMode'] == 'Constant':
        pows = [cfg_sweep['Power']]*len(freqs)
        pows_corr = _n.interp(freqs, f_calib, p_factor)**-1.0
    elif cfg_sweep['PowerMode'] == 'Calibrated':
        #print "I didn't break everything" #Adrian, it doesn't seem to always load p_factor properly
        p_precalib  = -10*_n.log10(_n.multiply(p_factor, p_factor_mixed))
        pows = _n.interp(freqs, f_calib, p_precalib) + cfg_sweep['Power']
        if max(pows) > cfg_hardware['SaturationPow']: #Edited to include SatPow, Adrian
           print "Power saturated (max calibrated power = " + str(max(pows)) + "dBm)" 
           pows[ _n.where(pows>cfg_hardware['SaturationPow']) ] = cfg_hardware['SaturationPow'] #Adrian
        pows_corr = _n.array([1]*len(freqs))**-1.0
    elif cfg_sweep['PowerMode'] == 'Mixed':
        pows = _n.interp(freqs, f_calib, p_calib) + cfg_sweep['Power']
        pows_corr = _n.interp(freqs, f_calib, p_factor_mixed)**-1.0
    
    return ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr
    
# Do it once when the software is initialized    
init_sweep_arrays()


def init_iv_arrays():
    """
    Initializes the various waveforms and arrays required for the IV curve acquisition
    and analysis
    """ 
    global ao_fourpt, cos_iv, sin_iv, curs
    
    ao_fourpt, cos_iv, sin_iv = init_fourpt(cfg_fourpt["Start"], cfg_fourpt["ModAmp"])    
    
    # Build current vector
    if cfg_fourpt["Mode"] in ['LIN'] :
        curs  = _n.linspace(cfg_fourpt["Start"], cfg_fourpt["Stop"] , cfg_fourpt["Steps"])
    else:
        logbounds = [_n.log10(cfg_fourpt["Start"]), _n.log10(cfg_fourpt["Stop"])]
        curs  = _n.logspace(logbounds[0], logbounds[1], cfg_fourpt["Steps"]) 
    
    curs = _n.concatenate((curs,curs[::-1]),1)
    curs = _n.roll(curs,len(curs)/4+1)[-1::-1]
    curs = _n.concatenate(([0], curs), 1)
    
    return ao_fourpt, cos_iv, sin_iv, curs
    
# Do it once when the software is initialized    
init_iv_arrays()



def init_bsweep_arrays():
    """
    Initializes the magnet path and B-field list
    """ 
    global bfields, conex_path

    # Build the magnetic field list
    if cfg_bsweep["Mode"] in ['LIN'] :
        bfields  = _n.linspace(cfg_bsweep["Start"], cfg_bsweep["Stop"] , cfg_bsweep["Steps"])
    else:
        logbounds = [_n.log10(cfg_bsweep["Start"]), _n.log10(cfg_bsweep["Stop"])]
        bfields  = _n.logspace(logbounds[0], logbounds[1], cfg_bsweep["Steps"])     
    
    # Build the path needed to create the aformentioned B-field list
    conex_path = []
    conex_path = _n.zeros( (cfg_bsweep["Steps"],3) )   
    try:             
        # Import the trajectory file giving the radial and depth coordinates of the magnet 
        # for equally stepped B-fields
        csvfile = open(cfg_bsweep["TrajectoryFile"], 'r')
        csvpath = csv.reader(csvfile)   
       
        magnet_path = []
        for row in csvpath:
            magnet_path.append(map(float,row))            
        magnet_path = _n.flipud(_n.array(magnet_path))
        
        # Transform the imported array in a set of xyz coordinates taking into account
        # the desired radial angle and the xyz offsets of the centre of the magnet WRT
        # to the sample, for each B-field step
        
        # Interpolate the magnet path, rotate to the desired radial angle, add offsets       
        conex_path[:,0] = _n.interp(bfields, magnet_path[:,2], magnet_path[:,0]) * _n.cos(cfg_bsweep["Offsets/Alpha"]) + cfg_bsweep["Offsets/x0"]
        conex_path[:,1] = _n.interp(bfields, magnet_path[:,2], magnet_path[:,0]) * _n.sin(cfg_bsweep["Offsets/Alpha"]) + cfg_bsweep["Offsets/y0"]
        conex_path[:,2] = _n.interp(bfields, magnet_path[:,2], magnet_path[:,1])                                       + cfg_bsweep["Offsets/z0"]

    except IOError:
        if cfg_bsweep["TrajectoryFile"] != "": print "No such file or directory: " + cfg_bsweep["TrajectoryFile"]
        
# Do it once when the software is initialized    
init_bsweep_arrays()


def init_isweep_arrays():
    """
    Initializes the current arrays for an FMR vs I sweep
    """ 
    global currentvals    
    currentvals = _n.linspace(cfg_isweep["Start"], cfg_isweep["Stop"], cfg_isweep["Steps"])

# Do it once when the software is initialized    
init_isweep_arrays()


def init_sa_calib_arrays(freq_ref = 0, power = -100):
    """
    Initializes the various waveforms and arrays required for the spectrum analyzer
    power calibration curve acquisition and analysis
    """ 
    global freqs_calib, pows_calib
    
    freqs_calib = _n.linspace(freq_ref-cfg_pcalib["SA_Calib/Range"], freq_ref+cfg_pcalib["SA_Calib/Range"], cfg_pcalib["SA_Calib/Steps"])
    nsteps = len(freqs_calib)    
    pows_calib  = nsteps*[power]

# Do it once when the software is initialized    
init_sa_calib_arrays()


def reinit_sweep_databoxes():
    """
    Reinitializes the sweep databoxes for a new acquisition.
    """    
    ### clear out the old data, reinitialize databoxes
    plot_raw.clear()
        
    # Reinitialize raw databox
    plot_raw['t']    = ao_t
    plot_raw["V0"]   = [None]*len(ao_t)
    #plot_raw["PULM"] = [None]*len(ao_t)

    if cfg_sweep["Experiment"] == "Driven Dynamics":
        # Reinitialize lock-in databox
        plot_sweep['f']      = []
        plot_sweep['P']      = []
        plot_sweep['P_f']    = []
        plot_sweep['V_sin']  = []
        plot_sweep['V_cos']  = []   
        plot_sweep['V_sinc'] = []
        plot_sweep['V_cosc'] = []

    if cfg_sweep["Experiment"] == "Spectrum Analyzer":
        plot_sweep['f']      = []
        plot_sweep['V']      = []
        plot_sweep['V_rms']  = []
        plot_sweep['P_out']  = []
 

def reinit_fmr_databoxes():
    """
    Reinitializes the preallocated FMR databoxes for a new acquisition.
    """              
    if cfg_sweep["Experiment"] == "Driven Dynamics":
        plot_fmr['f']  = freqs
        plot_fmr['P']  = pows
        if cfg_sweep['PowerMode'] == 'Constant':
            plot_fmr['I_ac'] = plot_calib['I_f'] * 10**(cfg_sweep['Power']/20.0) 
        else:
            plot_fmr['I_ac'] = plot_calib['I_f'] * 10**(cfg_sweep['Power']/20.0) / _n.sqrt(p_factor)
        plot_fmr['dR'] = [0]*cfg_sweep['Steps']
        plot_fmr['dR_c'] = [0]*cfg_sweep['Steps']   
        plot_fmr['std_dR'] = [0]*cfg_sweep['Steps']
        plot_fmr['std_dR_c'] = [0]*cfg_sweep['Steps']
    elif cfg_sweep["Experiment"] == "Spectrum Analyzer":
        plot_fmr['f']  = freqs
        plot_fmr['P_out'] = [0]*cfg_sweep['Steps']
        plot_fmr['P_in']  = [0]*cfg_sweep['Steps']
        plot_fmr['std_dR'] = [0]*cfg_sweep['Steps']
        plot_fmr['std_dR_c'] = [0]*cfg_sweep['Steps']     

def reinit_iv_databoxes():  
    """
    Reinitializes the IV curve databoxes for a new acquisition.
    """       
    ### clear out the old data, reinitialize databoxes
    plot_fourpt.clear()   
    plot_iv.clear()      
    
    # Reinitialize raw databox
    plot_iv['I']        = []
    plot_iv['V']        = []
    plot_iv['dVdI_cos'] = []
    plot_iv['dVdI_sin'] = []
    plot_iv['dVdI']     = []
    plot_iv['dVdI_namp']    = []
    #dVdI_phi  = [0]*2*cfg_fourpt["Steps"]
    
    
def reinit_calib_databoxes():
    """
    Reinitializes the preallocated power calibration databoxes for a new acquisition.    
    """
    if cfg_sweep["Experiment"] == "Driven Dynamics":
        # Change plot scripts
        sweep_script = "x = [ d[0], d[0] ]\ny = [ d[5], d[6] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=\'V\'"
        plot_sweep.script.set_text(sweep_script)
        
        calib_script = "x = [ d[0], d[0], d[0] ]\ny = [ d[1], d[3], d[4] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=[\'P_corr (dBm)\', \'P_fmix\', \'I_f (mA)\']"
        plot_calib.script.set_text(calib_script)
        plot_calib['f']     = freqs
        plot_calib['P']     = [cfg_sweep['Power']]*cfg_sweep['Steps']
        plot_calib['P_f']   = [1]*cfg_sweep['Steps']
        plot_calib['P_f_m'] = [1]*cfg_sweep['Steps']     
        plot_calib['I_f']   = [1]*cfg_sweep['Steps']
        
    if cfg_sweep["Experiment"] == "Spectrum Analyzer":
        # Change plot scripts
        sweep_script = "x = [ d[0], d[0] ]\ny = [ d[1], d[3] ]"
        plot_sweep.script.set_text(sweep_script)

        calib_script = "x = [ d[1], d[1] ]\ny = [ d[2], d[3] ]"
        plot_calib.script.set_text(calib_script)
        plot_calib['P_in']   = [] 
        plot_calib['f']      = []
        plot_calib['P_out']  = [] 
        plot_calib['Gain']   = [] 
        plot_calib['V_base'] = []
        
        global f_calib_sa
        f_calib_sa = lambda x,y: y
        
def reinit_bsweep_databoxes():
    """
    Reinitializes the B-field sweep databoxes for a new acquisition.    
    """
    plot_bsweep.set_size_and_bounds(size = [cfg_sweep['Steps'], cfg_bsweep['Steps']], \
        bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_bsweep['Start'], cfg_bsweep['Stop']])
    plot_bsweep.clear_databoxes()
        
def reinit_isweep_databoxes():
    """
    Reinitializes the DC current sweep databoxes for a new acquisition.    
    """
    plot_isweep.set_size_and_bounds(size = [cfg_sweep['Steps'], cfg_isweep['Steps']], \
        bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_isweep['Start'], cfg_isweep['Stop']])
    plot_isweep.clear_databoxes()

def cfg_sweep_changed(*b):
    """
    Called whenever we change the settings Dictionary.
    """
    global a, freqs, pows
    
    # Reinitialize acquisition arrays (mostly for developing / troubleshooting)
    init_sweep_arrays()
    init_iv_arrays()        
    
    #a.list_man_setup(freqs,pows)
      
    ##### save the settings for the next boot-up
    cfg_sweep.save()

# load some initial settings (will trigger the settings changed)
# or just fire the settings changed to take care of initialization
cfg_sweep.load()

# connect this function and get some initial settings
cfg_sweep.connect_any_signal_changed(cfg_sweep_changed)
cfg_sweep_changed()


def cfg_fourpt_changed(*b):
    """
    Called whenever we change the settings Dictionary.
    """
    global a, freqs, pows
    
    # Reinitialize acquisition arrays (mostly for developing / troubleshooting)
    init_iv_arrays()        
      
    ##### save the settings for the next boot-up
    cfg_fourpt.save()

# load some initial settings (will trigger the settings changed)
# or just fire the settings changed to take care of initialization
cfg_fourpt.load()
cfg_fourpt.connect_any_signal_changed(cfg_fourpt_changed)
cfg_fourpt_changed()


def cfg_pcalib_changed(*b):
    """
    Called whenever we change the settings Dictionary.
    """
    global a, freqs, pows
    
    # Reinitialize acquisition arrays (mostly for developing / troubleshooting)
    init_sweep_arrays()        
      
    ##### save the settings for the next boot-up
    cfg_pcalib.save()

# load some initial settings (will trigger the settings changed)
# or just fire the settings changed to take care of initialization
cfg_pcalib.load()
cfg_pcalib.connect_any_signal_changed(cfg_pcalib_changed)
cfg_pcalib_changed()



def cfg_bsweep_changed(*b):
    """
    Called whenever we change the settings Dictionary.
    """
    global bfields, conex_path
    
    # Reinitialize acquisition arrays (mostly for developing / troubleshooting)
    init_bsweep_arrays()        

    # Update trajectory plots
    plot_traj['B'] = bfields
    plot_traj['X'] = conex_path[:,0]
    plot_traj['Y'] = conex_path[:,1]
    plot_traj['Z'] = conex_path[:,2]
    
    plot_traj.plot()
      
    ##### save the settings for the next boot-up
    cfg_bsweep.save()
   

# load some initial settings (will trigger the settings changed)
# or just fire the settings changed to take care of initialization
cfg_bsweep.load()
cfg_bsweep.connect_any_signal_changed(cfg_bsweep_changed)
cfg_bsweep_changed()
    
def cfg_isweep_changed(*b):
    """
    Called whenever we change the settings Dictionary.
    """
    global currentvals
    
    # Reinitialize acquisition arrays
    init_isweep_arrays()
    
    ##### save the settings for the next boot-up
    cfg_isweep.save()

# load some initial settings (will trigger the settings changed)
# or just fire the settings changed to take care of initialization
cfg_isweep.load()
cfg_isweep.connect_any_signal_changed(cfg_isweep_changed)
#print "who needs pants"
cfg_isweep_changed()

def cfg_hardware_changed(*a):
    """
    Called whenever we change the settings Dictionary.
    """
    # Reinitialize acquisition arrays (mostly for developing / troubleshooting)
    init_sweep_arrays()    

    ##### save the settings for the next boot-up
    cfg_hardware.save()

# load some initial settings (will trigger the settings changed)
cfg_hardware.load()
     
# connect this function and get some initial settings
cfg_hardware.connect_any_signal_changed(cfg_hardware_changed)
cfg_hardware_changed()


def RL_changed(*a):
    """
    Called whenever we change a parameter that changes the device resistance
    """
    global R_L
    R_L = 0
    l_rliac_fmr.set_text("Device resistance = -----Ohm, AC current = -----mA")
#cfg_sweep.connect_signal_changed('DC_Current', RL_changed)


def Iac_changed(*a):
    """
    Called whenever the RF power is changed.
    """
    if I_ac != 1.0:
        get_Iac()
        a = l_rliac_fmr.get_text()
        b = a.split("AC current = ")[0]
        l_rliac_fmr.set_text(b + "AC current = " + str(round(I_ac*100000)/100) + "mA")
cfg_sweep.connect_signal_changed('Power', Iac_changed)    

    
def after_load_file(p):
    """
    To be called whenever one of the DataboxPlot instances loads a file.
    Used here to update the daq settings whenever we load a file.
    """
    cfg_sweep.update(p.headers, ignore_errors=True)

# overload the existing after_load_file functions
plot_raw    .after_load_file = after_load_file
plot_sweep   .after_load_file = after_load_file


def experiment_change(*a):
    """
    Reinitializes a few settings if we change the experiment
    """
    exp = cfg_sweep["Experiment"]   
    plot_calib.clear()
    
    reinit_calib_databoxes()

    if exp == "Driven Dynamics":
        # Change signal generator level
        cfg_sweep['Power'] = 0   
        
        # Change plot scripts
        sweep_script = "x = [ d[0], d[0] ]\ny = [ d[5], d[6] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=\'V\'"
        plot_sweep.script.set_text(sweep_script)
        
        calib_script = "x = [ d[0], d[0], d[0] ]\ny = [ d[1], d[3], d[4] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=[\'P_corr (dBm)\', \'P_fmix\', \'I_f\']"

        plot_calib.script.set_text(calib_script)
        plot_calib['f']     = freqs
        plot_calib['P']     = [cfg_sweep['Power']]*cfg_sweep['Steps']
        plot_calib['P_f']   = [1]*cfg_sweep['Steps']
        plot_calib['P_f_m'] = [1]*cfg_sweep['Steps']  
                
        fmr_script = "x = [ d[0], d[0] ]\ny = [ d[3], d[4] ]\n\nxlabels=\'Frequency (Hz)\'\nylabels=\'deltaR (Ohms)\'"
        plot_fmr.script.set_text(fmr_script)
        
        # Turn PULM button on
        b_pulm.set_checked(True)
        b_pulm.set_text("PULM On")        
        
    if exp == "Spectrum Analyzer":
        # Change signal generator level
        cfg_sweep['Power'] = 10   
        cfg_sweep['DC_Current'] = 0
        
        # Change plot scripts
        sweep_script = "x = [ d[0], d[0] ]\ny = [ d[1], d[3] ]"
        plot_sweep.script.set_text(sweep_script)

        calib_script = "x = [ d[1], d[1] ]\ny = [ d[2] , d[3] ]"
        plot_calib.script.set_text(calib_script)

        fmr_script   = "x = [ d[0], d[0] ]\ny = [ d[1], d[2] ]"
        plot_fmr.script.set_text(fmr_script)

        # Load the calibration data
        global f_calib_sa
        f_calib_sa = lambda x,y: y
        
        global sa_p_out, sa_gain
        sa_p_out = _n.ones((cfg_pcalib['SA_Calib/PowerSteps'], cfg_sweep['Steps']))
        sa_gain  = _n.ones((cfg_pcalib['SA_Calib/PowerSteps'], cfg_sweep['Steps']))

        # Turn PULM button off
        b_pulm.set_checked(False)
        b_pulm.set_text("PULM Off")
cfg_sweep.connect_signal_changed("Experiment", experiment_change)    

# Initialize experiment-dependent settings
experiment_change()



#########################
### LABVIEW COMMUNICATION
#########################

def labview_import_arrays():
    """
    Imports arrays from Labview for plotting and further analysis in Python.
    """
    labview = win32com.client.Dispatch("Labview.Application")
    vi = labview.getvireference(cfg_hardware["Labview/VI_Path"].replace("/", "\\"))  # Path to LabVIEW VI
    vi._FlagAsMethod("Call")  # Flag "Call" as Method
    xyzcounts = _n.array(vi.getcontrolvalue('XYZCounts'))

    return xyzcounts




##########################
### DAQ I/O TASKS
##########################

def ao_task_pulm():
    """
    Sets up the analog output task for pulse modulation.
    """    
    # Get the channel names
    ao_channel_names = _daqmx.get_ao_channel_names(cfg_hardware['AO/Device'])  

    # Set up the pulse train output    
    ao = _daqmx.ao_task(ao_trigger_source = "/"+cfg_hardware['AO/Device']+"/100kHzTimebase",
                        ao_channels  = [ao_channel_names[cfg_hardware['AO/PULM']], ao_channel_names[cfg_hardware['AO/DC']]],
                        ao_waveforms = [ao_sw, ao_dc],
                        ao_rate      =  cfg_hardware['AO/Rate'],
                        ao_min       = -cfg_hardware['AO/Range'],
                        ao_max       =  cfg_hardware['AO/Range'],
                        ao_StartTrigger_terminal = "/"+cfg_hardware['AO/Device']+"/PFI0")       
                        
    return ao
    

def ao_task_fourpt():
    """
    Sets up the analog output task for the four point measurements
    """   
    # Get the channel names
    ao_channel_names = _daqmx.get_ao_channel_names(cfg_hardware['AO/Device'])      
        
    # Set up the pulse train output    
    ao = _daqmx.ao_task(ao_trigger_source = "/"+cfg_hardware['AO/Device']+"/100kHzTimebase",
                        ao_channels  = [ao_channel_names[cfg_hardware['AO/DC']]],
                        ao_waveforms = [ao_fourpt],
                        ao_rate      =  cfg_hardware['AO/Rate'],
                        ao_min       = -cfg_hardware['AO/Range'],
                        ao_max       =  cfg_hardware['AO/Range'],
                        ao_StartTrigger_terminal = "/"+cfg_hardware['AO/Device']+"/PFI0")                         
    return ao
    
    
    
def ai_task_sweep():
    """
    Sets up the analog input task for sweeps.
    """        
    
    # Choose the appropriate analog input depending on the experiment
    if cfg_sweep["Experiment"] == "Driven Dynamics":
        ai_channels = eval(cfg_hardware['AI/Sweep'])
    elif cfg_sweep["Experiment"] == "Spectrum Analyzer":
        ai_channels = eval(cfg_hardware['AI/Crystal_Detector'])
        
    if not type(ai_channels) in [list, tuple]: ai_channels = [ai_channels]
    ai_channels = list(ai_channels)
    ai_channel_names = _daqmx.get_ai_channel_names(cfg_hardware['AI/Device'])

    # convert the channel numbers to channel names
    selected_ai_channel_names = []
    for n in range(len(ai_channels)): selected_ai_channel_names.append(ai_channel_names[ai_channels[n]])    
        
    # set up the range
    try:    
        ai_max = eval(cfg_hardware['AI/Range'])
    except:
        print "ERROR: "+repr(cfg_hardware['AI/Range'])+" cannot be evaluated."
        b_sweep.set_checked(False)
        return

    # if it's a list
    if hasattr(ai_max, "__len__"):
        # make a list of ai_min's
        ai_min = []
        for x in ai_max: ai_min.append(-x)

    # otherwise it's just a number
    else: ai_min = -ai_max


    # Define analog input task               
    ai = _daqmx.ai_task(ai_trigger_source   = "/"+cfg_hardware['AO/Device']+"/ao/StartTrigger",
                    ai_channels         =  selected_ai_channel_names,
                    ai_samples          =  int(cfg_sweep['PULM/Period']*cfg_sweep['PULM/Repeats']*cfg_hardware['AI/Rate']),
                    ai_rate             =  cfg_hardware['AI/Rate'],
                    ai_input_couplings  =  None,
                    ai_min              =  ai_min,
                    ai_max              =  ai_max) 
                    
    return ai
      

def ai_task_fourpt():
    """
    Sets up the analog input task for a four-point measurement.
    """        
    ai_channels = [ eval(cfg_hardware['AI/FourPT']), eval(cfg_hardware['AI/dVdI']) ]
    if not type(ai_channels) in [list, tuple]: ai_channels = [ai_channels]
    ai_channels = list(ai_channels)
    ai_channel_names = _daqmx.get_ai_channel_names(cfg_hardware['AI/Device'])

    # convert the channel numbers to channel names
    selected_ai_channel_names = []
    for n in range(len(ai_channels)): selected_ai_channel_names.append(ai_channel_names[ai_channels[n]])    

    # set up the range
    try:    
        ai_max = eval(cfg_hardware['AI/Range'])
    except:
        print "ERROR: "+repr(cfg_hardware['AI/Range'])+" cannot be evaluated."
        b_sweep.set_checked(False)
        return

    # if it's a list
    if hasattr(ai_max, "__len__"):
        # make a list of ai_min's
        ai_min = []
        for x in ai_max: ai_min.append(-x)

    # otherwise it's just a number
    else: ai_min = -ai_max
    
    # Define analog input task               
    ai = _daqmx.ai_task(ai_trigger_source   = "/"+cfg_hardware['AO/Device']+"/ao/StartTrigger",
                    ai_channels         =  selected_ai_channel_names,
                    ai_samples          =  int(round(cfg_fourpt["ModPeriod"] * cfg_fourpt["ModRepeats"] *cfg_hardware['AI/Rate'])),
                    ai_rate             =  cfg_hardware['AI/Rate'],
                    ai_min              =  ai_min,
                    ai_max              =  ai_max)     
    
    return ai


def current_ramp(cur_i, cur_f, nsamps = 10000):
    """
    Sets up the analog input task to ramp the DC current from one value to another. The process takes (nsamps/10000)s.
    """        
    global R_L    
  
    ramp_time = _n.abs(cur_i*1000.0-cur_f*1000.0) * cfg_hardware["AO/DCRampTime"]
    
    #t1 = time.clock() 
    if round(_n.abs(cur_i*100000.0-cur_f*100000.0)) != 0: 
        print "Ramping from " + str(round(cur_i*100000)/100) + "mA to " + str(round(cur_f*100000)/100) + "mA (" + str(float(ramp_time*10)/10) + "s)"
    
    # Get the channel names
    ao_channel_names = _daqmx.get_ao_channel_names(cfg_hardware['AO/Device'])      
    ao_rate = cfg_hardware["AO/Rate"]

    ao_ramp = _n.linspace(cur_i * (R_L + cfg_hardware["DCSourceR"]), cur_f * (R_L + cfg_hardware["DCSourceR"]), int(ao_rate*ramp_time)) 
 
    if round(_n.abs(cur_i*100000.0-cur_f*100000.0)) != 0:
        # Set up the pulse train output    
        ao = _daqmx.ao_task(ao_trigger_source = "/"+cfg_hardware['AO/Device']+"/100kHzTimebase",
                            ao_channels  = [ao_channel_names[cfg_hardware['AO/DC']]],
                            ao_waveforms = [ao_ramp],
                            ao_rate      =  ao_rate,
                            ao_min       = -cfg_hardware['AO/Range'],
                            ao_max       =  cfg_hardware['AO/Range'],
                            ao_StartTrigger_terminal = "/"+cfg_hardware['AO/Device']+"/PFI0")       
                            
        ao.start()
        ao.wait_and_clean()
        
    #t2 = time.clock()    
    #print cfg_hardware["AO/DCRampTime"], cfg_hardware["AO/Rate"], str(t2-t1)




########################
### ACQUISITION ROUTINES
########################
def sweep_and_plot():
    """
    Self-explanatory.
    """      
    # Reinitialize databoxes and add offset mixdown voltage in header
    reinit_sweep_databoxes()   
 
    # Insert headers
    plot_sweep.insert_header('Experiment', cfg_sweep['Experiment'])
    plot_sweep.insert_header('SigGen', cfg_hardware['RFSource'])
    plot_sweep.insert_header('AmpGain', cfg_hardware['Amp/Gain'])
    plot_sweep.insert_header('I_dc', cfg_sweep['DC_Current']) 
    plot_sweep.insert_header('PULM_Period', cfg_sweep['PULM/Period'])
    plot_sweep.insert_header('PULM_Repeats', cfg_sweep['PULM/Repeats'])
    plot_sweep.insert_header('R_L', R_L)
    plot_sweep.insert_header('eta', eta)
  
    # Set up analog output and input tasks (and activate pulse modulation)
    ao = ao_task_pulm()
    ai = ai_task_sweep()
    
    # Time the full process
    t_full = _n.array([time.clock(), 0])
    t = _n.zeros([cfg_sweep["Steps"],6])  

    # Loop over frequencies
    for n in range(0,cfg_sweep["Steps"]):     
        
        #if not b_quicksweep.is_checked() and not b_sweep.is_checked(): break
        
        # End acquisition is acquire button is unchecked
        # if not b_sweep.is_checked() and not b_calib.is_checked(): break
        
        # Time the various steps: checkpoint 0
        t[n,0] = time.clock()
        
        # Jump to next list point 
        a.list_change_index(n)
        
        # Settling        
        time.sleep(cfg_hardware['SettlingTime'])
        
        # Start analog input
        ai.start()
 
        # Time the various steps: checkpoint 1
        t[n,1] = time.clock()
        
        # Start analog output            
        ao.start()   

        # Time the various steps: checkpoint 2
        t[n,2] = time.clock()        
        
        # Wait until acquisition ends
        ao.wait_and_clean()
        
        # Time the various steps: checkpoint 3
        t[n,3] = time.clock()        
                
        # retrieve the data
        y = ai.read_and_clean() 
         
        # Time the various steps: checkpoint 4
        t[n,4] = time.clock()
    
        # store the data
        plot_raw["V0"] = y[0]
        #plot_raw["PULM"] = y[1]     
        
        # Plot raw data
        plot_raw.plot()

        # Analyze and store the data
        if cfg_sweep["Experiment"] == "Driven Dynamics":        
            # Perform boxcar analysis        
            V_sin, V_cos = boxcar_analysis(y) 
        
            # Calibrated values
            V_sinc = (V_sin-V_off) * pows_corr[n]
            V_cosc = (V_cos) * pows_corr[n]
        
            if n == 0: 
                plot_sweep.insert_header("V_nonres", V_sinc)
                plot_sweep.insert_header("V_off", V_off)     
            
            # Store into the lock-in databox and plot
            plot_sweep.append_data_point([freqs[n], pows[n], pows_corr[n], V_sin, V_cos, V_sinc, V_cosc])            
            plot_sweep.plot()
            
        elif cfg_sweep["Experiment"] == "Spectrum Analyzer":
            V_mean = _n.mean(y[0])
            V_rms  = _n.sqrt(_n.mean(y[0]**2))
            P_mean = 10**(_n.interp(_n.mean(y[0]), cd_calib[:,1], cd_calib[:,0]) / 10)  # power in milliwatts
            plot_sweep.append_data_point([freqs[n], V_mean, V_rms, P_mean])
            if n == 0: 
                plot_sweep.insert_header("V_off", V_off)     
            plot_sweep.plot()

        # Calculate resistance oscillations and plot
#        plot_fmr.append_data_point([freqs[n], pows[n], 2/(I_ac) * (V_sinc-V_sinc_ref)/cfg_hardware['Amp/Gain']])
#        plot_fmr.plot()
        
        # Time the various steps: checkpoint 5
        t[n,5] = time.clock()
        
        # let the user interface update
        window.process_events()
        
    # Time the full process
    t_full[1] = time.clock()
    duty_cycle_full = cfg_sweep['PULM/Period']*cfg_sweep['PULM/Repeats']*cfg_sweep['Steps']/(t_full[1]-t_full[0])
    
    t_acq = [];
    for n in range(1,6): t_acq.append(round(sum(t[:,n]-t[:,n-1])/cfg_sweep['Steps']*1000)/1000)

    # Update performance labels and plot
    l_perf_acqtime.set_text("Acquisition time: " + str(round((t_full[1]-t_full[0])*1000)/1000) +"s")
    l_perf_dutycyc.set_text("Duty cycle: " + str(round((duty_cycle_full)*100)) +"%")
    breakdown_str = "ai Start: " + str(t_acq[0]) + "s | " + \
                    "ao Start: " + str(t_acq[1]) + "s | " + \
                    "ao Wait: "  + str(t_acq[2]) + "s | " + \
                    "ai Read: "  + str(t_acq[3]) + "s | " + \
                    "Analysis: " + str(t_acq[4]) + "s"
    l_perf_breakdown.set_text(breakdown_str)


def iv_and_plot():
    """
    Performs a four-point measurement on the device.
    """
    global ao_fourpt, cos_iv, sin_iv, curs, R_L, eta
    
    reinit_iv_databoxes()

    # Insert headers
    plot_iv.insert_header('Experiment', cfg_sweep['Experiment'])
    plot_iv.insert_header('SigGen', cfg_hardware['RFSource'])
    plot_iv.insert_header('AmpGain', cfg_hardware['Amp/Gain'])
    plot_iv.insert_header('ModPeriod', cfg_fourpt['ModPeriod'])
    plot_iv.insert_header('ModRepeats', cfg_fourpt['ModRepeats'])
    plot_iv.insert_header('ModAmp', cfg_fourpt['ModAmp'])
    
    # Measure the device resistance
    R_L, offset = get_RL()

    ### Actual measurement
    # Forward loop over frequencies
    for n in range(0,2*cfg_fourpt["Steps"]+1):  
        if not b_fourpt.is_checked(): break;
    
        # Initialize acquisition and analysis arrays
        ao_fourpt, cos_iv, sin_iv = init_fourpt(curs[n], cfg_fourpt["ModAmp"])
        
        if n == 0:
            current_ramp(0,curs[0])
        else:
            current_ramp(curs[n-1],curs[n])
        time.sleep(1)
        
        # Prepare analog output and input tasks
        ao = ao_task_fourpt()
        ai = ai_task_fourpt()
    
        # Start analog input
        ai.start()
        
        # Start analog output
        ao.start()
    
        # Wait until acquisition ends
        ao.wait_and_clean()
        
        # retrieve the data
        y = ai.read_and_clean()
        
        # Get mean resistance
        V_iv, dVdI_cos, dVdI_sin, dVdI, dVdI_namp = fourpt_analysis(y, curs[n], offset)
        
        # Plot and store
        plot_fourpt['t'] = list(_n.arange(0, len(ao_fourpt)) / cfg_hardware['AI/Rate'])
        plot_fourpt['V'] = y[0]
        plot_fourpt['dVdI'] = y[1]
        
        plot_fourpt.plot()
        plot_iv.append_data_point([curs[n], V_iv, dVdI_cos, dVdI_sin, dVdI, dVdI_namp])
        plot_iv.plot()

        # Record current in case the acquisition gets stopped		
        curtemp = curs[n]          
        
        # let the user interface update        
        window.process_events()
    
    # Return current to 0
    current_ramp(curtemp,0)  
    
    # Fit resistance distribution to a parabola
    popt, perr = fit_dVdI(curs, plot_iv[4])
    eta = popt[1]
    
    # Update resistance label
    l_imp_iv.set_text("Offset Voltage = " + str(_n.round(100000*offset)/100) + "mV | Mean Resistance = " + str(_n.round(100*R_L)/100) + "Ohm, eta = " + "%.2f" % (eta/1e6) + "Ohm/(mA)^2")        
    
    # Insert headers
    plot_iv.insert_header("R_L", R_L)
    plot_iv.insert_header("eta", eta)
    
    # uncheck the fire button
    b_fourpt.set_checked(False) 


def sa_calib_step_and_plot():
    """
    Takes the spectrum for a single frequency step for the two-source calibration
    in the spectrum analyzer experiment
    """      
    # Reinitialize databoxes and add offset mixdown voltage in header
    reinit_sweep_databoxes()   
   
    # Set up analog output and input tasks (and activate pulse modulation)
    ao = ao_task_pulm()
    ai = ai_task_sweep()

    # Loop over frequencies
    for n in range(0,cfg_pcalib["SA_Calib/Steps"]):     
        
        # Jump to next list point 
        b.list_change_index(n)
        
        # Settling        
        time.sleep(cfg_hardware['SettlingTime'])
        
        # Start analog input
        ai.start()
        
        # Start analog output            
        ao.start()      
        
        # Wait until acquisition ends
        ao.wait_and_clean() 
                
        # retrieve the data
        y = ai.read_and_clean() 
    
        # store the data
        plot_raw["V0"] = y[0]
        #plot_raw["PULM"] = y[1]     
        
        # Plot raw data
        plot_raw.plot()
          
        # Store the data in the databox, and plot sweep  
        V_mean = _n.mean(y[0])
        V_rms  = _n.sqrt(_n.mean(y[0]**2))
        P_mean = 10**(_n.interp(_n.mean(y[0]), cd_calib[:,1], cd_calib[:,0]) / 10)  # power in milliwatts
        plot_sweep.append_data_point([freqs_calib[n], V_mean, V_rms, P_mean])
        if n == 0: 
            plot_sweep.insert_header("V_off", V_off)     
        plot_sweep.plot()

        # let the user interface update
        window.process_events()
        

def get_RL(dc_cur = 0.0):
    """
    Measures the resistance of the device in a single step.
    """
    global ao_fourpt, cos_iv, sin_iv, R_L, I_ac

    RL_changed()    
    
    ao_fourpt, cos_iv, sin_iv = init_fourpt(dc_cur, cfg_fourpt["ModAmp"]) 
    
    # Prepare analog output and input tasks
    ao = ao_task_fourpt()
    ai = ai_task_fourpt()

    # Start analog input
    ai.start()
    
    # Start analog output
    ao.start()

    # Wait until acquisition ends
    ao.wait_and_clean()
    
    # retrieve the data
    y = ai.read_and_clean()    

    # Get voltage offset
    offset = _n.mean(y[0])

    # Get mean resistance
    V_iv, dVdI_cos, dVdI_sin, dVdI, dVdI_namp = fourpt_analysis(y, 0, offset)
    
    # Calculate resistance (from unamplified channel)
    R_L = (dVdI_namp / (cfg_hardware['DCSourceR'] - dVdI_namp)) * cfg_hardware['DCSourceR']
    get_Iac()

    #print dVdI, R_L, I_ac   
    
    l_rliac_fmr.set_text("Device resistance = " + str(round(R_L*10)/10) + "Ohm, AC current = " + str(round(I_ac*100000)/100) + "mA")
    
    window.process_events()
    
    return R_L, offset


def fire_quicksweep(enabled):
    """
    Called whenever the "Quick sweep" button is pressed.
    This just takes infinite sweeps without background subtraction nor resistance readout.
    """    
    if tabs_plots.get_current_tab() not in [2,3]: 
        tabs_plots.set_current_tab(4)

    # don't do anything if we're disabling the acquisition (the loop will
    # finish and end gracefully)
    if not enabled: return

    # Insert headers
    plot_sweep.insert_header('I_dc', cfg_sweep['DC_Current']) 
    plot_sweep.insert_header('PULM_Period', cfg_sweep['PULM/Period'])
    plot_sweep.insert_header('PULM_Repeats', cfg_sweep['PULM/Repeats'])
    plot_sweep.insert_header('R_L', R_L)
    plot_sweep.insert_header('eta', eta)
        
    # Initialize acquisition and analysis arrays
    ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays()

    # Do the slow current ramp
    current_ramp(0,cfg_sweep["DC_Current"])    
    time.sleep(3)    
    
    # Turn on Anapico and set up for manual list mode
    a.list_man_setup(freqs,pows)
    time.sleep(1)   

    # Set loop to be infinite
    itertemp = cfg_sweep["Iterations"]
    cfg_sweep["Iterations"] = 0

    # start the loop. This will wait until the button is deactivated  
    n_fire_count.set_value(0)
    while b_quicksweep.is_checked():

        # see if we should stop
        if n_fire_count.get_value() >= cfg_sweep['Iterations'] \
        and not cfg_sweep['Iterations'] == 0: break

        # acquire and plot data
        sweep_and_plot()

        # increment the counters
        n_fire_count.increment(1)

        # update the window
        window.process_events()

    # Do the slow current ramp
    current_ramp(cfg_sweep["DC_Current"],0)    
    
    # Turn output off
    a.quit_list_mode()
    
    # Reset the previously specified number of acquisitions    
    cfg_sweep["Iterations"] = itertemp
    b_quicksweep.set_checked(False)
        
# connect this function
window.connect(b_quicksweep.signal_clicked, fire_quicksweep)
        

def fmr_and_plot(enabled, init_cur = 0.0, final_cur = 0.0, cur_sweep = False, getRL = True):
    """
    Routine recording an FMR plot safely (with current ramps)
    This takes a full sweep, including background subtraction and resistance readout.
    """
    
    if tabs_plots.get_current_tab() not in [2,3]: 
        tabs_plots.set_current_tab(4)

    # don't do anything if we're disabling the acquisition (the loop will
    # finish and end gracefully)
    if not enabled: return

    # Get the resistance of the device and the AC current (if available from previous measurements)
    # First do the resistance measurement on a single point   

    if getRL and b_resis.get_value():
        global R_L    
        time.sleep(3) 
        R_L, V_off = get_RL(init_cur)
              
    # Insert headers
    plot_fmr.insert_header('Experiment', cfg_sweep['Experiment'])
    plot_fmr.insert_header('SigGen', cfg_hardware['RFSource'])
    plot_fmr.insert_header('AmpGain', cfg_hardware['Amp/Gain'])
    plot_fmr.insert_header('I_dc', final_cur) 
    plot_fmr.insert_header('PULM_Period', cfg_sweep['PULM/Period'])
    plot_fmr.insert_header('PULM_Repeats', cfg_sweep['PULM/Repeats'])
    plot_fmr.insert_header('R_L', R_L)
    plot_fmr.insert_header('eta', eta)
          
    # Do the slow current ramp
    current_ramp(init_cur, final_cur) 

    # Initialize acquisition and analysis arrays
    ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(dc_cur = final_cur)    
    
    # Turn on Anapico and set up for manual list mode
    a.list_man_setup(freqs,pows)
    time.sleep(3) 
    
    ### FIRST STEP - REFERENCE SPECTRUM    
    
    # Reinitialize reference spectrum and FMR databox
    global fmr_ref
    fmr_ref = [0.0]*cfg_sweep['Steps']
    fmr_stdref = _n.array([0.0]*cfg_sweep['Steps'])

    reinit_fmr_databoxes()  
        
    def fmr_plot():
        if cfg_sweep["Experiment"] == "Driven Dynamics": 
            fmr_all = fmr_dd_plot()
        elif cfg_sweep["Experiment"] == "Spectrum Analyzer": 
            fmr_all = fmr_sa_plot()
        return fmr_all

    if b_sweepref.is_checked(): 
        # Prompt user to move magnet to n*pi/2 radians
        w_magnet_1 = egg.gui.Window("MAGNET!")
        w_magnet_1.place_object(egg.gui.GridLayout(False))
        w_magnet_1.place_object(egg.gui.Label("Please move the magnet to the reference position."))
        w_magnet_1.new_autorow()
        b_magnet_1 = w_magnet_1.place_object(egg.gui.Button("OKAY."))
        b_magnet_1.set_checkable(True); b_magnet_1.set_checked(False)
        w_magnet_1.show()
        while not b_magnet_1.is_checked():   
            time.sleep(0.1); w_magnet_1.process_events()
        w_magnet_1.close()

        # start the loop. This will wait until the button is deactivated  
        n_fire_count.set_value(0)
        while b_sweep.is_checked():
    
            # see if we should stop
            if n_fire_count.get_value() >= cfg_sweep['Iterations'] \
            and not cfg_sweep['Iterations'] == 0: break
    
            # acquire and plot data
            sweep_and_plot()
    
            # Update the FMR plot
            fmr_all = fmr_plot()
    
            # increment the counters
            n_fire_count.increment(1)
    
            # update the window
            window.process_events()
    
        # Get the reference spectrum and its standard deviations
        fmr_ref = plot_fmr['dR']
        fmr_stdref = _n.array(_n.std(fmr_all, axis = 0)) / _n.sqrt(n_fire_count.get_value())


    # SECOND STEP - Actual acquisition

    # Prompt user to move magnet to the acquisition position
    if b_sweepref.is_checked(): 
        w_magnet_2 = egg.gui.Window("MAGNET!")
        w_magnet_2.place_object(egg.gui.GridLayout(False))
        w_magnet_2.place_object(egg.gui.Label("Please move the magnet to the measurement position"))
        w_magnet_2.show()
        w_magnet_2.new_autorow()
        b_magnet_2 = w_magnet_2.place_object(egg.gui.Button("OKAY."))
        b_magnet_2.set_checkable(True); b_magnet_2.set_checked(False)
        w_magnet_2.show()
        while not b_magnet_2.is_checked():
            time.sleep(0.1); w_magnet_2.process_events()
        w_magnet_2.close()

    # start the loop. This will wait until the button is deactivated  
    n_fire_count.set_value(0)
    while b_sweep.is_checked():

        # see if we should stop
        if n_fire_count.get_value() >= cfg_sweep['Iterations'] \
        and not cfg_sweep['Iterations'] == 0: break

        # acquire and plot data
        sweep_and_plot()

        # Update the FMR plot
        fmr_all = fmr_plot()

        # increment the counters
        n_fire_count.increment(1)

        # update the window
        window.process_events()
      
    # Turn output off
    a.quit_list_mode()

    # Return to zero current if we are not sweeping the current
    if not cur_sweep: current_ramp(final_cur,init_cur)
    
    # Compute and store the standard deviations of the sweep
    if cfg_sweep['Iterations'] != 1:
        fmr_std = _n.array(_n.std(fmr_all, axis = 0)) / _n.sqrt(n_fire_count.get_value())
        plot_fmr['std_dR'] = fmr_std 
        fmr_std_c = _n.sqrt(fmr_std**2 + fmr_stdref**2)
        plot_fmr['std_dR_c'] = fmr_std_c
    
    # uncheck the fire button
    b_sweep.set_checked(False)
    
    # Insert headers
    plot_fmr.insert_header('iterations', n_fire_count.get_value())
    
# connect this function
def fire_sweep(enabled):
    """    
    Called whenever the "Sweep" button is pressed.
    """
    fmr_and_plot(enabled, init_cur = 0.0, final_cur = cfg_sweep["DC_Current"])    
    
window.connect(b_sweep.signal_clicked, fire_sweep)



def fire_bsweep(enabled):
    """
    Fires up a magnetic field sweep using the motorized actuators.
    """    
    global conex_path
    
    # Reinitialize databox
    reinit_bsweep_databoxes()
   
    # Loop over magnetic field steps
    for k in range(cfg_bsweep['Steps']):
                
        # Continue loop if bsweep button remains checked
        if b_bsweep_fire.is_checked():

            # Check something that was always checked in the loop (why? -Adrian)
            b_sweep.set_checked(True)     
                             
            # Move the magnet
            conex_move.SetPosition( tuple( conex_path[k,:] ) )
            
            # Sweep!
            fire_sweep(True)
    
            # Store and plot
            plot_bsweep.append_data_row(0, plot_fmr[2])
            plot_bsweep.append_data_row(1, plot_fmr[3])
            plot_bsweep.append_data_row(2, plot_fmr[4])
            plot_bsweep.append_data_row(3, plot_fmr[5])
          
            # Fill in color plot
            plot_bsweep.update_plots()
      
        # Leave loop gracefully if bsweep button unchecked
        else:
            #I don't know if this is everything I need to leave gracefully -Adrian
            break
            b_sweep.set_checked(False)
        
    # uncheck the fire button
    b_bsweep_fire.set_checked(False)

window.connect(b_bsweep.signal_clicked, fire_bsweep)


def fire_isweep(enabled):
    
    global currentvals
    
    plot_isweep.set_size_and_bounds(size = [cfg_sweep['Steps'], cfg_isweep['Steps']], \
    bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_isweep['Start'], cfg_isweep['Stop']])
 
    global R_L    
    R_L, V_off = get_RL()
   
   # Loop over magnetic field steps
    for k in range(cfg_isweep['Steps']):
                
        # Continue loop if bsweep button remains checked
        if b_isweep.is_checked():

            # Check something that was always checked in the loop (why? -Adrian)
            b_sweep.set_checked(True)     
                             
            ### ACTUAL CURRENT SWEEPING           
            # Sweep!
            if k == 0:
                fmr_and_plot(True, init_cur = 0, final_cur = currentvals[k], cur_sweep = True, getRL = False)
            else:
                fmr_and_plot(True, init_cur = currentvals[k-1], final_cur = currentvals[k], cur_sweep = True, getRL = False)
    
            # Store and plot
            plot_isweep.append_data_row(0, plot_fmr[2])
            plot_isweep.append_data_row(1, plot_fmr[3])
            plot_isweep.append_data_row(2, plot_fmr[4])
            plot_isweep.append_data_row(3, plot_fmr[5])
          
            # Fill in color plot
            plot_isweep.update_plots()
            
            index_cur = k

        # Leave loop gracefully if bsweep button unchecked
        else:         
            index_cur = k
            break
            b_sweep.set_checked(False)         
            
        
    # uncheck the fire button
    b_isweep.set_checked(False)

    # Ramp current back to zero    
    current_ramp(currentvals[index_cur], 0)

window.connect(b_isweep.signal_clicked, fire_isweep)          


def fire_IV(enabled):
    """
    Called whenever the "IV Curve" button is pressed.
    """          
    
    if tabs_plots.get_current_tab() != 5: 
        tabs_plots.set_current_tab(6)
    
    # Reinitialize databoxes
    RL_changed()
    
    while b_fourpt.is_checked():
        # acquire and plot data
        iv_and_plot()
        
        # update the window
        window.process_events()
        
    # uncheck the fire button
    b_fourpt.set_checked(False)

# connect the function
window.connect(b_fourpt.signal_clicked, fire_IV)


def fire_calib_dd():
    """
    Routine doing a power calibration by using the joule heating mixdown term.
    """
    b_calib.set_checked(True)
    n_fire_count.set_value(0)
    
    # Reset the calibration and reinitialize the power calibration databox
    ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays()
    b_calibreset_clicked()    
    global calib_factors		
    calib_factors = _n.array([[1.0]*cfg_sweep['Steps']] * 3)

    # Insert headers    
    plot_calib.insert_header('Experiment', cfg_sweep['Experiment'])
    plot_calib.insert_header('SigGen', cfg_hardware['RFSource'])
    plot_calib.insert_header('AmpGain', cfg_hardware['Amp/Gain'])
    plot_calib.insert_header('Mode', cfg_pcalib["DD_Calib/Mode"])
    plot_calib.insert_header('Power', cfg_sweep["Power"])
    plot_calib.insert_header('I_DC', cfg_pcalib["DD_Calib/DC_Current"])
    plot_calib.insert_header('PULM_Period', cfg_sweep["PULM/Period"])
    plot_calib.insert_header('PULM_Repeats', cfg_sweep["PULM/Repeats"])
    plot_calib.insert_header('Iterations', cfg_pcalib["DD_Calib/Iterations"])
    plot_calib.insert_header('dVdI_Current', cfg_pcalib["DD_Calib/dVdI_Current"])
    plot_calib.insert_header('dVdI_Steps', cfg_pcalib["DD_Calib/dVdI_Steps"])
    
    # Ramp the current back down    
    RL_changed()
    current_ramp(0,0)    

    # Prompt user to move magnet to n*pi/2 radians    
    w_magnet_1 = egg.gui.Window("MAGNET!")
    w_magnet_1.place_object(egg.gui.GridLayout(False))
    w_magnet_1.place_object(egg.gui.Label("Please move the magnet to the measurement position."))
    w_magnet_1.new_autorow()
    b_magnet_1 = w_magnet_1.place_object(egg.gui.Button("OKAY."))
    b_magnet_1.set_checkable(True); b_magnet_1.set_checked(False)
    w_magnet_1.show()
    try:
        while not b_magnet_1.is_checked():   
            time.sleep(0.1); 
            w_magnet_1.process_events()
        w_magnet_1.close()
    except RuntimeError:
        b_calib.set_checked(False)
        print "Calibration aborted."
    
    # Break if calibration button is not checked
    if not b_calib.is_checked(): return
 
    ### START CALIBRATION        

    # STEP 0 : Get DC offset by doing a sweep at minimum power and no DC current
    print "0: Get device resistance, heating factor and voltage offset..."

    # Get the resistance and heating factor of the device        
    # Start by storing temporary values for the IV curve
    ivstart_temp = cfg_fourpt["Start"]
    ivstop_temp  = cfg_fourpt["Stop"]
    ivsteps_temp = cfg_fourpt["Steps"]
    
    # Replace them by values specific to the calibration
    cfg_fourpt["Start"] = -cfg_pcalib["DD_Calib/dVdI_Current"]
    cfg_fourpt["Stop"]  =  cfg_pcalib["DD_Calib/dVdI_Current"]
    cfg_fourpt["Steps"] =  cfg_pcalib["DD_Calib/dVdI_Steps"]
    
    # Get the IV curve, as well as R_L and eta in the process
    b_fourpt.set_checked(True)
    tabs_plots.set_current_tab(6)
    iv_and_plot()
    b_fourpt.set_checked(False)
    
    # Store R_L and eta in header
    plot_calib.insert_header('R_L', R_L)
    plot_calib.insert_header('eta', eta)    
    
    # Reinitialize fourpt parameters
    cfg_fourpt["Start"] = ivstart_temp
    cfg_fourpt["Stop"] = ivstop_temp
    cfg_fourpt["Steps"] = ivsteps_temp
        
    # Now get the DC offset voltage
    global V_off

    # Set the values for the RF and DC currents for the 1st step
    curtemp = cfg_sweep['DC_Current']
    powtemp = cfg_sweep['Power']    
    powmtemp= cfg_sweep['PowerMode']
    stepstemp=cfg_sweep['Steps']     

    cfg_sweep['Power'] = -30 
    cfg_sweep['Steps'] = 1        

    # Initialize acquisition and analysis arrays
    ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(dc_cur = 0.0)
    
    # Turn on signal generator and get the offset
    a.list_man_setup(freqs,pows)
    time.sleep(1)  
    sweep_and_plot()
    global V_off
    V_off = _n.mean(plot_sweep['V_sin'])       
                
    # Reinitialize values for the second step of the power calibration    
    cfg_sweep['Power'] = powtemp
    cfg_sweep['PowerMode'] = 'Constant'
    cfg_sweep['Steps'] = stepstemp
    cfg_sweep_changed()
    
    # Iterate over the next steps as specified
    for k in range(cfg_pcalib["DD_Calib/Iterations"]): 
     
        # STEP 1 : Do a first calibration by working at constant power    
        print "1: Calibrate at constant output power..."
 
        # Do the slow current ramp
        current_ramp(0,cfg_pcalib["DD_Calib/DC_Current"]) 
        time.sleep(3)     
       
        # Break if calibration button is not checked
        if not b_calib.is_checked(): return
        
        # Switch tab
        if tabs_plots.get_current_tab() != 2: tabs_plots.set_current_tab(3)
        
    
        # Initialize acquisition and analysis arrays
        ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(dc_cur = cfg_pcalib["DD_Calib/DC_Current"])
        
        # Set up Anapico for manual list mode with updated frequencies
        a.list_man_setup(freqs,pows)
        time.sleep(1)
        
        # acquire and plot data
        sweep_and_plot() 
        
        # Store as the positive current calibration signal
        V_calib_pos = plot_sweep['V_sin']        
        V_calib_neg = 0.0   
       
        
        # If we do background subtraction, do the same procedure for the opposite current direction
        # First define a variable that tells us if we'll be at +ve or -ve current after the full step 1                      
        invert_cur  = 1    
        
        if cfg_pcalib["DD_Calib/Mode"] == "Background Subtraction":  
            invert_cur = -1            
            
            # Do the slow current ramp
            current_ramp(cfg_pcalib["DD_Calib/DC_Current"], -cfg_pcalib["DD_Calib/DC_Current"]) 
            time.sleep(3)

            # Break if calibration button is not checked
            if not b_calib.is_checked(): return      

            # Switch tab
            if tabs_plots.get_current_tab() != 2: tabs_plots.set_current_tab(3)

            # Initialize acquisition and analysis arrays
            ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(-cfg_pcalib["DD_Calib/DC_Current"])
            
            # Set up Anapico for manual list mode with updated frequencies
            a.list_man_setup(freqs,pows)
        
            # acquire and plot data
            sweep_and_plot()
            
            # Store as the negative current calibration signal
            V_calib_neg = plot_sweep['V_sin']
       
        # Get the calibration data
        global f_calib, p_calib, p_factor, p_factor_mixed
        V_calib = (V_calib_pos - V_calib_neg)/2
        f_calib, p_calib, p_factor = power_calibration(freqs,V_calib)
        p_factor_mixed = [1]*len(p_factor)  
        
        # Get the output current (calculated for an output power of 0dBm)
        plot_calib['I_f'] = _n.sqrt( 2/(3*eta*cfg_pcalib["DD_Calib/DC_Current"]) * V_calib / cfg_hardware["Amp/Gain"]) * 10**(-cfg_sweep["Power"]/20) 
        
        # Break if calibration button is not checked
        if not b_calib.is_checked(): return
  
  
        # Move on to the fine calibration step after adjusting the output powers. 
        # Stop here if we work in constant output power mode.
        if powmtemp != 'Constant': 
            
            # STEP 2: Recalibrate to account for non-linearities in the output power.
            print "2: Recalibrate with adjusted output power levels..."              
                                  
            # Break if calibration button is not checked
            if not b_calib.is_checked(): return     
            
            # Set up power mode and initialize acquisition and analysis arrays      
            cfg_sweep['PowerMode'] = 'Calibrated'
            ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(invert_cur*cfg_pcalib["DD_Calib/DC_Current"])
            a.list_man_setup(freqs,pows)
    
            # Run.             
            sweep_and_plot()
            V_calib_pos = plot_sweep['V_sin']
            V_calib_neg = 0

            # If we do background subtraction, do the same procedure for the opposite current direction
            if cfg_pcalib["DD_Calib/Mode"] == "Background Subtraction":   
                    
                # Break if calibration button is not checked
                if not b_calib.is_checked(): break
                       
                # Do the slow current ramp
                current_ramp(-cfg_pcalib["DD_Calib/DC_Current"], cfg_pcalib["DD_Calib/DC_Current"]) 
                time.sleep(3) 

                # Set up power mode and initialize acquisition and analysis arrays      
                cfg_sweep['PowerMode'] = 'Calibrated'
                ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(cfg_pcalib["DD_Calib/DC_Current"])
                a.list_man_setup(freqs,pows)
        
                # Run.             
                sweep_and_plot()
                V_calib_neg = plot_sweep['V_sin']                   
                   
            # Get the resulting power factor
            V_calib = invert_cur*(V_calib_pos - V_calib_neg)/2
            f_calib, p_c2, p_factor_mixed = power_calibration(freqs,V_calib)

            # Recalculate output current (obtained for an output power of 0dBm)
            #V_c = V_calib[0]*_n.ones(cfg_sweep["Steps"])
            #plot_calib['I_f'] = _n.sqrt( 2/(3*eta*cfg_pcalib["DD_Calib/DC_Current"]) * V_c / cfg_hardware["Amp/Gain"]) * 10**(-cfg_sweep["Power"]/20.0)         
        
        # update the window
        window.process_events()        
        
        # POWER CALIBRATION FINISHED. 
        
        # don't do anything if we're disabling the acquisition (the loop will
        # finish and end gracefully)
        if not b_calib.is_checked(): return
            
        # Plot the results.
        tabs_plots.set_current_tab(1)
        calib_dd_plot()
    
        # update the window and increment counter
        window.process_events()            
        n_fire_count.increment(1)        
        
        # Ramp the current back down    
        current_ramp(cfg_pcalib["DD_Calib/DC_Current"],0) 

    # Reset RF and DC currents values
    cfg_sweep['DC_Current'] = curtemp
    cfg_sweep['Power'] = powtemp
    cfg_sweep['PowerMode'] = powmtemp
    
    
    # FULL CALIBRATION FINISHED.    
    # uncheck the fire button
    b_calib.set_checked(False)

    # Turn Anapico output off    
    a.quit_list_mode()
    
    # Reset the counters
    n_fire_count.set_value(0)    
   

#def fire_calib_dd_bckg(): 
#    """
#    Routine doing a power calibration by doing a subtraction of the resonant features
#    by taking calibration curves at +x mA and -x mA.
#    """    
#    b_calib.set_checked(True)
#    n_fire_count.set_value(0)
#    
#    # Reset the calibration and reinitialize the power calibration databox
#    ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays()
#    b_calibreset_clicked()
#    plot_calib['f']     = freqs
#    plot_calib['P']     = [cfg_sweep['Power']]*cfg_sweep['Steps']
#    plot_calib['P_f']   = [1]*cfg_sweep['Steps']
#    plot_calib['P_f_m'] = [1]*cfg_sweep['Steps']    
#    global calib_factors		
#    calib_factors = _n.array([[1.0]*cfg_sweep['Steps']] * 3)
#    
#    # Ramp the current back down    
#    RL_changed()
#    current_ramp(0,0) 
#
#    # Set the values for the RF and DC currents for the 1st step
#    curtemp = cfg_sweep['DC_Current']
#    powtemp = cfg_sweep['Power']    
#    powmtemp= cfg_sweep['PowerMode']
#    stepstemp=cfg_sweep['Steps']        
#
#    # Prompt user to move magnet to n*pi/2 radians    
#    w_magnet_1 = egg.gui.Window("MAGNET!")
#    w_magnet_1.place_object(egg.gui.GridLayout(False))
#    w_magnet_1.place_object(egg.gui.Label("Please move the magnet to the measurement position."))
#    w_magnet_1.new_autorow()
#    b_magnet_1 = w_magnet_1.place_object(egg.gui.Button("OKAY."))
#    b_magnet_1.set_checkable(True); b_magnet_1.set_checked(False)
#    w_magnet_1.show()
#    try:
#        while not b_magnet_1.is_checked():   
#            time.sleep(0.1); 
#            w_magnet_1.process_events()
#        w_magnet_1.close()
#    except RuntimeError:
#        b_calib.set_checked(False)
#        print "Calibration aborted."
#    
#    # Break if calibration button is not checked
#    if not b_calib.is_checked(): return
# 
#    ### START CALIBRATION        
#    
#    for k in range(cfg_pcalib["DD_Calib/Iterations"]):
#
#        # STEP 0 : Get DC offset by doing a sweep at minimum power and no DC current
#        print "Step 0..."
#    
#        global V_off
#        R_L, offset = get_RL()
#        #print "V_off = " + str(round(V_off,3))  
#        
#        cfg_sweep['Power'] = -30 
#        cfg_sweep['Steps'] = 1        
#        
#        # Initialize acquisition and analysis arrays
#        ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(dc_cur = 0.0)
#        
#        # Turn on Anapico and set up for manual list mode
#        a.list_man_setup(freqs,pows)
#        time.sleep(1)  
#        sweep_and_plot()
#        global V_off
#        V_off = _n.mean(plot_sweep['V_sin'])   
#        print "V_off = " + str(round(V_off,3))       
##        
##        
#        # Reinitialize values for the second step of the power calibration    
#        cfg_sweep['Power'] = powtemp
#        cfg_sweep['PowerMode'] = 'Constant'
#        cfg_sweep['Steps'] = stepstemp
#        cfg_sweep_changed()
#        
#        # Do the slow current ramp
#        current_ramp(0,cfg_pcalib["DD_Calib/DC_Current"]) 
#        time.sleep(3)        
#        
#        # Reset calibration
#        b_calibreset_clicked()
#        
#        # STEP 1a : Do a first calibration by working at constant power    
#        print "Step 1a..."
#        
#        # Break if calibration button is not checked
#        if not b_calib.is_checked(): return
#        
#        # Switch tab
#        if tabs_plots.get_current_tab() != 2: 
#            tabs_plots.set_current_tab(3)
#        
#    
#        # Initialize acquisition and analysis arrays
#        ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(dc_cur = cfg_pcalib["DD_Calib/DC_Current"])
#        
#        # Set up Anapico for manual list mode with updated frequencies
#        a.list_man_setup(freqs,pows)
#        time.sleep(1)
#        
#        # acquire and plot data
#        sweep_and_plot()
#        
#        # Store as the positive current calibration signal
#        V_calib_pos = plot_sweep['V_sin']
#    
#    #        # Get the calibration data
#    #        global f_calib, p_calib, p_factor, p_factor_mixed
#    #        f_calib, p_calib, p_factor = power_calibration(freqs,plot_sweep['V_sin'])
#    #        p_factor_mixed = [1]*len(p_factor)
#
#     
#        # Do the slow current ramp
#        current_ramp(cfg_pcalib["DD_Calib/DC_Current"], -cfg_pcalib["DD_Calib/DC_Current"]) 
#        time.sleep(3)
#    
#        # STEP 1b : Same as step 1a, with negative DC current
#        print "Step 1b..."
#
#        # Break if calibration button is not checked
#        if not b_calib.is_checked(): return        
#        
#        # Switch tab
#        if tabs_plots.get_current_tab() != 2: 
#            tabs_plots.set_current_tab(3)
#        
#    
#        # Initialize acquisition and analysis arrays
#        ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(-cfg_pcalib["DD_Calib/DC_Current"])
#        
#        # Set up Anapico for manual list mode with updated frequencies
#        a.list_man_setup(freqs,pows)
#    
#        # acquire and plot data
#        sweep_and_plot()
#        
#        # Store as the negative current calibration signal
#        V_calib_neg = plot_sweep['V_sin']
#    
#        # Get the calibration data
#        global f_calib, p_calib, p_factor, p_factor_mixed
#        V_calib = (V_calib_pos - V_calib_neg)/2
#        f_calib, p_calib, p_factor = power_calibration(freqs,V_calib)
#        p_factor_mixed = [1]*len(p_factor)  
#        
#        # Break if calibration button is not checked
#        if not b_calib.is_checked(): return
#    
#        # Stop if we work in constant mode
#        if powmtemp != 'Constant':
#            
#            # STEP 2a : Same thing with corrected powers (for pre-calibrated / mixed modes) to flatten out generator nonlinearities
#            print "Step 2a..."
#            
#            # Break if calibration button is not checked
#            if not b_calib.is_checked(): return        
#            
#            # Set up power mode and initialize acquisition and analysis arrays      
#            if b_calibthird.is_checked(): 
#                cfg_sweep['PowerMode'] = 'Calibrated'
#                ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(-cfg_pcalib["DD_Calib/DC_Current"])
#                a.list_man_setup(freqs,pows)
#        
#                # Run.    
#                sweep_and_plot()
#                V_calib_neg = plot_sweep['V_sin']
#    #                else:
#    #                    V_calib_neg = (N*V_calib_neg + plot_sweep['V_sin'])/(N+1)
#            
#            # Break if calibration button is not checked
#            if not b_calib.is_checked(): return
#        
#            
#            # STEP 2b : Same as step 2a, with positive current.
#            print "Step 2b..."
#            
#            # Break if calibration button is not checked
#            if not b_calib.is_checked(): return        
#            
#            # Do the slow current ramp
#            current_ramp(-cfg_pcalib["DD_Calib/DC_Current"], cfg_pcalib["DD_Calib/DC_Current"]) 
#            time.sleep(3)        
#            
#            # Set up power mode and initialize acquisition and analysis arrays      
#            if b_calibthird.is_checked(): 
#                cfg_sweep['PowerMode'] = 'Calibrated'
#                ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays(cfg_pcalib["DD_Calib/DC_Current"])
#                a.list_man_setup(freqs,pows)
#        
#                # Run.             
#                sweep_and_plot()
#                V_calib_pos = plot_sweep['V_sin']
#    #            else:
#    #                V_calib_pos = (N*V_calib_pos + plot_sweep['V_sin'])/(N+1)
#                
#                # Break if calibration button is not checked
#                if not b_calib.is_checked(): break
#                           
#                # Get the resulting power factor
#                V_calib = (V_calib_pos - V_calib_neg)/2
#                f_calib, p_c2, p_factor_mixed = power_calibration(freqs,V_calib)
#        #            p_factor = _n.multiply(p_factor, p_factor_mixed)
#        #            p_calib  = -10*_n.log10(p_factor)       
#        
#        # update the window
#        window.process_events()        
#        
#        # POWER CALIBRATION FINISHED. 
#        
#        # don't do anything if we're disabling the acquisition (the loop will
#        # finish and end gracefully)
#        if not b_calib.is_checked(): return
#            
#        # Plot the results.
#        tabs_plots.set_current_tab(1)
#        calib_dd_plot()
#    
#        # update the window and increment counter
#        window.process_events()            
#        n_fire_count.increment(1)        
#        
#        # Ramp the current back down    
#        current_ramp(cfg_pcalib["DD_Calib/DC_Current"],0) 
#
#    # Reset RF and DC currents values
#    cfg_sweep['DC_Current'] = curtemp
#    cfg_sweep['Power'] = powtemp
#    cfg_sweep['PowerMode'] = powmtemp
#    
#    
#    # FULL CALIBRATION FINISHED.    
#    # uncheck the fire button
#    b_calib.set_checked(False)
#
#    # Turn Anapico output off    
#    a.quit_list_mode()
#    
#    # Reset the counters
#    n_fire_count.set_value(0)    


def fire_calib_sa():
    b_calib.set_checked(True)
    # Ask for save path
    w_calib_sa = egg.gui.Window("Calibration Save Path")
    w_calib_sa.place_object(egg.gui.GridLayout(False))
    w_calib_sa.place_object(egg.gui.Label("Please specify the save path for the calibration files."))
    w_calib_sa.new_autorow()
    t_calib_sa = w_calib_sa.place_object(egg.gui.TextBox("M:/Data - Magnetometry/Probe Station Data/"))
    b_calib_sa = w_calib_sa.place_object(egg.gui.Button("OKAY."))
    b_calib_sa.set_checkable(True); b_calib_sa.set_checked(False)
    w_calib_sa.show()
    while not b_calib_sa.is_checked():   
        time.sleep(0.1); w_calib_sa.process_events()
    save_folder = t_calib_sa.get_text()
    try:    
        os.mkdir(save_folder)
    except WindowsError:
        save_folder += "_1"
        os.mkdir(save_folder)
    w_calib_sa.close()

    # Set up databox plot
    plot_calib_sa.set_size_and_bounds(size = [cfg_sweep['Steps'], cfg_pcalib['SA_Calib/PowerSteps']], \
    bnds = [cfg_sweep['Start'], cfg_sweep['Stop'], cfg_pcalib['SA_Calib/PowerStart'], cfg_pcalib['SA_Calib/PowerStop']])

    
    # Reset the calibration and reinitialize the power calibration databox
    ao_t, ao_sw, ao_dc, cos_li, sin_li, freqs, pows, pows_corr = init_sweep_arrays()
    b_calibreset_clicked()
    reinit_calib_databoxes()
 
    # Set up signal generators for manual list mode with updated frequencies
    a.list_man_setup(freqs,pows)
    powers = _n.linspace(cfg_pcalib["SA_Calib/PowerStart"], cfg_pcalib["SA_Calib/PowerStop"], cfg_pcalib["SA_Calib/PowerSteps"])    
    b_sweep.set_checked(True)    
    
    # Loop over main signal generator frequencies
    for k in range(len(powers)):
        # Break if calibration button is not checked
        if not b_sweep.is_checked(): break

        # Reinitialize data row
        p_out  = _n.array([0.0]*cfg_sweep['Steps'])
        v_base = _n.array([0.0]*cfg_sweep['Steps'])
        gain   = _n.array([1.0]*cfg_sweep['Steps'])
            
        for n in range(0, cfg_sweep['Steps']):
            
            # Break if calibration button is not checked
            if not b_sweep.is_checked(): break
                
            # Take some sweeps    
            print n
            # Change reference frequency
            a.list_change_index(n)
            init_sa_calib_arrays(freqs[n], powers[k])
            
            # Initialize frequency and power arrays on calibration signal generator, plot and save
            b.list_man_setup(freqs_calib, pows_calib)
            sa_calib_step_and_plot()
            plot_sweep.save_file(save_folder+"/"+str(powers[k])+'dBm_'+str(n)+".swp")  
            b.quit_list_mode()
            
            # Calculate the output power, store in the databox and plot
            p_out[n], v_base[n] =  sa_calib_integration(plot_sweep['P_out'], freqs[n]) 
            gain[n]  =  p_out[n] / 10**(powers[k]/10)
            if _n.isnan(gain[n]): gain[n] = gain[n-1] # Make it safe
            if n == 0:
                plot_calib_sa.append_data_row(0, p_out)
                plot_calib_sa.append_data_row(1, v_base)
                plot_calib_sa.append_data_row(2, gain)
                plot_calib_sa.append_data_row(3, 10.0*_n.log10(gain))
            else:
                plot_calib_sa.update_data_row(0, p_out, k)
                plot_calib_sa.update_data_row(1, v_base, k)
                plot_calib_sa.update_data_row(2, gain, k)
                plot_calib_sa.update_data_row(3, 10.0*_n.log10(gain), k)
            
            plot_calib.append_data_point( [powers[k], freqs[n], p_out, gain, v_base] )
            plot_calib.plot()
            
            plot_calib_sa.update_plots()

    # FULL CALIBRATION FINISHED.    
    # uncheck the fire button
    b_calib.set_checked(False)

    # uncheck the fire button
    b_sweep.set_checked(False)
      
    # Turn output off
    a.quit_list_mode()
    b.quit_list_mode()

    # Save the calibration data and store the interpolation function
    plot_calib.save_file(save_folder+"/calib_data.calib")
    
    # Store the calibration in the local variables
    global sa_p_out, sa_v_base, sa_gain
    sa_p_out  = plot_calib_sa.fetch_data(db = 0)[0]
    sa_v_base = plot_calib_sa.fetch_data(db = 1)[0]
    sa_gain   = plot_calib_sa.fetch_data(db = 2)[0]
    b_calibload_clicked()

    # Reset the counters
    n_fire_count.set_value(0)   


def fire_calib():
    if cfg_sweep["Experiment"] == "Driven Dynamics":
        fire_calib_dd()
#        if cfg_pcalib["DD_Calib/Mode"] == "High Damping":
#            fire_calib_dd_damp()
#        elif cfg_pcalib["DD_Calib/Mode"] == "Background Subtraction":
#            fire_calib_dd_bckg()
    elif cfg_sweep["Experiment"] == "Spectrum Analyzer":
        fire_calib_sa()
    
# connect the function
window.connect(b_calib.signal_clicked, fire_calib)

##########################
### DATA PLOTTING ROUTINES
##########################   

def fmr_dd_plot():
    """
    Plots the FMR spectrum for the driven dynamics experiment (with averaging)
    """
    # Get the sweep and FMR databoxes
    d_sw  = plot_sweep
    d_fmr = plot_fmr
    
    # Calculate the FMR signal from the last sweep (including background subtraction)
    global fmr_ref  
    V_nonres = plot_sweep.headers['V_nonres']
    fmr_sw   = 2/I_ac * (d_sw['V_sinc']-V_nonres) / cfg_hardware['Amp/Gain'] 
    
    # Do the averaging
    N = n_fire_count.get_value()
    d_fmr['dR'] = ( N*d_fmr['dR'] + fmr_sw ) / (N+1)
    d_fmr['dR_c'] = d_fmr['dR_c'] - fmr_ref
    
    # Store sweep in the databox (without offset correction) and plot
    global fmr_all    
    if N == 0:    
        fmr_all = 2/I_ac * d_sw['V_sinc'] / cfg_hardware['Amp/Gain'] 
    else:
        fmr_all = _n.row_stack([fmr_all, 2/I_ac * d_sw['V_sinc'] / cfg_hardware['Amp/Gain']])
    plot_fmr.plot()
    
    # let the user interface update        
    window.process_events()    
    
    return fmr_all


def fmr_sa_plot():
    """
    Plots the FMR spectrum for the spectrum analysis experiment (with averaging)
    """
    # Get the sweep and FMR databoxes
    d_sw  = plot_sweep
    #d_ca  = plot_calib
    d_fmr = plot_fmr
    
    # Substract baseline level
    # d_sw['P_out'] += -_n.interp(freqs, d_ca['f'], d_ca['P_base'])
    
    # Calculate the FMR signal from the last sweep (including background subtraction)    
    fmr_sw = _n.array([])    
    for k in range(0,cfg_sweep['Steps']):
        # Calculate the gain associated to the each data point
        #print f_calib_sa(freqs[k],d_sw['P_out'][k])
        #fmr_sw = _n.append( fmr_sw, f_calib_sa(freqs[k],d_sw['P_out'][k]) )   
        
        # Calculate the input power from the gain look-up table obtained during the calibration        
        p_in = d_sw['P_out'][k] / _n.interp(d_sw['P_out'][k], sa_p_out[:,k], sa_gain[:,k])
#        print p_in, d_sw['P_out'][k]
#        print _n.interp(d_sw['P_out'][k], sa_p_out[:,k], sa_gain[:,k]), sa_p_out[:,k], sa_gain[:,k]
        fmr_sw = _n.append( fmr_sw,  p_in )
    
    # Do the averaging
    N = n_fire_count.get_value()
    d_fmr[1] = ( N*d_fmr[1] + d_sw['P_out']) / (N+1)
    d_fmr[2] = ( N*d_fmr[2] + fmr_sw ) / (N+1)
    d_fmr[3] = d_fmr[2] - fmr_ref
    
    # Store sweep in the databox (without offset correction) and plot
    global fmr_all    
    if N == 0:    
        fmr_all = fmr_sw
    else:
        fmr_all = _n.row_stack([fmr_all, fmr_sw])
    plot_fmr.plot()
    
    # let the user interface update        
    window.process_events()    
    
    return fmr_all


def calib_dd_plot():
    """
    Plots the power calibration curve (with averaging) for the driven dynamics experiment
    """
    global calib_factors, p_factor, p_factor_mixed

#    print p_factor
#    print p_factor_mixed    
#    
    # Do the averaging
    N = n_fire_count.get_value()
    if N == 0:
        calib_factors[1] = p_factor
        calib_factors[2] = p_factor_mixed
    else:
        calib_factors[1] = (N*calib_factors[1] + p_factor)/(N+1)
        calib_factors[2] = (N*calib_factors[2] + p_factor_mixed)/(N+1)
    calib_factors[0] = -10*_n.log10(calib_factors[1])
    
    # Store the databox and plot
    plot_calib[1] = calib_factors[0]
    plot_calib[2] = calib_factors[1]
    plot_calib[3] = calib_factors[2]
    
    plot_calib.plot()    
#
#    print p_factor
#    print p_factor_mixed    

    # let the user interface update        
    window.process_events()          
    

def calib_sa_plot():
    """
    Plots the reference spectrum (with averaging) for the spectrum analyzer experiment
    """
    global plot_sweep, plot_calib
    # Get the sweep and FMR databoxes
    d_sw = plot_sweep
    d_ca = plot_calib    
    
    # Do the averaging
    N = n_fire_count.get_value()
    if N == 0:
        d_ca['V'] = d_sw['V']
        d_ca['V_rms'] = d_sw['V_rms']
    else:
        d_ca['V'] = (N*d_ca['V'] + d_sw['V'])/(N+1)
        d_ca['V_rms'] = (N*d_ca['V_rms'] + d_sw['V_rms'])/(N+1)
    
    # Store the databox and plot
    plot_calib = d_ca
    plot_calib.plot()    

    # let the user interface update        
    window.process_events()      


#################
### DATA ANALYSIS
#################     
    
def boxcar_analysis(y):
    """
    Performs the lock-in (boxcar) analysis on the acquired data
    """      

    y_ar = _n.array(y[0])
        
    # Calculate quadratures of the lock-in measurement (additional pi factor = 1st fourier component of square wave)     
    V_cos = _n.sum(_n.multiply(y_ar,cos_li)) / (cfg_sweep['PULM/Period']*cfg_hardware['AO/Rate']) * _n.pi / cfg_sweep["PULM/Repeats"]
    V_sin = _n.sum(_n.multiply(y_ar,sin_li)) / (cfg_sweep['PULM/Period']*cfg_hardware['AO/Rate']) * _n.pi / cfg_sweep["PULM/Repeats"]
        
    # Get amplitude and phase of the measurement
    # V_mag   = _n.sqrt(V_cos**2 + V_sin**2)
    # phi_li = _n.angle(V_sin+ 1j*V_cos)

    return V_sin, V_cos


def fourpt_analysis(y, cur, offset=0):
    """
    Performs the lock-in (boxcar) analysis on the acquired data
    namp = non-amplified
    """ 
    
    # Separate DC and AC parts    
    V_dc = y[0]
    V_ac = y[1]
    
    # Truncate the dataset to remove the transient from the analysis
    trunc_pt = cfg_hardware['AI/Rate'] * cfg_fourpt["PeriodCutoff"] * cfg_fourpt["ModPeriod"]
    V_dc = y[0][trunc_pt::]
    V_ac = y[1][trunc_pt::]
    cos_iv_trunc = cos_iv[trunc_pt::]
    sin_iv_trunc = sin_iv[trunc_pt::]
    
    
    # Calculate voltage and local resistance 
    V = _n.mean(V_dc[0:(len(ao_fourpt)-1)]) - offset      
     
    # Lock-in analysis for dV/dI 
    integration_index = cfg_fourpt["ModPeriod"]*cfg_hardware['AO/Rate']
    num_periods = (cfg_fourpt["ModRepeats"] - cfg_fourpt["PeriodCutoff"])
    dVdI_namp_cos = _n.sum(_n.multiply(V_dc,cos_iv_trunc)) / integration_index / (cfg_fourpt["ModAmp"]/2) / num_periods
    dVdI_namp_sin = _n.sum(_n.multiply(V_dc,sin_iv_trunc)) / integration_index / (cfg_fourpt["ModAmp"]/2) / num_periods
    dVdI_cos  = _n.sum(_n.multiply(V_ac,cos_iv_trunc))     / integration_index / (cfg_fourpt["ModAmp"]/2) / num_periods / cfg_hardware["Amp/Gain"]  * 10**(-cfg_hardware["Amp/Attenuation"]/10)
    dVdI_sin  = _n.sum(_n.multiply(V_ac,sin_iv_trunc))     / integration_index / (cfg_fourpt["ModAmp"]/2) / num_periods / cfg_hardware["Amp/Gain"]  * 10**(-cfg_hardware["Amp/Attenuation"]/10)
    dVdI_namp = _n.sqrt(dVdI_namp_cos**2 + dVdI_namp_sin**2)
    dVdI  = _n.sqrt(dVdI_cos**2 + dVdI_sin**2)
    
    return V, dVdI_cos, dVdI_sin, dVdI, dVdI_namp
 

def sa_calib_integration(y, f_lo, pad = 0.2):
    """
    Performs the integration over the calibration "step functions"
    """
    # Get the index numbers for all points that we will use to get the baseline level of each plot    
    f_c = cfg_hardware['SA/LowPass'] # corner frequency of the low-pass filter
    left_index  = ( _n.abs( freqs_calib - (f_lo - (1+pad)*f_c) ) ).argmin()    
    right_index = ( _n.abs( freqs_calib - (f_lo + (1+pad)*f_c) ) ).argmin()    
    print left_index, right_index
    power_baseline = _n.mean( _n.concatenate( (y[:left_index:], y[right_index::] ) ) )
    print ("Baseline = " + str(power_baseline))
    
    # Do the integration over the rest of the array
    p_integrated = _n.sum( y-power_baseline )  * (freqs_calib[1]-freqs_calib[0]) / (2*f_c)
    print p_integrated        
    
    return p_integrated, power_baseline


def fit_dVdI(I, dVdI):
    """
    Fits dV/dI curves to a parabola, getting the base resistance and the heating factor eta
    (R = R_0 + eta I_dc^2)
    """
    
    # Define function... convert I_dc to mA so that scipy does not hate you
    def fit_func(I_dc, R_0, eta):
        return R_0 + eta*(I_dc*1e3)**2
        
    # Get fit and errors    
    popt, pcov = spop.curve_fit(fit_func, I, dVdI, p0 = [dVdI[0], 1.0])
    perr = _n.sqrt(_n.diag(pcov))
    
    # Rescale eta to Ohms/A^2
    popt[1] = popt[1] * 1e6
    perr[1] = perr[1] * 1e6

    return popt, perr
    

### SHUTDOWN PROCEDURE
def shutdown():  
    #print "\n\n ELDERLY COUPLE BEATEN WITH HAMMER \n\n"
    print "Goodbye"
    a.quit_list_mode(force_quit = True)
    a.write("OUTP OFF")
    button_stream.set_checked(False) 

# overwrite the existing one
window.event_close = shutdown
       
# show the window
window.show()