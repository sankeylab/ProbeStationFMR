# -*- coding: utf-8 -*-
"""
Created on Tue Oct 06 19:19:14 2015

@author: CHILDRESSLAB

Requires that files CONEX_controller.py, CONEX_SMC_common.py,
and motrion_controller_main.py be located in the same folder as
this script.

Use autocompletion to view possible function, or view them under
the motion_controller_main.py script.
"""

from CONEX_controller import CONEX
#from numpy import matrix
#import math
import time

#Initialize the x,y,z-axis TRB actuators by feeding the relevant COM port ("COM3")
X = CONEX("COM5")
Y = CONEX("COM3")
Z = CONEX("COM4")

#Some useful commands:

#Gets the firmware version (must specify for axis 1)
#print X.get_firmware_version(1)
#print Y.get_firmware_version(1)
#print Z.get_firmware_version(1)

def SetState(IO):
    """
    Set each axis state to IO = 1 == "enable" (or IO = 0 == "disable")
    """
    X.set_enable_disable_state(1,IO)
    Y.set_enable_disable_state(1,IO)
    Z.set_enable_disable_state(1,IO)

def GetPosition():
    """
    Returns a list of three floats representing the (x,y,z) vector
    """
    xpos = float(X.get_current_position(1)[3:])
    ypos = float(Y.get_current_position(1)[3:])
    zpos = float(Z.get_current_position(1)[3:])
    return (xpos, ypos, zpos)

def SetPosition(r):
    """
    Assumes a 3 entry array for r = (x,y,z)
    """
    X.move_absolute(1, r[0])
    Y.move_absolute(1, r[1])
    Z.move_absolute(1, r[2])

def GetState():
    """
    Returns the state of the XYZ actuators
    """
    state = X.get_enable_disable_state(1)+Y.get_enable_disable_state(1)\
        +Z.get_enable_disable_state(1)
    return state
        

def WaitForMove():
    """
    Waits for the actuators to stop moving (ie, all actuators out of 1MM28
    moving state).
    Returns True if upon completion of the move all actuators are in command
    ready state (ie, 1MM33).
    Returns False if any of the actuator is not in the command ready state
    (use GetState to figure out which actuator is in which state).
    """
    while (X.get_enable_disable_state(1) == u'1MM28\r\n' \
    and Y.get_enable_disable_state(1) == u'1MM28\r\n' \
    and Z.get_enable_disable_state(1) == u'1MM28\r\n'):
        #print (GetState())
        time.sleep(0.1)

    #Exit program if a problem occured
    if (X.get_enable_disable_state(1) != u'1MM33\r\n' \
    and Y.get_enable_disable_state(1) != u'1MM33\r\n' \
    and Z.get_enable_disable_state(1) != u'1MM33\r\n'):
        #print (GetState())
        return False
    else:
        return True

