# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:32:50 2015

@author: CHILDRESSLAB

Attempting to replicate the Following Errors. Note that the CONEX_controller,
CONEX_SMC_common and motrion_controller_main libraries have been slightly
altered to allow for each object to be constructed based on and input COMport.
"""

from CONEX_controller import CONEX
import time
import sys

X = CONEX("COM5")
Y = CONEX("COM3")

Xmin = 2.26994357
Ymin = 20.3500223

#Typically, I have found this command causes a following error with
#one of my controllers but not the other. This may or may not cause the
#controller to enter "disabled from moving" state from a Following Error
Y.move_absolute(1, Ymin)

#Force program to wait until actuator gets to position 
while Y.get_enable_disable_state(1) == u'1MM28\r\n':
    print Y.get_enable_disable_state(1)   
    time.sleep(0.1)

#Exit program if a problem occured
if Y.get_enable_disable_state(1) != u'1MM33\r\n':
    print "Exited since Y = "+Y.get_enable_disable_state(1) 
    sys.exit()


X.move_absolute(1, Xmin)

#Force program to wait until actuator gets to position 
while (X.get_enable_disable_state(1) == u'1MM28\r\n'):
    print X.get_enable_disable_state(1)    
    time.sleep(0.1)

#Exit program if a problem occured
if X.get_enable_disable_state(1) != u'1MM33\r\n':
    print "Exited since X = "+X.get_enable_disable_state(1) 
    sys.exit()
    
    
#More reproducibly, sending commands to two controllers simultaneously tends
#to cause the second to become disabled from moving due to a following error.
X.move_absolute(1,25)
Y.move_absolute(1,25)

X.get_enable_disable_state(1)
Y.get_enable_disable_state(1)