# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 12:20:09 2016

@author: admin
"""

import numpy              as _n
import time               as _time 
import visa               as _visa
from visa import constants


class LakeShore121():
        
    def __init__(self, device_string="ASRL3"):
        """
        Visa-based connection to the ls121 generator.
        """

        # defaults
        self._ls121 = None
        self._data    = dict()

        # figure out what's connected        
        self._resource_manager = _visa.ResourceManager()
        self._device_names     = list(self._resource_manager.list_resources()); 
        
        # update the user and find the ls121
        print "Connected VISA devices:"
        i = -1
        for n in range(len(self._device_names)): 
            
            # for easy coding
            d = self._device_names[n]            
            
            # if it contains the magic string
            if d.find(device_string) >= 0:
                print " ", str(n)+": *"+d
                i = n
            else:
                print " ", str(n)+":  "+d

        # if we didn't find it, bomb out
        if i < 0: 
            print "ERROR: Could not find supplied Lake Shore 121 string '"+device_string+"'"
            return

        print "Connecting to device", i
        self._ls121 = self._resource_manager.get_instrument(self._device_names[i])
        print "  Success!"
        
        # Set baud rate and other subtleties        
        self._ls121.baud_rate = 57600
        self._ls121.data_bits = 7
        self._ls121.parity = constants.Parity.odd
        self._ls121.stop_bits = constants.StopBits.one
        
        # Set to custom current mode
        self._ls121.write("RANGE 13")


    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        s = "lakeshore121() instance:"
        if self._ls121 == None: return s+" (not connected)"
        
        ks = self._data.keys()
        ks.sort()
        for k in ks: s += "\n  "+k+" = "+self._data[k]
        
        return s
        
        
    def write(self, command_string):
        """
        Safe send to ls121 (if connected).
        """
        if self._ls121 == None: return "Not connected."
        
        # send the command        
        self._ls121.write(command_string)
        
        # make sure it worked by querying the same parameter
        # this also serves to ensure the ls121 is ready 
        # for the next command. I couldn't break this method.
        q = command_string.split(' ')[0]
        a = self.query(q+"?")
        
        # nicely, the ls121 returns an empty string for nonsense queries        
        if len(a): 
            print "    "+q+" = "+a
            self._data[q] = a
            
        print "am I sleeping or what"    
        _time.sleep(0.3)
      
      
    def read(self):
        """
        Reads data from the ls121 (if connected)
        """
        if self._ls121 == None: return "Not connected."
        
        # read it and remove the newline character        
        return self._ls121.read().strip()
        
        _time.sleep(0.3)
    
    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        if self._ls121 == None: return "Not connected."
        
        # use the built-in query (works!) 
        a = self._ls121.query(query_string).strip()
        
        # store this data for reference        
        q = query_string.replace("?","")
        if len(a): self._data[q] = a
            
        _time.sleep(0.3)
        
        return a
     
     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        self._ls121.write("*RST")
        self._ls121.query("*IDN?")
        self._data = dict()        
        
        
    def set_current(self, I):
        """
        Set current (if connected).
        """
        if self._ls121 == None: return "Not connected."
        
        # Get current polarity
        polarity = 0; current = abs(I)
        if I<0: polarity = 1
        
        # send the command      
        self._ls121.write("IPOL " + str(polarity))
        self._ls121.write("SETI " + str(current))
       
    def turn_on(self):
        """
        Turns on the ls121 output
        """
        self._ls121.write("IENBL 1")        
      
    def turn_off(self):
        """
        Turns off the ls121 output
        """
        self._ls121.write("IENBL 0")


ls121 = LakeShore121()

def current_sweep():
    Is = _n.linspace(-10e-6, 10e-6, 21) 
    
    ls121.turn_on()    
    
    for k in range(len(Is)):
        print k
        ls121.set_current(Is[k])
        
    ls121.turn_off()