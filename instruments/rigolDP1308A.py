# -*- coding: utf-8 -*-
"""
Created on Wed Feb 03 15:46:25 2016

@author: admin
"""

import numpy              as _n
import time               as _time 
import visa               as _visa
from visa import constants

class Rigol():
        
    def __init__(self, device_string="DP1A140300033"):
        """
        Visa-based connection to the rigol power supply.
        """

        # defaults
        self._rigol = None
        self._data    = dict()

        # figure out what's connected        
        self._resource_manager = _visa.ResourceManager()
        self._device_names     = list(self._resource_manager.list_resources()); 
        
        # update the user and find the rigol
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
        self._rigol = self._resource_manager.get_instrument(self._device_names[i])
        print "  Success!"

        self.write("OUTP:OCP P25V, 0.1")

    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        s = "lakeshore121() instance:"
        if self._rigol == None: return s+" (not connected)"
        
        ks = self._data.keys()
        ks.sort()
        for k in ks: s += "\n  "+k+" = "+self._data[k]
        
        return s
        
        
    def write(self, command_string):
        """
        Safe send to rigol (if connected).
        """
        if self._rigol == None: return "Not connected."
        
        # send the command        
        self._rigol.write(command_string)
        
        # make sure it worked by querying the same parameter
        # this also serves to ensure the rigol is ready 
        # for the next command. I couldn't break this method.
        q = command_string.split(' ')[0]
        a = self.query(q+"?")
        
        # nicely, the rigol returns an empty string for nonsense queries        
        if len(a): 
            print "    "+q+" = "+a
            self._data[q] = a
      
      
    def read(self):
        """
        Reads data from the rigol (if connected)
        """
        if self._rigol == None: return "Not connected."
        
        # read it and remove the newline character        
        return self._rigol.read().strip()
        
        _time.sleep(0.3)

    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        if self._rigol == None: return "Not connected."
        
        # use the built-in query (works!) 
        a = self._rigol.query(query_string).strip()
        
        # store this data for reference        
        q = query_string.replace("?","")
        if len(a): self._data[q] = a
        
        return a
        
     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        self._rigol.write("*RST")
        self._rigol.query("*IDN?")
        self._data = dict()      
        
  
    # Commands specific to the instrument
      
    def turn_on(self, chan):
        self.write("OUTP " + chan + ",ON")
        
    def turn_off(self, chan):
        self.write("OUTP " + chan + ",OFF")        
        
    def set_ocp(self, chan, cur = 0.1):
        self.write("OUTP:OCP " + chan + "," + str(cur))
        self.write("OUTP:OCP:STAT " + chan + ",ON")

    def set_ovp(self, chan, volt = None):
        if volt == None and chan in ("P25V, N25V"): volt = 25
        elif volt == None and chan == "P6V": volt = 6     
        self.write("OUTP:OVP " + chan + "," + str(volt))
        self.write("OUTP:OVP:STAT " + chan + ",ON")        
        
    def set_current(self, chan, cur, volt = None):
        if volt == None and chan in ("P25V, N25V"): volt = 25
        elif volt == None and chan == "P6V": volt = 6      
        self.write("APPL " + chan + ", " + str(volt) + "," + str(cur))
        
    def measure_current(self, chan = None):
        if chan != None: chan = " " + chan
        cur = self.query("MEAS:CURR?" + chan)
        return cur 
        
    def measure_voltage(self, chan = None):
        if chan != None: chan = " " + chan
        volt = self.query("MEAS:VOLT?" + chan)
        return volt    
        
    def measure_power(self, chan = None):
        if chan != None: chan = " " + chan
        powe = self.query("MEAS:POWE?" + chan)
        return powe       