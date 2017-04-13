# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:55:44 2017

@author: admin
"""
import visa as _visa

###############################################
# Anapico class for talking to the MW generator
###############################################
class Anapico():
        
    def __init__(self, device_string="121-227450000-0385", device_string_tcp = "169.254", timeout = 100000):
        """
        Visa-based connection to the Anapico generator.
        """

        # defaults
        self._anapico = None
        self._data    = dict()

        # figure out what's connected        
        self._resource_manager = _visa.ResourceManager()
        self._device_names     = list(self._resource_manager.list_resources()); 
        
        # update the user and find the Anapico
        print "Connected VISA devices:"
        i = -1
        for n in range(len(self._device_names)): 
            
            # for easy coding
            d = self._device_names[n]            
            
            # if it contains the magic string
            if d.find(device_string) >= 0:
                print " ", str(n)+": *"+d
                i = n
            elif d.find(device_string_tcp) >= 0:
                print " ", str(n)+": *"+d
                i = n
            else:
                print " ", str(n)+":  "+d

        # if we didn't find it, bomb out
        if i < 0: 
            print "ERROR: Could not find supplied Anapico string '"+device_string+"' or '" +device_string_tcp+"'"
            return

        print "Connecting to device", i
        self._anapico = self._resource_manager.get_instrument(self._device_names[i])
        print "  Success!"
        
        self._anapico.timeout = timeout

    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        s = "Anapico() instance:"
        if self._anapico == None: return s+" (not connected)"
        
        ks = self._data.keys()
        ks.sort()
        for k in ks: s += "\n  "+k+" = "+self._data[k]
        
        return s


    # Basic commands    
    def write(self, command_string, print_outp=True):
        """
        Safe send to anapico (if connected).
        """
        if self._anapico == None: return "Not connected."
        
        # send the command with proper formatting 
        self._anapico.write(command_string)
        
        # make sure it worked by querying the same parameter
        # this also serves to ensure the anapico is ready 
        # for the next command. I couldn't break this method.
    
        # FIRMWARE        
        q = command_string.split(' ')[0]
        a = self.query(q+"?")
        
        # nicely, the anapico returns an empty string for nonsense queries        
        if len(a) & print_outp: 
           print "    "+q+" = "+a
           self._data[q] = a
        
    def read(self):
        """
        Reads data from the anapico (if connected)
        """
        if self._anapico == None: return "Not connected."
        
        # read it and remove the newline character        
        return self._anapico.read().strip()
    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        if self._anapico == None: return "Not connected."
        
        # use the built-in query (works!) 
        a = self._anapico.query(query_string).strip()
        
        # store this data for reference        
        q = query_string.replace("?","")
        if len(a): self._data[q] = a
        
        return a
     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        self._anapico.write("*RST")
        self._anapico.query("*IDN?")
        self._data = dict()
    
    
    # Complex commands
    def list_man_setup(self,freqs,pows,dwels=0,dels=0):
        """
        Sets up the device for a manual list sweep. Supports both single values
        or lists for powers and frequencies.
        """
      
        # Set power and frequency to fixed mode so we can rewrite the sweep parameters
        self.write("POW:MODE FIX")       
        self.write("FREQ:MODE FIX")

        
        # Find number of steps in the sweep
        numsteps = len(freqs)
       
        # Set the power for all list points
        command = "LIST:POW "        
        try:
            iter(pows)
        except TypeError:
            command += str(pows)            
            for n in range (0,numsteps-1):
                command += ","+str(pows)
        else:
            for p in pows: command += str(p)+","
        self.write(command, False)
      
        # Set the frequency for all list points
        command = "LIST:FREQ "
        try:
            iter(freqs)
        except TypeError:
            command += str(freqs)            
            for n in range (0,numsteps-1):
                command += ","+str(freqs)
            self.write(command)
        else:
            for f in freqs: command += str(f)+","          
        self.write(command, False)
            
        # Set dwell and delay
        self.write("LIST:DWEL " + str(dwels))
        self.write("LIST:DEL " + str(dels))

        # Set sweep to manual mode    
        self.write("LIST:MODE MAN")
        
        # Set power and frequency to either fixed or list mode
        self.write("POW:MODE LIST")       
        self.write("FREQ:MODE LIST")

        # Turn output on          
        self.write("OUTP ON")

    def list_change_index(self, index):
        self.write("LIST:MAN " + str(index+1), print_outp = False)  
        
    def quit_list_mode(self, force_quit = False):
        #if cfg_sweep["PowerMode"] != "Mixed" or force_quit:        
        if force_quit:            
            self.write("FREQ:MODE FIX")
            self.write("POW:MODE FIX")
            self.write("OUTP OFF")
        




###############################################
# SMB100A class for talking to the MW generator
###############################################
class SMB100A_CMC():
        
    def __init__(self, device_string="176142", device_string_tcp = "169.254", timeout = 30000):
        """
        Visa-based connection to the Rohde and Schwarz generator.
        """

        # defaults
        self._smb100a = None
        self._data    = dict()

        # figure out what's connected        
        self._resource_manager = _visa.ResourceManager()
        self._device_names     = list(self._resource_manager.list_resources()); 
        
        # update the user and find the Anapico
        print "Connected VISA devices:"
        i = -1
        for n in range(len(self._device_names)): 
            
            # for easy coding
            d = self._device_names[n]            
            
            # if it contains the magic string
            if d.find(device_string) >= 0:
                print " ", str(n)+": *"+d
                i = n
            elif d.find(device_string_tcp) >= 0:
                print " ", str(n)+": *"+d
                i = n
            else:
                print " ", str(n)+":  "+d

        # if we didn't find it, bomb out
        if i < 0: 
            print "ERROR: Could not find supplied SMB100A string '"+device_string+"' or '" +device_string_tcp+"'"

        print "Connecting to device", i
        self._smb100a = self._resource_manager.get_instrument(self._device_names[i])
        print "  Success!"       
        
        self._smb100a.timeout = timeout

    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        s = "Anapico() instance:"
        if self._smb100a == None: return s+" (not connected)"
        
        ks = self._data.keys()
        ks.sort()
        for k in ks: s += "\n  "+k+" = "+self._data[k]
        
        return s


    # Basic commands    
    def write(self, command_string, print_outp=True):
        """
        Safe send to SMB100A (if connected).
        """
        if self._smb100a == None: return "Not connected."
        
        # send the command with proper formatting 
        self._smb100a.write(command_string)
        
        # make sure it worked by querying the same parameter
        # this also serves to ensure the anapico is ready 
        # for the next command. I couldn't break this method.
    
        # FIRMWARE        
        q = command_string.split(' ')[0]
        a = self.query(q+"?")
        
        # nicely, the anapico returns an empty string for nonsense queries        
        if len(a) & print_outp: 
           print "    "+q+" = "+a
           self._data[q] = a
        
    def read(self):
        """
        Reads data from the SMB100A (if connected)
        """
        if self._smb100a == None: return "Not connected."
        
        # read it and remove the newline character        
        return self._smb100a.read().strip()
    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        if self._smb100a == None: return "Not connected."
        
        # use the built-in query (works!) 
        a = self._smb100a.query(query_string).strip()
        
        # store this data for reference        
        q = query_string.replace("?","")
        if len(a): self._data[q] = a
        
        return a
     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        self._smb100a.write("*RST")
        self._smb100a.query("*IDN?")
        self._data = dict()
    
    
    # Complex commands
    def list_man_setup(self,freqs,pows,dwels=0.001,dels=0):
        """
        Sets up the device for a manual list sweep. Supports both single values
        or lists for powers and frequencies.
        """
      
        # Set power and frequency to fixed mode so we can rewrite the sweep parameters
        self.write("POW:MODE CW", False)       
        self.write("FREQ:MODE CW", False)
     
        # Find number of steps in the sweep
        numsteps = len(freqs)
       
        # Create list
        self.write("LIST:SEL \"test\"", False) 
        self.query("SYST:ERR:ALL?")
      
        # Set the frequency for all list points
        command = "LIST:FREQ "
        try:
            iter(freqs)
        except TypeError:
            command += str(int(freqs))            
            for n in range (0,numsteps-1):
                command += ","+str(int(freqs))
            self.write(command)
        else:
            for f in freqs: command += str(f)+","   
        command = command[:-1]
        self.write(command, False)
        self.query("LIST:FREQ?")    
        
        # Set the power for all list points
        command = "LIST:POW "        
        try:
            iter(pows)
        except TypeError:
            command += str(pows)            
            for n in range (0,numsteps-1):
                command += ","+str(pows)
        else:
            for p in pows: command += str(p)+","
        command = command[:-1]
        self.write(command, False)
        self.query("LIST:POW?")     
            
        # Set dwell and delay
        self.write("LIST:DWEL " + str(dwels))

        # Set sweep to manual mode    
        self.write("LIST:MODE STEP")
        
        # Turn output on
        self.write("OUTP 1")        
        
        # Set power and frequency to either fixed or list mode
        self.write("FREQ:MODE LIST")
        self.write("POW:MODE LIST")       
        
        # Turn output on       
        self.write("OUTP ON")
        
    def list_change_index(self, index):
        self.write("LIST:IND " + str(index))

    def quit_list_mode(self, force_quit = False):
        self.write("FREQ:MODE CW")
        self.write("POW:MODE CW")
        #if cfg_sweep["PowerMode"] != "Mixed" or force_quit:
        if force_quit:
            self.write("OUTP OFF")




###############################################
# SMB100A class for talking to the MW generator
###############################################
class SMB100A_RnS():
        
    def __init__(self, device_string="178288", device_string_tcp = "169.254", timeout = 100000):
        """
        Visa-based connection to the Rohde and Schwarz generator.
        """

        # defaults
        self._smb100a = None
        self._data    = dict()

        # figure out what's connected        
        self._resource_manager = _visa.ResourceManager()
        self._device_names     = list(self._resource_manager.list_resources()); 
        
        # update the user and find the Anapico
        print "Connected VISA devices:"
        i = -1
        for n in range(len(self._device_names)): 
            
            # for easy coding
            d = self._device_names[n]            
            
            # if it contains the magic string
            if d.find(device_string) >= 0:
                print " ", str(n)+": *"+d
                i = n
            elif d.find(device_string_tcp) >= 0:
                print " ", str(n)+": *"+d
                i = n
            else:
                print " ", str(n)+":  "+d

        # if we didn't find it, bomb out
        if i < 0: 
            print "ERROR: Could not find supplied SMB100A string '"+device_string+"' or '" +device_string_tcp+"'"
            return

        print "Connecting to device", i
        self._smb100a = self._resource_manager.get_instrument(self._device_names[i])
        print "  Success!"
        
        self._smb100a.timeout = timeout

    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        s = "Anapico() instance:"
        if self._smb100a == None: return s+" (not connected)"
        
        ks = self._data.keys()
        ks.sort()
        for k in ks: s += "\n  "+k+" = "+self._data[k]
        
        return s


    # Basic commands    
    def write(self, command_string, print_outp=True):
        """
        Safe send to SMB100A (if connected).
        """
        if self._smb100a == None: return "Not connected."
        
        # send the command with proper formatting 
        self._smb100a.write(command_string)
        
        # make sure it worked by querying the same parameter
        # this also serves to ensure the anapico is ready 
        # for the next command. I couldn't break this method.
    
        # FIRMWARE        
        q = command_string.split(' ')[0]
        a = self.query(q+"?")
        
        # nicely, the anapico returns an empty string for nonsense queries        
        if len(a) & print_outp: 
           print "    "+q+" = "+a
           self._data[q] = a
        
    def read(self):
        """
        Reads data from the SMB100A (if connected)
        """
        if self._smb100a == None: return "Not connected."
        
        # read it and remove the newline character        
        return self._smb100a.read().strip()
    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        if self._smb100a == None: return "Not connected."
        
        # use the built-in query (works!) 
        a = self._smb100a.query(query_string).strip()
        
        # store this data for reference        
        q = query_string.replace("?","")
        if len(a): self._data[q] = a
        
        return a
     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        self._smb100a.write("*RST")
        self._smb100a.query("*IDN?")
        self._data = dict()
    
    
    # Complex commands
    def list_man_setup(self,freqs,pows,dwels=0.001,dels=0):
        """
        Sets up the device for a manual list sweep. Supports both single values
        or lists for powers and frequencies.
        """
      
        # Set power and frequency to fixed mode so we can rewrite the sweep parameters
        self.write("POW:MODE CW", False)       
        self.write("FREQ:MODE CW", False)
     
        # Find number of steps in the sweep
        numsteps = len(freqs)
       
        # Create list
        self.write("LIST:SEL \"test\"", False) 
        self.query("SYST:ERR:ALL?")
      
        # Set the frequency for all list points
        command = "LIST:FREQ "
        try:
            iter(freqs)
        except TypeError:
            command += str(int(freqs))            
            for n in range (0,numsteps-1):
                command += ","+str(int(freqs))
            self.write(command)
        else:
            for f in freqs: command += str(f)+","   
        command = command[:-1]
        self.write(command, False)
        self.query("LIST:FREQ?")    
        
        # Set the power for all list points
        command = "LIST:POW "        
        try:
            iter(pows)
        except TypeError:
            command += str(pows)            
            for n in range (0,numsteps-1):
                command += ","+str(pows)
        else:
            for p in pows: command += str(p)+","
        command = command[:-1]
        self.write(command, False)
        self.query("LIST:POW?")     
            
        # Set dwell and delay
        self.write("LIST:DWEL " + str(dwels))

        # Set sweep to manual mode    
        self.write("LIST:MODE STEP")
        
        # Turn output on
        self.write("OUTP 1")        
        
        # Set power and frequency to either fixed or list mode
        self.write("FREQ:MODE LIST")
        self.write("POW:MODE LIST")       
        
        # Turn output on       
        self.write("OUTP ON")
        
    def list_change_index(self, index):
        self.write("LIST:IND " + str(index))

    def quit_list_mode(self, force_quit = False):
        self.write("FREQ:MODE CW")
        self.write("POW:MODE CW")
        if force_quit:
        #if cfg_sweep["PowerMode"] != "Mixed" or force_quit:
            self.write("OUTP OFF")



###############################################
# SMB100A class for talking to the MW generator
###############################################
class SMB100A_6GHz():
        
    def __init__(self, device_string="113055", device_string_tcp = "169.254", timeout = 100000):
        """
        Visa-based connection to the Rohde and Schwarz generator.
        """

        # defaults
        self._smb100a = None
        self._data    = dict()

        # figure out what's connected        
        self._resource_manager = _visa.ResourceManager()
        self._device_names     = list(self._resource_manager.list_resources()); 
        
        # update the user and find the Anapico
        print "Connected VISA devices:"
        i = -1
        for n in range(len(self._device_names)): 
            
            # for easy coding
            d = self._device_names[n]            
            
            # if it contains the magic string
            if d.find(device_string) >= 0:
                print " ", str(n)+": *"+d
                i = n
            elif d.find(device_string_tcp) >= 0:
                print " ", str(n)+": *"+d
                i = n
            else:
                print " ", str(n)+":  "+d

        # if we didn't find it, bomb out
        if i < 0: 
            print "ERROR: Could not find supplied SMB100A string '"+device_string+"' or '" +device_string_tcp+"'"
            return

        print "Connecting to device", i
        self._smb100a = self._resource_manager.get_instrument(self._device_names[i])
        print "  Success!"
        
        self._smb100a.timeout = timeout

    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        s = "Anapico() instance:"
        if self._smb100a == None: return s+" (not connected)"
        
        ks = self._data.keys()
        ks.sort()
        for k in ks: s += "\n  "+k+" = "+self._data[k]
        
        return s


    # Basic commands    
    def write(self, command_string, print_outp=True):
        """
        Safe send to SMB100A (if connected).
        """
        if self._smb100a == None: return "Not connected."
        
        # send the command with proper formatting 
        self._smb100a.write(command_string)
        
        # make sure it worked by querying the same parameter
        # this also serves to ensure the anapico is ready 
        # for the next command. I couldn't break this method.
    
        # FIRMWARE        
        q = command_string.split(' ')[0]
        a = self.query(q+"?")
        
        # nicely, the anapico returns an empty string for nonsense queries        
        if len(a) & print_outp: 
           print "    "+q+" = "+a
           self._data[q] = a
        
    def read(self):
        """
        Reads data from the SMB100A (if connected)
        """
        if self._smb100a == None: return "Not connected."
        
        # read it and remove the newline character        
        return self._smb100a.read().strip()
    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        if self._smb100a == None: return "Not connected."
        
        # use the built-in query (works!) 
        a = self._smb100a.query(query_string).strip()
        
        # store this data for reference        
        q = query_string.replace("?","")
        if len(a): self._data[q] = a
        
        return a
     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        self._smb100a.write("*RST")
        self._smb100a.query("*IDN?")
        self._data = dict()
    
    
    # Complex commands
    def list_man_setup(self,freqs,pows,dwels=0.001,dels=0):
        """
        Sets up the device for a manual list sweep. Supports both single values
        or lists for powers and frequencies.
        """
      
        # Set power and frequency to fixed mode so we can rewrite the sweep parameters
        self.write("POW:MODE CW", False)       
        self.write("FREQ:MODE CW", False)
     
        # Find number of steps in the sweep
        numsteps = len(freqs)
       
        # Create list
        self.write("LIST:SEL \"test\"", False) 
        self.query("SYST:ERR:ALL?")
      
        # Set the frequency for all list points
        command = "LIST:FREQ "
        try:
            iter(freqs)
        except TypeError:
            command += str(int(freqs))            
            for n in range (0,numsteps-1):
                command += ","+str(int(freqs))
            self.write(command)
        else:
            for f in freqs: command += str(f)+","   
        command = command[:-1]
        self.write(command, False)
        self.query("LIST:FREQ?")    
        
        # Set the power for all list points
        command = "LIST:POW "        
        try:
            iter(pows)
        except TypeError:
            command += str(pows)            
            for n in range (0,numsteps-1):
                command += ","+str(pows)
        else:
            for p in pows: command += str(p)+","
        command = command[:-1]
        self.write(command, False)
        self.query("LIST:POW?")     
            
        # Set dwell and delay
        self.write("LIST:DWEL " + str(dwels))

        # Set sweep to manual mode    
        self.write("LIST:MODE STEP")
        self.write("LIST:TRIG:SOUR SING")
        
        # Turn output on
        self.write("OUTP 1")        
        
        # Set power and frequency to either fixed or list mode
        self.write("FREQ:MODE LIST")   
        
        # Turn output on       
        self.write("OUTP ON")
        
    def list_change_index(self, index):
        self.write("LIST:IND " + str(index))

    def quit_list_mode(self, force_quit = False):
        self.write("FREQ:MODE CW")
        self.write("POW:MODE CW")
        if force_quit:
        #if cfg_sweep["PowerMode"] != "Mixed" or force_quit:
            self.write("OUTP OFF")


###############################################
# NoneGen class acting as a dummy generator
###############################################
class NoneGen():
        
    def __init__(self, device_string="113055", device_string_tcp = "169.254", timeout = 100000):
        """
        Visa-based connection to the Rohde and Schwarz generator.
        """
        print "No signal generator selected"

    def __repr__(self):
        """
        Returns a python representation string for the object."
        """
        print "No signal generator selected"


    # Basic commands    
    def write(self, command_string, print_outp=True):
        """
        Safe send to SMB100A (if connected).
        """
        print "No signal generator selected"
        
    def read(self):
        """
        Reads data from the SMB100A (if connected)
        """
        print "No signal generator selected"
    
    def query(self, query_string):
        """
        Sends the supplied command and waits for data from the device (if connected).
        """
        print "No signal generator selected"

     
    def reset(self):
        """
        Resets the device to factory defaults (and clears the internal memory
        of this python object).
        """
        print "No signal generator selected"
    
    
    # Complex commands
    def list_man_setup(self,freqs,pows,dwels=0.001,dels=0):
        print "No signal generator selected"
        
    def list_change_index(self, index):
        print "No signal generator selected"

    def quit_list_mode(self, force_quit = False):
        print "No signal generator selected"