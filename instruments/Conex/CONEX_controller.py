# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 14:33:37 2015

@author: Ali.Eftekhari
"""
from CONEX_SMC_common import CONEXSMC

class CONEX(CONEXSMC):
    def __init__(self, COMPort):
        """
        Give COM port of device you want as a string, ie
        "COMX", where X = 1,2,3...
        """
        super(CONEX,self).__init__(COMPort)
        self.connect=self.rm.open_resource(COMPort, baud_rate=921600, timeout=2000, data_bits=8, write_termination='\r\n')
        #device_key=raw_input("What is the COM port associated with CONEX controller? ")
        #self.connect=self.rm.open_resource(device_key, baud_rate=921600, timeout=2000, data_bits=8, write_termination='\r\n')

    def get_velocity_feedforward(self):
        return self.write_read("1KV?")
        
    def set_velocity_feedforward(self,value):
        self.value=str(value)
        return self.send("1KV"+self.value)   
        
    def move_absolute_axisnumber(self,axis,value):
        self.axis=str(axis)
        self.value=str(value)
        return self.send(self.axis+"SE"+self.value)
       
    def move_relative(self, axis,value):
        self.axis=str(axis)
        self.value=str(value)
        return self.send(self.axis+"PR"+self.value)