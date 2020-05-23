# -*- coding: utf-8 -*-

"""
* @file pwm.py
* @author Gustavo Diaz H.
* @date 22 May 2020
* @brief Simple test pwm signal
"""

import numpy as np
from signalbase import SignalBase

class PWM(SignalBase):
    def __init__(self, ti, tf, dt, amplitude, freq, duty):
    	SignalBase.__init__(self, ti, tf, dt)
    	self.amplitude = amplitude			#[]
    	self.freq = freq 					#[Hz]
    	self.T = 1.0/freq 					#[s]
    	self.duty = duty 					#[%]
    	self.Doff = (1-duty/100.0)/freq
    	self.signal = amplitude*np.zeros(self.N)
    	self.generate()

    def __call__(self, *args, **kwargs):
        return self.signal

    def generate(self):
    	t0 = self.ti
    	for i in range(0, self.N):
    		dt = self.time[i]-t0
    		if dt>=self.Doff and dt<self.T:
    			self.signal[i] = self.amplitude
    		elif dt>=self.T:
    			t0 = self.time[i]