# -*- coding: utf-8 -*-

"""
* @file step.py
* @author Gustavo Diaz H.
* @date 22 May 2020
* @brief Simple test step signal
"""

import numpy as np
from signalbase import SignalBase

class Step(SignalBase):
    def __init__(self, ti, tf, dt, amplitude, delay):
    	SignalBase.__init__(self, ti, tf, dt)
    	self.amplitude = amplitude			#[]
    	self.delay = delay 					#[s]
    	self.i_on = int(delay/dt)			#[i-timestep]
    	self.signal = amplitude*np.ones(self.N)
    	self.signal[0:self.i_on] *= 0

    def __call__(self, *args, **kwargs):
        return self.signal