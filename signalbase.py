# -*- coding: utf-8 -*-

"""
* @file signal.py
* @author Gustavo Diaz H.
* @date 22 May 2020
* @brief Simple test signals handler
"""

import numpy as np

class SignalBase(object):
    def __init__(self, ti, tf, dt):
    	self.ti = ti
    	self.tf = tf
    	self.dt = dt
    	self.N = int(tf/dt)
    	self.time = np.linspace(ti, tf, self.N)

    def function():
    	pass