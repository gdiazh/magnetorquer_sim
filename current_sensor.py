# -*- coding: utf-8 -*-

"""
* @file current_sensor.py
* @author Gustavo Diaz H.
* @date 23 May 2020
* @brief Simple current sensor model
"""

import numpy as np

class CurrentSensor(object):
    def __init__(self, bias, nrs_std):
        self.std = nrs_std                          #[A]
        self.bias = bias                            #[A]
        self.nrs = np.random.normal(0, self.std)    #[A]

    def measure(self, true):
        self.nrs = np.random.normal(0, self.std)    #[A]
        measurement = true + self.bias + self.nrs
        return measurement