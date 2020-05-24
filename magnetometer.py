# -*- coding: utf-8 -*-

"""
* @file magnetometer.py
* @author Gustavo Diaz H.
* @date 23 May 2020
* @brief Simple magnetometer model
"""

import numpy as np

class Magnetometer(object):
    def __init__(self, bias, nrs_std):
        self.std = nrs_std                          #[A]
        self.bias = bias                            #[A]
        self.nrs = np.random.normal(0, self.std)    #[A]

    def measure(self, z, I_true, S, N):
        Bz_true = self.calc(z, I_true, S, N)        #[uT]
        self.nrs = np.random.normal(0, self.std)    #[A]
        measurement = Bz_true + self.bias + self.nrs
        return measurement

    def calc(self, z, I, S, N):
        u0 = 4*np.pi*10e-7      #[H/m]
        Bz = (u0*I*N*S**2)/(2*np.pi*(z**2+(0.5*S)**2)*np.sqrt(z**2+0.5*(S)**2))
        return Bz*1e6   #[uT]