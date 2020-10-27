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

    def measure_circular(self, z, I_true, b, N):
        Bz_true = self.calc_circular(z, I_true, b, N)   #[uT]
        self.nrs = np.random.normal(0, self.std)        #[A]
        measurement = Bz_true + self.bias + self.nrs
        return measurement

    def measure_rod(self, z, I_true, b, N, L, A):
        Bz_true = self.calc_rod(z, I_true, b, N, L, A)   #[uT]
        self.nrs = np.random.normal(0, self.std)         #[A]
        measurement = Bz_true + self.bias + self.nrs
        return measurement

    def getMagneticMoment(self, z, I, b, N, L, A):
        R = z + L/2         #[m]
        u0 = 4*np.pi*1e-7   #[H/m]
        a1 = (R/L - 0.5)/((R**2-R*L+0.25*L**2)**1.5)    #[m^-3]
        a2 = (R/L + 0.5)/((R**2+R*L+0.25*L**2)**1.5)    #[m^-3]
        dMz = (1/(a1 - a2))                             #[m^3]
        return dMz                                      #[]

    def calc(self, z, I, S, N):
        u0 = 4*np.pi*1e-7      #[H/m]
        Bz = (u0*I*N*S**2)/(2*np.pi*(z**2+(0.5*S)**2)*np.sqrt(z**2+0.5*(S)**2))
        return Bz*1e6   #[uT]

    def calc_circular(self, z, I, b, N):
        u0 = 4*np.pi*1e-7      #[H/m]
        Bz = (0.5*u0*I*N)*((b**2)/((b**2+z**2)**1.5))
        return Bz*1e6   #[uT]

    def calc_rod(self, z, I, b, N, L, A):
        R = z + L/2         #[m]
        u0 = 4*np.pi*1e-7   #[H/m]
        a1 = (R/L - 0.5)/((R**2-R*L+0.25*L**2)**1.5)    #[m^-3]
        a2 = (R/L + 0.5)/((R**2+R*L+0.25*L**2)**1.5)    #[m^-3]
        Bz = (0.25*u0*I*N*A/np.pi)*(a1 - a2)            #[T]
        return Bz*1e6       #[uT]