# -*- coding: utf-8 -*-

"""
* @file mtq.py
* @author Gustavo Diaz H.
* @date 22 May 2020
* @brief Simple model of a Magnetorquer actuator
"""

import numpy as np

class MTQ(object):
    def __init__(self, R, L, A, N):
        self.R = R
        self.L = L
        self.A = A
        self.N = N
        self.i = 0
        self.m = 0

    def main_routine(self, V, dt):
    	self.rungeKutta(self.i, V, dt)

    def dynamics(self, i, V):
    	di = V/self.L - (self.R/self.L)*i
    	return di

    def rungeKutta(self, i, V, dt):
        """
        * Runge-Kutta method to integrate the electric current differential equation
        *
        * @param i float Actual electric current state
        * @param V float Actual input Voltage
        * @param dt float integration time step
        * @update i float Next electric current state
        """
        k1 = self.dynamics(i, V)
        k2 = self.dynamics(i + 0.5*dt*k1, V)
        k3 = self.dynamics(i + 0.5*dt*k2, V)
        k4 = self.dynamics(i + dt*k3, V)

        i_next = i + dt*(k1 + 2*(k2+k3)+k4)/6.0
        self.i = i_next
        self.m = self.N * self.i * self.A

if __name__ == '__main__':
    # TEST
    import numpy as np
    import matplotlib.pyplot as plt
    import time
    from step import Step
    from pwm import PWM
    from current_sensor import CurrentSensor

    # MTQ Model
    R = 145.9 		#[Ohm]
    L = 10.08e-3	#[H]
    A = 46.82e-4	#[m2]
    N_coil = 25*10	#[number of turns]
    tau = L/R    	#[s]
    mtq = MTQ(R, L, A, N_coil)

    # Time parameters for simulation
    tf = 1
    dt = 0.1*tau
    N = int(tf/dt)
    t = np.linspace(0, tf, N)

    # Voltage Signal parameters
    #Step Response
    V = 3.3 			#[V]
    delay = 0.2*tf		#[s]
    # voltage = Step(0, tf, dt, V, delay).signal
    #PWM input
    freq = 100			#[Hz]
    duty = 30			#[%]
    voltage = PWM(0, tf, dt, V, freq, duty).signal
    
    # MTQ Data Storage
    i_data = np.zeros(N)
    m_data = np.zeros(N)

    # Sensors
    i_std = 22e-5		#[A]
    i_bias = 22e-6		#[A]
    current_sensor = CurrentSensor(i_bias, i_std)

    for i in range(0, N):
        # Process data
        mtq.main_routine(voltage[i], dt)
        
        # Data of interest
        i_data[i] = current_sensor.measure(mtq.i)
        m_data[i] = mtq.m

    # Data Visualization
    from monitor import Monitor
    
    v_mon = Monitor([t], [voltage], "MTQ input voltage", "V[V]", "time[s]", sig_name = ["V"])
    v_mon.plot()

    i_mon = Monitor([t], [i_data], "MTQ electric current", "i[A]", "time[s]", sig_name = ["i"])
    i_mon.plot()

    m_mon = Monitor([t], [m_data], "MTQ magnetic moment", "m[Am2]", "time[s]", sig_name = ["m"])
    m_mon.plot()
    m_mon.show()