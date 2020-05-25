# -*- coding: utf-8 -*-

"""
* @file duty_test.py
* @author Gustavo Diaz H.
* @date 24 May 2020
* @brief Test: Mesurement of magnetic field at different PWM duty cicles of input voltage
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtq import MTQ
from step import Step
from pwm import PWM
from current_sensor import CurrentSensor
from magnetometer import Magnetometer

# MTQ Model
R = 145.9 		#[Ohm]
L = 10.08e-3	#[H]
A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m]
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
#PWM input
freq = 100          #[Hz]
d0 = 0              #[%]
dd = 10             #[%]
df = 90             #[%]
Nd = int((df-d0)/dd+1) #[number of points]
duty = np.linspace(d0, df, Nd) #[%]
voltage = np.zeros((Nd, N))
v_rms = np.zeros(Nd)
v_med = np.zeros(Nd)
for i in range(0, Nd):
    voltage[i,:] = PWM(0, tf, dt, V, freq, duty[i]).signal
    v_rms[i] = np.sqrt(np.mean(voltage[i,:]**2))
    v_med[i] = np.mean(voltage[i,:])

# Sensors
# Electric current
i_std = 22e-5		#[A]
i_bias = 22e-6		#[A]
current_sensor = CurrentSensor(i_bias, i_std)
# Magnetometer
B_std = 1			#[uT]
B_bias = 1e-1		#[uT]
magnetometer = Magnetometer(B_bias, B_std)
z = 16e-3			#[m]

# MTQ Data Storage
i_data = np.zeros((Nd, N))
m_data = np.zeros((Nd, N))
B_data = np.zeros((Nd, N))
m_calc = np.zeros((Nd, N))

i_rms = np.zeros(Nd)
m_rms = np.zeros(Nd)
B_rms = np.zeros(Nd)
m_calc_rms = np.zeros(Nd)

i_med = np.zeros(Nd)
m_med = np.zeros(Nd)
B_med = np.zeros(Nd)
m_calc_med = np.zeros(Nd)

# Run test
for j in range(0, Nd):
    print("PWM duty cycle: ", duty[j], "\t", j+1, "/", Nd)
    for i in range(0, N):
        # Process data
        mtq.main_routine(voltage[j][i], dt)
        
        # Data of interest
        i_data[j][i] = current_sensor.measure(mtq.i)
        m_data[j][i] = mtq.m
        B_data[j][i] = magnetometer.measure(z, mtq.i, l, N_coil)
        m_calc[j][i] = (5)*B_data[j][i]*(z**2+(0.5*l)**2)*np.sqrt(z**2+0.5*l**2)

    # rms values
    i_rms[j] = np.sqrt(np.mean(i_data[j]**2))
    m_rms[j] = np.sqrt(np.mean(m_data[j]**2))
    B_rms[j] = np.sqrt(np.mean(B_data[j]**2))
    m_calc_rms[j] = np.sqrt(np.mean(m_calc[j]**2))

    # med values
    i_med[j] = np.mean(i_data[j])
    m_med[j] = np.mean(m_data[j])
    B_med[j] = np.mean(B_data[j])
    m_calc_med[j] = np.mean(m_calc[j])

    # reset state
    mtq.reset()

# Save data to csv
from pathlib import Path
import datetime
folder = "data/osc/"
Path(folder).mkdir(parents=True, exist_ok=True)
data = {"time[s]": t, "voltage[V]": voltage[0], "current[A]": i_data[0]}
df = pd.DataFrame(data, columns=["time[s]", "voltage[V]", "current[A]"])
test_name = "-test_pwm[voltage-current]"
date = datetime.datetime.now().strftime('%Y-%m-%d %H-%M-%S')
path = Path(folder+date+test_name+".csv")
df.to_csv(path)

folder = "data/magnetometer/"+date+"/"
Path(folder).mkdir(parents=True, exist_ok=True)
for j in range(0, Nd):
    data = {"time[s]": t, "magField[uT]": B_data[j]}
    df = pd.DataFrame(data, columns=["time[s]", "magField[uT]"])
    test_name = "-test_pwm[magField]["+str(int(duty[j]))+"%]"
    path = Path(folder+test_name+".csv")
    df.to_csv(path)

# Data Visualization
from monitor import Monitor

v_mon = Monitor([t, t, t, t], [voltage[0], voltage[1], voltage[2], voltage[3]], "MTQ input voltage", "V[V]", "time[s]", sig_name = ["V0", "V1", "V2", "V3"])
v_mon.plot()

i_mon = Monitor([t, t, t, t], [i_data[0], i_data[1], i_data[2], i_data[3]], "MTQ electric current", "i[A]", "time[s]", sig_name = ["i0", "i1", "i2", "i3"])
i_mon.plot()

B_mon = Monitor([t, t, t, t], [B_data[0], B_data[1], B_data[2], B_data[3]], "MTQ magnetic field", "B[uT]", "time[s]", sig_name = ["B0", "B1", "B2", "B3"])
B_mon.plot()

mCalc_mon = Monitor([t, t, t, t], [m_calc[0], m_calc[1], m_calc[2], m_calc[3]], "MTQ calculated magnetic moment", "m[Am2]", "time[s]", sig_name = ["m0", "m1", "m2", "m3"])
mCalc_mon.plot()

m_mon = Monitor([t, t, t, t], [m_data[0], m_data[1], m_data[2], m_data[3]], "MTQ magnetic moment", "m[Am2]", "time[s]", sig_name = ["m0", "m1", "m2", "m3"])
m_mon.plot()

vrms_mon = Monitor([duty], [v_rms], "MTQ input rms voltage", "V[Vrms]", "time[s]", sig_name = ["Vrms"])
vrms_mon.plot()

vmed_mon = Monitor([duty], [v_med], "MTQ input med voltage", "V[Vmed]", "time[s]", sig_name = ["Vmed"])
vmed_mon.plot()

irms_mon = Monitor([duty], [i_rms], "MTQ input rms current", "i[irms]", "time[s]", sig_name = ["irms"])
irms_mon.plot()

Brms_mon = Monitor([duty], [B_rms], "MTQ magnetic field rms", "B[Brms]", "time[s]", sig_name = ["Brms"])
Brms_mon.plot()

mCalcrms_mon = Monitor([duty], [m_calc_rms], "MTQ calc magnetic moment rms", "mCalc[mrms]", "time[s]", sig_name = ["mrms"])
mCalcrms_mon.plot()

mrms_mon = Monitor([duty], [m_rms], "MTQ magnetic moment rms", "m[mrms]", "time[s]", sig_name = ["mrms"])
mrms_mon.plot()

imed_mon = Monitor([duty], [i_med], "MTQ input med current", "i[imed]", "time[s]", sig_name = ["imed"])
imed_mon.plot()

Bmed_mon = Monitor([duty], [B_med], "MTQ magnetic field med", "B[Bmed]", "time[s]", sig_name = ["Bmed"])
Bmed_mon.plot()

mCalcmed_mon = Monitor([duty], [m_calc_med], "MTQ calc magnetic moment med", "mCalc[mmed]", "time[s]", sig_name = ["mmed"])
mCalcmed_mon.plot()

mmed_mon = Monitor([duty], [m_med], "MTQ magnetic moment med", "m[mmed]", "time[s]", sig_name = ["mmed"])
mmed_mon.plot()
mmed_mon.show()