# -*- coding: utf-8 -*-

"""
* @file distance_test.py
* @author Gustavo Diaz H.
* @date 24 May 2020
* @brief Test: Mesurement of magnetic field at different distances
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
R = 28.7                #[Ohm]
L = 78.2e-3             #[H]
d = 8.96e-3             #[m]
l = 70.12e-3            #[m]
A = np.pi*(d/2)**2      #[m2]
e = 0.45e-3             #[m] wire thickness
tau = L/R               #[s]
L_rod = 20e-3           #[m]
r_rod = 5e-3            #[m]
N_coil = (L_rod/e)      #[number of turns]
A_rod = np.pi*r_rod**2  #[m2]
dz0 = 2e-3              #[m]
mtq = MTQ(R, L, A, N_coil)

# Time parameters for simulation
tf = 1
dt = 0.1*tau
N = int(tf/dt)
t = np.linspace(0, tf, N)

# Voltage Signal parameters
#Step Response
V = 3.3             #[V]
delay = 0.2*tf      #[s]
voltage = Step(0, tf, dt, V, delay).signal

# Sensors
# Electric current
i_std = 22e-5       #[A]
i_bias = 22e-6      #[A]
current_sensor = CurrentSensor(i_bias, i_std)
# Magnetometer
B_std = 1e-2        #[uT]
B_bias = 1e-3       #[uT]
magnetometer = Magnetometer(B_bias, B_std)
z0 = 0             #[mm]
dz = 1             #[mm]
zf = 40            #[mm]
Nz = int((zf-z0+1)/dz) #[number of points]
z = np.linspace(z0, zf, Nz)*1e-3 #[m]

# MTQ Data Storage
i_data = np.zeros((Nz, N))
m_data = np.zeros((Nz, N))
B_data = np.zeros((Nz, N))
m_calc = np.zeros((Nz, N))

# Run test
for j in range(0, Nz):
    print("Distance: ", z[j], "\t", j+1, "/", Nz)
    for i in range(0, N):
        # Process data
        mtq.main_routine(voltage[i], dt)
        
        # Data of interest
        i_data[j][i] = current_sensor.measure(mtq.i)
        m_data[j][i] = mtq.m
        B_data[j][i] = magnetometer.measure_rod(z[j]+dz0, mtq.i, d/2, N_coil, L_rod, A_rod)
        m_calc[j][i] = (10)*B_data[j][i]*magnetometer.getMagneticMoment(z[j]+dz0, mtq.i, d/2, N_coil, L_rod, A_rod)
    mtq.reset()

# Save data to csv
SAVE_DATA = True
if SAVE_DATA:
    from pathlib import Path
    import datetime
    folder = "data/gst-600/electric/sim/"
    Path(folder).mkdir(parents=True, exist_ok=True)
    data = {"time[s]": t, "voltage[V]": voltage, "current[A]": i_data[0]}
    df = pd.DataFrame(data, columns=["time[s]", "voltage[V]", "current[A]"])
    test_name = "-test_distances[voltage-current]"
    date = datetime.datetime.now().strftime('%Y-%m-%d %H-%M-%S')
    path = Path(folder+date+test_name+".csv")
    df.to_csv(path)

    folder = "data/gst-600/magnetic/sim/"+date+"/"
    Path(folder).mkdir(parents=True, exist_ok=True)
    for j in range(0, Nz):
        data = {"time[s]": t, "magField_z[uT]": B_data[j]}
        df = pd.DataFrame(data, columns=["time[s]", "magField_z[uT]"])
        test_name = "-test_distances[magField]["+str(int(z[j]*1e3))+"mm]"
        path = Path(folder+test_name+".csv")
        df.to_csv(path)

# Data Visualization
from monitor import Monitor

v_mon = Monitor([t], [voltage], "MTQ input voltage", "V[V]", "time[s]", sig_name = ["V"])
v_mon.plot()

i_mon = Monitor([t], [i_data[0]], "MTQ electric current", "i[A]", "time[s]", sig_name = ["i"])
i_mon.plot()

B_mon = Monitor([t, t, t, t, t, t], [B_data[0], B_data[1], B_data[2], B_data[3], B_data[9], B_data[10]], "MTQ magnetic field", "B[uT]", "time[s]", sig_name = ["B0", "B1", "B2", "B3", "B9", "B10"])
B_mon.plot()

mCalc_mon = Monitor([t, t, t, t], [m_calc[0], m_calc[1], m_calc[2], m_calc[3]], "MTQ calculated magnetic moment", "m[Am2]", "time[s]", sig_name = ["m0", "m1", "m2", "m3"])
mCalc_mon.plot()

m_mon = Monitor([t, t, t, t], [m_data[0], m_data[1], m_data[2], m_data[3]], "MTQ magnetic moment", "m[Am2]", "time[s]", sig_name = ["m0", "m1", "m2", "m3"])
m_mon.plot()
m_mon.show()