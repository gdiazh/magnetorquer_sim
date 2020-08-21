import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtq import MTQ 
from step import Step
from pwm import PWM
from current_sensor import CurrentSensor
from magnetometer import Magnetometer

# MTQ Model
R = 145.9 	#[Ohm]
L = 10.08e-3	#[H]
A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m] 
N_coil = 25*10	#[number of turns]
tau = L/R    	#[s]
mtq = MTQ(R, L, A, N_coil)

# Time parameters for simulation
tf = 1      
dt = 0.2*tau ##VER BIEN EL TIMEPO CON LA FREC MAX Y EL TAU (DELTA T PARECIDO A 10 MICRO SEUGUNDOS)) 
N = int(tf/dt)
t = np.linspace(0, tf, N)

# Voltage Signal parameters
V = 3.3 #[V]
duty=50 #[%]
#PWM input              #PWM nombre que se le da a la señal  o       
f0 = 1                  #[Hz] valor inicial #CAMBIE fo=0 a fo=1 para que no me tire "RuntimeWarning: divide by zero encountered in double_scalars"
fd = 100                #[Hz]pasos 
ff = 10000              #[Hz] valor final
Nf = int((ff-f0)/fd+1)  #[number of points] numero de valores que voy a tener
freq = np.linspace(f0, ff, Nf) #[%] 
voltage = np.zeros((Nf, N)) #matriz tabla de valores dependiendo de frecuencia y tiempo)
v_rms = np.zeros(Nf)
v_med = np.zeros(Nf)


for i in range (0,Nf):
    voltage[i,:] = PWM(0, tf, dt, V, freq[i], duty).signal #[i,:], filas
    v_rms[i] = np.sqrt(np.mean(voltage[i,:]**2)) #rms root mean sqr value
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


# MTQ Data Storage #define variables para guardar datos dsp, se inician todos los datos en cero
i_data = np.zeros((Nf, N))
m_data = np.zeros((Nf, N))
B_data = np.zeros((Nf, N))
m_calc = np.zeros((Nf, N))

i_rms = np.zeros(Nf)
m_rms = np.zeros(Nf)
B_rms = np.zeros(Nf)
m_calc_rms = np.zeros(Nf)

i_med = np.zeros(Nf)
m_med = np.zeros(Nf)
B_med = np.zeros(Nf)
m_calc_med = np.zeros(Nf)


#Run Test
for j in range(0, Nf):
    print("PWM frequency\t: ", freq[j], "\t", j+1, "/", Nf) 
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


# Data Visualization
from monitor import Monitor
v_mon= Monitor([t, t, t, t], [voltage[0], voltage[30], voltage[60], voltage[99]],"MTQ input freq", "V[V]","time[s]", sig_name = ["V0", "V1", "V2", "V3"])
v_mon.plot()

i_mon = Monitor([t, t, t, t], [i_data[0], i_data[30], i_data[60], i_data[99]], "MTQ electric current", "i[A]", "time[s]", sig_name = ["i0", "i1", "i2", "i3"])
i_mon.plot()

B_mon = Monitor([t, t, t, t], [B_data[0], B_data[30], B_data[60], B_data[99]], "MTQ magnetic field", "B[uT]", "time[s]", sig_name = ["B0", "B1", "B2", "B3"])
B_mon.plot()

mCalc_mon = Monitor([t, t, t, t], [m_calc[0], m_calc[30], m_calc[60], m_calc[99]], "MTQ calculated magnetic moment", "m[Am2]", "time[s]", sig_name = ["m0", "m1", "m2", "m3"])
mCalc_mon.plot()

m_mon = Monitor([t, t, t, t], [m_data[0], m_data[30], m_data[60], m_data[99]], "MTQ magnetic moment", "m[Am2]", "time[s]", sig_name = ["m0", "m1", "m2", "m3"])
m_mon.plot()

vrms_mon = Monitor([freq], [v_rms], "MTQ input rms voltage", "V[Vrms]", "freq[Hz]", sig_name = ["Vrms"])
vrms_mon.plot()

vmed_mon = Monitor([freq], [v_med], "MTQ input med voltage", "V[Vmed]", "freq[Hz]", sig_name = ["Vmed"])
vmed_mon.plot()

irms_mon = Monitor([freq], [i_rms], "MTQ input rms current", "i[irms]", "freq[Hz]", sig_name = ["irms"])
irms_mon.plot()

Brms_mon = Monitor([freq], [B_rms], "MTQ magnetic field rms", "B[Brms]", "freq[Hz]", sig_name = ["Brms"])
Brms_mon.plot() 

mCalcrms_mon = Monitor([freq], [m_calc_rms], "MTQ calc magnetic moment rms", "mCalc[mrms]", "time[s]", sig_name = ["mrms"])
mCalcrms_mon.plot()

mrms_mon = Monitor([freq], [m_rms], "MTQ magnetic moment rms", "m[mrms]", "time[s]", sig_name = ["mrms"])
mrms_mon.plot()

imed_mon = Monitor([freq], [i_med], "MTQ input med current", "i[imed]", "time[s]", sig_name = ["imed"])
imed_mon.plot()

Bmed_mon = Monitor([freq], [B_med], "MTQ magnetic field med", "B[Bmed]", "time[s]", sig_name = ["Bmed"])
Bmed_mon.plot()

mCalcmed_mon = Monitor([freq], [m_calc_med], "MTQ calc magnetic moment med", "mCalc[mmed]", "time[s]", sig_name = ["mmed"])
mCalcmed_mon.plot()

mmed_mon = Monitor([freq], [m_med], "MTQ magnetic moment med", "m[mmed]", "time[s]", sig_name = ["mmed"])
mmed_mon.plot()
v_mon.show()














