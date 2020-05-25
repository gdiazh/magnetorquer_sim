import pandas
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from monitor import Monitor

# read data
path = "data/magnetometer/"
test_date = "2020-05-25 01-06-13/"
file = "-test_pwm[magField][10%].csv"

duty = 10	#[%]
mag_offset = 10	#[uT]

z = 16e-3 #[m]
A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m]

file_ = Path(path+test_date+file)
data = pandas.read_csv(file_)

time = data['time[s]'].values
B_m = data["magField[uT]"].values

Br = np.abs(B_m - mag_offset)
mc = (5)*Br*(z**2+(0.5*l)**2)*np.sqrt(z**2+0.5*l**2)

# rms values
Br_rms = np.sqrt(np.mean(Br**2))
mc_rms = np.sqrt(np.mean(mc**2))

# mean values
Br_med = np.mean(Br)
mc_med = np.mean(mc)

# std values
Br_std = np.std(Br)
mc_std = np.std(mc)

print("Magnetic Field")
print("mean [uT]: ", Br_med)
print("rms [uT]: ", Br_rms)
# print("std [uT]: ", Br_std)

print("Magnetic Moment")
print("mean [Am2]: ", mc_med)
print("rms [Am2]: ", mc_rms)
# print("std [Am2]: ", mc_std)

# Data Visualization
mF_mon = Monitor([time], [B_m], "MTQ Magnetic Field raw", "B[uT]", "time[s]", sig_name = ["B"])
mF_mon.plot()

mF_r_mon = Monitor([time], [Br], "MTQ Magnetic Field relative", "B[uT]", "time[s]", sig_name = ["B"])
mF_r_mon.plot()

mF_r_mon = Monitor([time], [mc], "MTQ Magnetic Moment calculated", "m[Am2]", "time[s]", sig_name = ["m"])
mF_r_mon.plot()
mF_r_mon.show()