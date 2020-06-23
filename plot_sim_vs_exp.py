import pandas
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from monitor import Monitor

# read data
path = "data/magnetometer/"
test_exp = "2020-05-25 13-06-49/"
test_sim = "2020-05-27 18-43-10/"

file_exp = Path(path+test_exp+"results[B-m-mean-std][16-35mm].csv")
file_sim = Path(path+test_sim+"results[B-m-mean-std][16-35mm].csv")

data_exp = pandas.read_csv(file_exp)
data_sim = pandas.read_csv(file_sim)

z_exp = data_exp['z[m]'].values*1e3
B_mean_exp = data_exp['B_mean[uT]'].values
B_std_exp = data_exp['B_std'].values
m_mean_exp = data_exp['m_mean[Am2]'].values*1e3
m_std_exp = data_exp['m_std'].values

z_sim = data_sim['z[m]'].values*1e3
B_mean_sim = data_sim['B_mean[uT]'].values
B_std_sim = data_sim['B_std'].values
m_mean_sim = data_sim['m_mean[Am2]'].values*1e3
m_std_sim = data_sim['m_std'].values

print(B_mean_exp)
print(B_std_sim)

# Calc Least Square Regretion
n = np.sum((z_exp - np.mean(z_exp)) * (B_mean_exp - np.mean(B_mean_exp)))
d = np.sum((z_exp - np.mean(z_exp))**2)
m = n / d
b = np.mean(B_mean_exp) - m * np.mean(z_exp)
print("Experiment Coefficients:")
print("m = ", m, "b = ", b)
B_lsr = m * z_exp + b
Brmse_exp = np.sqrt(np.mean((B_lsr - B_mean_exp)**2))
print("B RMSE Experiment: ",Brmse_exp, "[uT]")
A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m]
M_est = np.abs(m) * 5 * np.sqrt(2)*(l/2)**3
# print("M estimated from Experiment LSR: ", M_est*1e3, "[mAm2]")

n2 = np.sum((z_sim - np.mean(z_sim)) * (B_mean_sim - np.mean(B_mean_sim)))
d2 = np.sum((z_sim - np.mean(z_sim))**2)
m2 = n2 / d2
b2 = np.mean(B_mean_sim) - m2 * np.mean(z_sim)
print("Simulation Coefficients:")
print("m = ", m2, "b = ", b2)
B_lsr2 = m2 * z_sim + b2
Brmse_sim = np.sqrt(np.mean((B_lsr2 - B_mean_sim)**2))
print("B RMSE Simulation: ",Brmse_sim, "[uT]")
M_est2 = np.abs(m2) * 5 * np.sqrt(2)*(l/2)**3
# print("M estimated from Simulation LSR: ", M_est2*1e3, "[mAm2]")

l_mm = l*1e3
fz_exp = 1 / ((z_exp**2+(l_mm/2)**2)*np.sqrt(z_exp**2+2*(l_mm/2)**2))
n3 = np.sum((fz_exp - np.mean(fz_exp)) * (B_mean_exp - np.mean(B_mean_exp)))
d3 = np.sum((fz_exp - np.mean(fz_exp))**2)
m3 = n3 / d3
b3 = np.mean(B_mean_exp) - m3 * np.mean(fz_exp)
print("Experiment Coefficients:")
print("m3 = ", m3, "b = ", b3)
M_est3 = np.abs(m3) * 5 * 1e-9
print("M estimated from Simulation LSR: ", M_est3*1e3, "[mAm2]")

fz_sim = 1 / ((z_sim**2+(l_mm/2)**2)*np.sqrt(z_sim**2+2*(l_mm/2)**2))
n4 = np.sum((fz_sim - np.mean(fz_sim)) * (B_mean_sim - np.mean(B_mean_sim)))
d4 = np.sum((fz_sim - np.mean(fz_sim))**2)
m4 = n4 / d4
b4 = np.mean(B_mean_sim) - m4 * np.mean(fz_sim)
print("Experiment Coefficients:")
print("m4 = ", m4, "b = ", b4)
M_est4 = np.abs(m4) * 5 * 1e-9
print("M estimated from Simulation LSR: ", M_est4*1e3, "[mAm2]")

B_mon = Monitor([z_sim, z_exp, z_exp, z_sim], [B_mean_sim, B_mean_exp, B_lsr, B_lsr2], "MTQ magnetic field", "B[uT]", "distance[mm]", color = ["b", "orange", "black", "black"], marker="**..", sig_name = ["B_sim", "B_exp", "B_LSR_exp", "B_LSR_sim"])
B_mon.plot()

m_mon = Monitor([z_sim, z_exp], [m_mean_sim, m_mean_exp], "MTQ magnetic moment", "m[mAm2]", "distance[mm]", marker="**", sig_name = ["m_sim", "m_exp"])
m_mon.plot()

m_mon.show()