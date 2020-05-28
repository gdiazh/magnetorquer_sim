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

B_mon = Monitor([z_sim, z_exp], [B_mean_sim, B_mean_exp], "MTQ magnetic field", "B[uT]", "distance[mm]", marker="**", sig_name = ["B_sim", "B_exp"])
B_mon.plot()

m_mon = Monitor([z_sim, z_exp], [m_mean_sim, m_mean_exp], "MTQ magnetic moment", "m[mAm2]", "distance[mm]", marker="**", sig_name = ["m_sim", "m_exp"])
m_mon.plot()

m_mon.show()