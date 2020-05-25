import pandas
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from monitor import Monitor

# read data
path = "data/magnetometer/"
test_date = "2020-05-24 20-15-36/"
file = "-test_distances[magField][17mm].csv"

z = 17e-3 #[m]
A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m]

file_ = Path(path+test_date+file)
data = pandas.read_csv(file_)

# process data
def find_step(t, step_response):
	n = len(step_response)
	delta = np.zeros(n)
	for i in range(0, n-1):
		delta[i] = step_response[i+1] - step_response[i]
	step_size = np.max(delta)
	dk = 50
	k = np.argmax(delta) + dk
	print("argmax_k: ", k, dk)
	offset = []
	steady = []
	t_steady = []
	if k == 0:
		stp = Monitor([t, t], [step_response, delta], "Step response", "step[]", "time[s]", sig_name = ["s"])
		stp.plot()
		stp.show()
		raise Exception('No Step Found')
	else:
		offset = step_response[0:k]
		steady = step_response[k:n]
		t_steady = t[k:n]
	return [offset, steady, t_steady]

time = data['time[s]'].values
magF = data["magField[uT]"].values

[m_offset, m_steady, t_steady] = find_step(time, magF)
magF_rel = m_steady-np.mean(m_offset)
t_magF_rel = t_steady
magF_mean = np.mean(magF_rel)
magF_std = np.std(magF_rel)

magM_calc = (5)*magF_rel*(z**2+(0.5*l)**2)*np.sqrt(z**2+0.5*l**2)
magM_calc_mean = np.mean(magM_calc)
magM_calc_std = np.std(magM_calc)

print("Magnetic Field")
print("mean [uT]: ", magF_mean)
print("std [uT]: ", magF_std)

print("Magnetic Moment")
print("mean [Am2]: ", magM_calc_mean)
print("std [Am2]: ", magM_calc_std)

# Data Visualization
mF_mon = Monitor([time], [magF], "MTQ Magnetic Field raw", "B[uT]", "time[s]", sig_name = ["B"])
mF_mon.plot()

mF_r_mon = Monitor([t_magF_rel], [magF_rel], "MTQ Magnetic Field relative", "B[uT]", "time[s]", sig_name = ["B"])
mF_r_mon.plot()

mF_r_mon = Monitor([t_magF_rel], [magM_calc], "MTQ Magnetic Moment calculated", "m[Am2]", "time[s]", sig_name = ["m"])
mF_r_mon.plot()
mF_r_mon.show()