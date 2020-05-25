import pandas
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from monitor import Monitor

# read data
path = "data/magnetometer/"
test_date = "2020-05-24 20-15-36/"

z0 = 16			   #[mm]
dz = 1             #[mm]
zf = 35            #[mm]
Nz = int((zf-z0+1)/dz) #[number of points]
z = np.linspace(z0, zf, Nz)*1e-3 #[m]

A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m]

data = []
for i in range(0, Nz):
	file = Path(path+test_date+"-test_distances[magField]["+str(int(z[i]*1e3))+"mm].csv")
	data.append(pandas.read_csv(file))

# process data
def find_step(t, step_response):
	n = len(step_response)
	delta = np.zeros(n)
	for i in range(0, n-1):
		delta[i] = step_response[i+1] - step_response[i]
	step_size = np.max(delta)
	k = np.argmax(delta) + 50
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

time = []
magF = []

magF_rel = []
t_magF_rel = []
magF_mean = []
magF_std = []

magM_calc = []
magM_calc_mean = []
magM_calc_std = []
for i in range(0, Nz):
	time.append(data[i]['time[s]'].values)
	magF.append(data[i]["magField[uT]"].values)
	[m_offset, m_steady, t_steady] = find_step(time[i], magF[i])
	magF_rel.append(m_steady-np.mean(m_offset))
	t_magF_rel.append(t_steady)
	magF_mean.append(np.mean(magF_rel[i]))
	magF_std.append(np.std(magF_rel[i]))

	magM_calc.append((5)*magF_rel[i]*(z[i]**2+(0.5*l)**2)*np.sqrt(z[i]**2+0.5*l**2))
	magM_calc_mean.append(np.mean(magM_calc[i]))
	magM_calc_std.append(np.std(magM_calc[i]))

magF_mean_mean = np.mean(magF_mean)
magF_std_std = np.std(magF_mean)
print("Magnetic Field")
print("mean [uT]: ", magF_mean_mean)
print("std [uT]: ", magF_std_std)

magM_mean_mean = np.mean(magM_calc_mean)
magM_std_std = np.std(magM_calc_mean)
print("Magnetic Moment")
print("mean [Am2]: ", magM_mean_mean)
print("std [Am2]: ", magM_std_std)

# Data Visualization
for i in range(0, 2):
	mF_mon = Monitor([time[i]], [magF[i]], "MTQ Magnetic Field raw", "B[uT]", "time[s]", sig_name = ["B"+str(i)])
	mF_mon.plot()

for i in range(0, 2):
	mF_r_mon = Monitor([t_magF_rel[i]], [magF_rel[i]], "MTQ Magnetic Field relative", "B[uT]", "time[s]", sig_name = ["B"+str(i)])
	mF_r_mon.plot()

for i in range(0, 2):
	mF_r_mon = Monitor([t_magF_rel[i]], [magM_calc[i]], "MTQ Magnetic Moment calculated", "m[Am2]", "time[s]", sig_name = ["m"+str(i)])
	mF_r_mon.plot()

magF_mean = np.array(magF_mean)
magF_std = np.array(magF_std)

magM_calc_mean = np.array(magM_calc_mean)
magM_calc_std = np.array(magM_calc_std)

mF_m = Monitor([z, z, z], [magF_mean-magF_std, magF_mean, magF_mean+magF_std], "MTQ Magnetic Field means", "B[uT]", "distance[mm]", marker = "*.*", sig_name = ["stdm", "mean", "stdp"])
mF_m.plot()

mM_c = Monitor([z, z, z], [magM_calc_mean-magM_calc_std, magM_calc_mean, magM_calc_mean+magM_calc_std], "MTQ Magnetic Moment means", "m[Am2]", "distance[mm]", marker = "*.*", sig_name = ["stdm", "mean", "stdp"])
mM_c.plot()
mM_c.show()