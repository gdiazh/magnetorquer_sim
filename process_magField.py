import pandas
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from monitor import Monitor

# read data
path = "data/magnetometer/"
# test_date = "2020-05-25 13-06-49/"
test_date = "2020-05-27 18-43-10/"

z0 = 16			   #[mm]
dz = 1             #[mm]
zf = 35            #[mm]
Nz = int((zf-z0+1)/dz) #[number of points]
z = np.linspace(z0, zf, Nz)*1e-3 #[m]

A = 46.82e-4	#[m2]
l = np.sqrt(A)	#[m]

debug = False

data = []
for i in range(0, Nz):
	file = Path(path+test_date+"-test_distances[magField]["+str(int(z[i]*1e3))+"mm].csv")
	data.append(pandas.read_csv(file))

# correction sign
d0 = data[0]["magField_z[uT]"].values
if max(d0)>abs(min(d0)):
	correction_sign = 1
else:
	correction_sign = -1

# process data
def find_step(t, step_response):
	n = len(step_response)
	delta = np.zeros(n)
	for i in range(0, n-2):
		delta[i] = step_response[i+2] - step_response[i]
	step_size = np.max(delta)
	k = np.argmax(delta) + 50
	k2 = np.argmin(delta) - 50
	offset = []
	t_offset = []
	steady = []
	t_steady = []
	if debug:
		print("delta=", delta)
		print("np.argmax(delta)=", np.argmax(delta))
		print("np.argmin(delta)=", np.argmin(delta))
		print("delta(argmax)=", delta[np.argmax(delta)])
		print("delta(argmin)=", delta[np.argmin(delta)])
		print("t(argmin)=", t[np.argmin(delta)])
		print("k=", k)
		print("k2=", k2)
	if k == 0:
		stp = Monitor([t, t], [step_response, delta], "Step response", "step[]", "time[s]", sig_name = ["raw", "steady"])
		stp.plot()
		stp.show()
		raise Exception('No Step Found')
	elif k2>k and min(delta)<-30:
		offset = step_response[0:k-50]
		t_offset = t[0:k-50]
		steady = step_response[k:k2]
		t_steady = t[k:k2]
		if debug:
			stp = Monitor([t, t_steady, t_offset], [step_response, steady, offset], "Step response", "step[]", "time[s]", sig_name = ["raw", "steady", "offset"])
			stp.plot()
			stp.show()
	else:
		offset = step_response[0:k-50]
		t_offset = t[0:k-50]
		steady = step_response[k:n]
		t_steady = t[k:n]
		if debug:
			stp = Monitor([t, t_steady, t_offset], [step_response, steady, offset], "Step response", "step[]", "time[s]", sig_name = ["raw", "steady", "offset"])
			stp.plot()
			stp.show()
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
	time.append(1e-3*data[i]['time[s]'].values)
	magF.append(correction_sign*data[i]["magField_z[uT]"].values)
	[m_offset, m_steady, t_steady] = find_step(time[i], magF[i])
	magF_rel.append(m_steady-np.mean(m_offset))
	t_magF_rel.append(t_steady)
	magF_mean.append(np.mean(magF_rel[i]))
	magF_std.append(np.std(magF_rel[i]))
	magM_calc.append((5)*magF_rel[i]*(z[i]**2+(0.5*l)**2)*np.sqrt(z[i]**2+0.5*l**2))
	magM_calc_mean.append(np.mean(magM_calc[i]))
	magM_calc_std.append(np.std(magM_calc[i]))

magF_mean = np.array(magF_mean)
magF_std = np.array(magF_std)
magM_calc_mean = np.array(magM_calc_mean)
magM_calc_std = np.array(magM_calc_std)

magM_calc_mean2 = (5)*np.multiply(magF_mean,np.multiply(z**2+(0.5*l)**2,np.sqrt(z**2+0.5*l**2)))
magM_calc_var2 = (5**2)*np.dot(magF_std**2,np.multiply(z**2+(0.5*l)**2,np.sqrt(z**2+0.5*l**2))**2)
magM_calc_std2 = np.sqrt(magM_calc_var2)
print("magF_mean: ", magF_mean)
print("magM_calc_mean: ", magM_calc_mean)
print("magM_calc_mean2: ", magM_calc_mean2)

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

# Save Main Results
data = {"z[m]": z, "B_mean[uT]": magF_mean, "m_mean[Am2]": magM_calc_mean, "B_std": magF_std, "m_std": magM_calc_std}
df = pandas.DataFrame(data, columns=["z[m]", "B_mean[uT]", "m_mean[Am2]", "B_std", "m_std"])
file = Path(path+test_date+"results[B-m-mean-std][16-35mm].csv")
df.to_csv(file)

# Data Visualization

mF_m = Monitor([z*1e3], [magF_mean, magF_std], "MTQ Magnetic Field means", "B[uT]", "distance[mm]", marker = "o-", sig_name = ["B_mean"])
mF_m.errorbar()

mM_c = Monitor([z*1e3], [magM_calc_mean*1e3, magM_calc_std*1e3], "MTQ Magnetic Moment means", "m[mAm2]", "distance[mm]", marker = "o-", sig_name = ["m_mean"])
mM_c.errorbar()

mM2_c = Monitor([z*1e3], [magM_calc_mean2*1e3, magM_calc_std2*1e3], "MTQ Magnetic Moment means", "m[mAm2]", "distance[mm]", marker = "o-", sig_name = ["m_mean"])
mM2_c.errorbar()
mM2_c.show()