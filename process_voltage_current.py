import pandas
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# read data
path = "data/osc/"
file = "2020-05-24 20-15-36-test_distances[voltage-current].csv"
file_ = Path(path+file)
data = pandas.read_csv(file_)

time = data['time[s]'].values
voltage = data["voltage[V]"].values
current = data["current[A]"].values

# Data Visualization
from monitor import Monitor

v_mon = Monitor([time], [voltage], "MTQ input voltage", "V[V]", "time[s]", sig_name = ["V"])
v_mon.plot()

i_mon = Monitor([time], [current], "MTQ electric current", "i[A]", "time[s]", sig_name = ["i"])
i_mon.plot()
i_mon.show()