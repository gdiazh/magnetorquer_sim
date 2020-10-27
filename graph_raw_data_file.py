import glob
import pandas
import numpy as np
import sys
sys.path.append("..")
from Visualization.monitor import Monitor

# read data
path = "../../DRV10987_Firmware/ros_uart_controller/data/inisteps/"
name = "2020-10-13 21-51-51[511]"
file = path+name+".csv"
data = pandas.read_csv(file)

time = data['time[s]'].values
voltage = data['voltage[V]'].values
current = data['current[mA]'].values

# Graph data
cmd_mon = Monitor([time], [voltage], "MTQ-GST600 Voltage", "V[V]", "time[s]", sig_name = ["V"])
cmd_mon.plot()

speed_mon = Monitor([time], [current], "MTQ-GST600 Current", "i[mA]", "time[s]", sig_name = ["i"])
speed_mon.plot()
speed_mon.show()