import numpy as np
import glob
import os
from pathlib import Path

path = "data/gst-600/2020-10-05 16-04-53_x/"
path_fmt = "data/gst-600/2020-10-05 16-04-53_x/fmt/"
files = glob.glob(path+"*.csv")

for file in files:
    print(file)
    with open(file, 'r') as f_in:
        with open(path_fmt+file[35:],'w+') as f_out:
            for line_no, line in enumerate(f_in, 1):
                if line_no == 1:
                    f_out.write(line[:-1]+',foo\n')
                else:
                    f_out.write(line)