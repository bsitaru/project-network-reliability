import matplotlib.pyplot as plt
import matplotlib.cbook as cbook

import csv
import numpy as np
import pandas as pd

file_object = open('plot_data/complete.dat', "r")
file_data = pd.read_csv('plot_data/complete.dat', delimiter=' ')

# print(file_data)

file_data.plot("#p", ["unrel_c10", "unrel_c20", "unrel_c30", "unrel_g3", "unrel_g5", "unrel_g7"], linestyle='-')
plt.show()