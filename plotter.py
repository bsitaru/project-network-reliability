import matplotlib.pyplot as plt
import matplotlib.cbook as cbook

import json
import csv
import numpy as np
import pandas as pd

from collections import defaultdict

# file_object = open('plot_data/complete.dat', "r")
# file_data = pd.read_csv('plot_data/complete.dat', delimiter=' ')

# print(file_data)

# file_data.plot("#p", ["unrel_c10", "unrel_c20", "unrel_c30", "unrel_g3", "unrel_g5", "unrel_g7"], linestyle='-')
# plt.show()

with open('experiment_results/complete_diff_p_exp.json', 'r') as f:
# with open('experiment_results/grids_diff_p_exp.json', 'r') as f:
    j = json.load(f)
# j = pd.read_json('experiment_results/complete_diff_p_exp.json')
data = defaultdict(lambda: dict())
# print(j)
for graph_id in j:
    # n = int(graph_id.removeprefix("grid_"))
    n = int(graph_id.removeprefix("complete_"))
    # n = n * n
    data['time'][n] = j[graph_id]['average_time'] / 1000000.0
    # print(graph_id)
    # for subgraph_id in j[graph_id]:
        # if subgraph_id.startswith(graph_id):
            # data[graph_id][j[graph_id][subgraph_id]['p']] = j[graph_id][subgraph_id]['try 1']['answer']
df = pd.DataFrame(data)
df.plot()
plt.show()
# print(df)
