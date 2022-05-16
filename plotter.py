import matplotlib.pyplot as plt
import matplotlib.cbook as cbook

import json
import csv
import numpy as np
import pandas as pd

import seaborn as sb

from collections import defaultdict
from itertools import chain

with open('experiment_results/empirical_eps_rel.json', 'r') as f:
    j = json.load(f)

# ps = ['0.001000', '0.010000', '0.100000', '0.200000', '0.250000', '0.500000', '0.750000', '0.900000']
ps = ['0.100000', '0.250000', '0.500000', '0.750000', '0.900000']
epses = ['0.050000', '0.100000', '0.200000']

def get_act_eps(e, ps):

    l = []
    for a in j:
        for b in ps:
            l.append( abs(j[a][b][e]['act_eps']) )
    return l

def plot_eps(e):
    l = get_act_eps(e, ps)
    sb.displot(l, kind='kde')

# plot_eps('0.050000')
# plot_eps('0.100000')
# plot_eps('0.200000')
# plt.show()

def get_empirical_eps(e, delta):
    l = sorted( get_act_eps(e, ps) )
    pos = int(len(l) * (1.0 - delta))
    return l[pos]

def get_empirical_delta(e, eps):
    l = sorted( get_act_eps(e, ps) )
    cnt = 0
    for x in l:
        if x <= eps:
            cnt += 1
    act_delta = 1.0 - float(cnt) / float(len(l))
    return act_delta

# for e in epses:
#     print( get_empirical_eps(e, 0.05) )

print( get_empirical_delta('0.050000', 0.05) )
print( get_empirical_delta('0.100000', 0.1) )
print( get_empirical_delta('0.200000', 0.2) )


# data = {
#     # '0.25': get_act_eps('0.250000', ps),
#     # '0.20': get_act_eps('0.200000', ps),
#     '0.10': get_act_eps('0.100000', ps),
#     # '0.05': get_act_eps('0.050000', ps)
#
# }
# df = pd.DataFrame(data)
# # print(l)
#
# sb.displot(data=df, kind='kde')
#
# plt.show()

# l = [y["0.010000"]["act_eps"] for y in chain(x for x in j.values())]
# print(l)

# print(list(chain([t.values() for t in j.values()])))

# file_object = open('plot_data/complete.dat', "r")
# file_data = pd.read_csv('plot_data/complete.dat', delimiter=' ')

# print(file_data)

# file_data.plot("#p", ["unrel_c10", "unrel_c20", "unrel_c30", "unrel_g3", "unrel_g5", "unrel_g7"], linestyle='-')
# plt.show()

# with open('experiment_results/complete_diff_p_exp.json', 'r') as f:
# with open('experiment_results/grids_diff_p_exp.json', 'r') as f:
#     j = json.load(f)
# j = pd.read_json('experiment_results/complete_diff_p_exp.json')
# data = defaultdict(lambda: dict())
# print(j)
# for graph_id in j:
# n = int(graph_id.removeprefix("grid_"))
# n = int(graph_id.removeprefix("complete_"))
# n = n * n
# data['time'][n] = j[graph_id]['average_time'] / 1000000.0
# print(graph_id)
# for subgraph_id in j[graph_id]:
# if subgraph_id.startswith(graph_id):
# data[graph_id][j[graph_id][subgraph_id]['p']] = j[graph_id][subgraph_id]['try 1']['answer']
# df = pd.DataFrame(data)
# df.plot()
# plt.show()
# print(df)
