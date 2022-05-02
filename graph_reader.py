import networkx as nx

g = nx.read_graphml('graphs/geant2009.graphml')

print(len(g.nodes()), len(g.edges()))

for u, v in g.edges():
    print(str(u) + " " + str(v))
