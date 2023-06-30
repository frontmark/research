from collections.abc import Iterable

import matplotlib as mpl
import matplotlib.collections  # call as mpl.collections
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx



INFILENAME = "sample/2011_04_19.gml"
OUTFILENAME = "sample/2011_04_19_out.txt"



edge_pos = []
xcoord = []
ycoord = []
with open(OUTFILENAME) as f:
    for line in f.readlines():
        #print(line)
        split = line.split('|')
        xcoord.append(float(split[0]))
        ycoord.append(float(split[1]))



if "txt" in INFILENAME: 
    with open(INFILENAME) as f:
        for line in f.readlines()[1:]:
            split = line.split(' ')
            edge_pos.append(((xcoord[int(split[0])], ycoord[int(split[0])]), (xcoord[int(split[1])], ycoord[int(split[1])])))
elif "gml" in INFILENAME:
    g = nx.read_gml(INFILENAME, label='id')
    for edge in g.edges():
        edge_pos.append(((xcoord[edge[0]], ycoord[edge[0]]), (xcoord[edge[1]], ycoord[edge[1]])))
    #print(len(g.nodes()))
    #print(nx.number_weakly_connected_components(g))

ax = plt.gca()

#print(xcoord)
#print(ycoord)

node_collection = ax.scatter(
    xcoord,
    ycoord,
    s=30,
    c="#1f78b4",
    marker="o",
)


edge_collection = mpl.collections.LineCollection(
    edge_pos,
    colors="k",
    linewidths=1.0,
    antialiaseds=(1,),
    linestyle="solid",
    alpha=None,
)

ax.add_collection(edge_collection)

plt.show()