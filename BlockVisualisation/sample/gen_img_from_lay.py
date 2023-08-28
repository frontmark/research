import matplotlib as mpl
import matplotlib.collections  # call as mpl.collections
import matplotlib.pyplot as plt
import networkx as nx

INFILENAME = "2011_05_19.gml"
OUTFILENAME = "2011_05_19_out.txt"


edge_pos = []
xcoord = []
ycoord = []
vertex_type = []
with open(OUTFILENAME) as f:
    for line in f.readlines():
        split = line.split("|")
        xcoord.append(float(split[0]))
        ycoord.append(float(split[1]))


if "txt" in INFILENAME:
    with open(INFILENAME) as f:
        for line in f.readlines()[1:]:
            split = line.split(" ")
            edge_pos.append(
                (
                    (xcoord[int(split[0])], ycoord[int(split[0])]),
                    (xcoord[int(split[1])], ycoord[int(split[1])]),
                )
            )
elif "gml" in INFILENAME:
    g = nx.read_gml(INFILENAME, label="id")
    for edge in g.edges():
        edge_pos.append(
            ((xcoord[edge[0]], ycoord[edge[0]]), (xcoord[edge[1]], ycoord[edge[1]]))
        )

    with open(INFILENAME) as f:
        for line in f.readlines()[1:]:
            if "address" in line:
                vertex_type.append("b")
            if "timestamp" in line:
                vertex_type.append("r")


ax = plt.gca()


plt.scatter(
    xcoord,
    ycoord,
    s=0.3,
    c=vertex_type[: len(xcoord)],
    marker=".",
)


edge_collection = mpl.collections.LineCollection(
    edge_pos,
    colors="k",
    linewidths=0.1,
    antialiaseds=(1,),
    linestyle="solid",
    alpha=None,
)

ax.add_collection(edge_collection)
plt.axis("off")
plt.savefig("2011_05_19_pic.png", dpi=1000)
plt.show()
