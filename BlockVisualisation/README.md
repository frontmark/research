# GraphVisualisation

Two fundamental graph visualizing ideas are implemented in this code.
## Kamada-Kawai Algorithm

https://www.sciencedirect.com/science/article/abs/pii/0020019089901026

Here an optimal euclidean distance between two nodes is calculated based on the corresponding graph distance.
Then an energy function based on the variance of node distances and the optimal distances is optimized with newton's method.

## Clustering
We also implemented a simple clustering and unclustering algorithm, that allows Kamada Kawai to be run on bigger graphs, but this is not currently used, instead we just use the forceatlas2 algorithm on big graphs.

Similar to https://github.com/visjs/vis-network

## ForceAtlas2

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679

An force-based layout algortihm.

## Fast Multipole
The forces are calculated with the fast multipole algorithm. Also includes an own quadtree implementation.

https://cpsc.yale.edu/sites/default/files/files/tr533.pdf

## Putting it Together
Layouts for components are computed separately. For smaller components we use Kamada-Kawai and for bigger components ForceAtlas2. Then we construct a new graph consisting of the components as vertices, and calculate a kamada kawai layout on the result, to place the components.

# Execution
## build

g++ -Ofast main.cpp forces.cpp quadtree.cpp layout.cpp block_visual.cpp cluster.cpp -o calc_layout.exe

## run
adjust visual_settings.txt
calc_layout.exe will now calculate a layout

## draw
adjust filenames in, and run gen_img_from_lay.py to draw the calculated layout
