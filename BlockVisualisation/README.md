# BlockVisualisation

Two fundamental graph visualizing ideas are implemented in this code.

## Kamada-Kawai Algorithm

https://www.sciencedirect.com/science/article/abs/pii/0020019089901026

An optimal euclidean distance between two nodes is calculated based on the corresponding graph distance. Then, an energy function based on the variance of node distances and the optimal distances is optimized with Newton's method.

## Clustering

We also implemented a simple clustering and unclustering algorithm, that allows Kamada-Kawai to be run on bigger graphs, but this is not currently used. Instead, we just use the ForceAtlas2 algorithm on big graphs.

Similar to https://github.com/visjs/vis-network

## ForceAtlas2

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679

A force-based layout algortihm.

## Fast Multipole

The forces are calculated with the fast multipole algorithm. Also includes an own quadtree implementation.

https://cpsc.yale.edu/sites/default/files/files/tr533.pdf

## Putting it together

Layouts for components are computed separately. For smaller components we use Kamada-Kawai and for bigger components ForceAtlas2. Then we construct a new graph consisting of the components as vertices, and calculate a Kamada-Kawai layout on the result, to place the components.

## Execution

### build

To build the executable, run:

```sh
g++ -Ofast main.cpp forces.cpp quadtree.cpp layout.cpp block_visual.cpp cluster.cpp -o calc_layout.exe
```

or on macOS with Homebrew's `gcc@13`:

```sh
/opt/homebrew/bin/g++-13 -Ofast main.cpp forces.cpp quadtree.cpp layout.cpp block_visual.cpp cluster.cpp -o calc_layout.exe
```

or use Docker's official `gcc:13` image to build the executable:

```sh
docker run -it --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp gcc:13 g++ -Ofast main.cpp forces.cpp quadtree.cpp layout.cpp block_visual.cpp cluster.cpp -o calc_layout.exe
```

### run

Run

```sh
./calc_layout.exe
```

on the command line to calculate a layout

or use Docker's official `gcc:13` image to do so:

```sh
docker run -it --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp gcc:13 ./calc_layout.exe
```

You should see something like this:

```sh
Reading file sample/2011_05_19.gml.
Calculating fa2 layout for graph with 14867 vertices:
...

Ran through in: 33.7754 seconds...
```

You should see a new file `sample/2011_05_19_out.txt` being created.

Adjust `sample/visual_settings.txt` to define the input and output files.

### draw

Run

```sh
make -C sample reinit
```

on the command line

- to build a Docker image with all necessary Python dependencies,
- to run the Docker image to draw the calculated layout.

Alternatively, execute `sample/gen_img_from_lay.py` directly.

You should see a new file `sample/2011_05_19_pic.png` being created.

Adjust `sample/gen_img_from_lay.py` to define the input and output files.
