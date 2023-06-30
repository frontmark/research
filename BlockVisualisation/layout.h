// ffm.h
#ifndef LAYOUT_H
#define LAYOUT_H

#include "forces.h"
#include "graph.h"

#include <vector>
#include <tuple>
#include <complex>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <bitset>
#include <string>

#include <iostream>
#include <fstream>

class LayoutGraph
{

public:
    LayoutGraph();
    void calc_fa2(Graph& graph, int iterations, int _terms, int _thresh, double _edge_force, double _gravity, double _jitter_tol, double _speed, double _speed_efficiency, int verbosity);
    void calc_kamada_kawai(Graph& graph, double edge_strength, double thresh, int maxIter, int verbosity);
    std::vector<std::vector<Compl>> calc_components_layout(std::vector<std::vector<Compl>> positions, int iterations=5);



private:
};

#endif