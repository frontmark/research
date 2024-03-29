// ffm.h
#ifndef LAYOUT_H
#define LAYOUT_H

#include "forces.h"
#include "graph.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <complex>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <fstream>
#include <iostream>

class LayoutGraph {

public:
  LayoutGraph();
  void calc_fa2(Graph &graph, int iterations, int _terms, int _thresh,
                double _edge_force, double _gravity, double _jitter_tol,
                double _speed, double _speed_efficiency, int verbosity);
  void calc_kamada_kawai(Graph &graph, double edge_strength, double thresh,
                         int maxIter, int verbosity);
  std::vector<std::vector<Compl>>
  calc_components_layout(std::vector<std::vector<Compl>> positions,
                         int iterations = 5);
  std::vector<std::vector<Compl>>
  calc_many_components_layout(std::vector<std::vector<Compl>> positions,
                              int iterations = 5, int comp_per_calc = 500);

private:
};

#endif