// block_visual.h
#ifndef BLOCK_VISUAL_H
#define BLCOK_VISUAL_H

#include "forces.h"

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

struct Settings {

  int componentwise = 0;
  int iterations = 100;
  double edge_force = 0.04;
  double gravity = 0.04;
  double jitter_tol = 1.0;
  double speed = 10.0;
  double speed_efficiency = 1.0;
  int layout_verbosity = 1;
  int terms = 2;
  int thresh = -1;
  double kk_edge_strength = 1;
  double kk_thresh = 0.00001;
  int kk_iterations = 100;
  int kk_verbosity = 2;
  int max_kamada_kawai_nodes = 200;
};

class BlockVisual {

public:
  BlockVisual(std::string file_settings);
  void calc_visual();
  void fa2_visual();

private:
  std::string _file_path_in;
  std::string _file_path_out;
  int _timestamp;
  Settings _settings;
};

#endif