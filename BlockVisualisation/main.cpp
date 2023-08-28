// square.cpp (Compute the Square of an Integer)

#include "block_visual.h"
#include "forces.h"
#include "layout.h"
#include "quadtree.h"
#include <iostream>

#include <chrono>
#include <complex>
#include <iostream>
#include <vector>

using namespace std::chrono;

// hardcoded that terms can be at most 16 since we used bitsets for the binomial
// maps g++ main.cpp forces.cpp quadtree.cpp g++ main.cpp forces.cpp quadtree.cpp
// layout.cpp

/*
Reading file sample/2019_06_15.gml.
Calculating layout for graph with 1297265 vertices:
100%Ran through in: 2714.97 seconds...

100 iterations*/

int main() {
  auto start = high_resolution_clock::now();

  BlockVisual BV = BlockVisual("sample/visual_settings.txt");
  BV.calc_visual();
  // BV.fa2_visual();

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout << "\nRan through in: " << double(duration.count()) / 1000000.0
            << " seconds...\n";
}
