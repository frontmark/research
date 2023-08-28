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
