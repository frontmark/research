#include "layout.h"
#include "quadtree.h"

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using Compl = std::complex<double>;

LayoutGraph::LayoutGraph(){};

bool is_normal(Compl a) {
  return std::isnormal(a.real()) && std::isnormal(a.imag());
}

void LayoutGraph::calc_kamada_kawai(Graph &graph, double edge_strength,
                                    double thresh, int maxIter, int verbosity) {

  if (graph.get_particles().size() == 1) {
    return;
  }

  // initial rescaling so average edge size is correct.
  double distance = 0;

  for (auto e : graph.get_edges()) {
    distance += std::abs(graph.get_particles()[e.source].get_pos() -
                         graph.get_particles()[e.target].get_pos());
  }
  double rescale = graph.get_edges().size() / (distance * 1.2);

  std::vector<Particle> &_particles = graph.get_particles();
  for (auto p : _particles) {
    p.mul_pos(rescale);
  }

  if (verbosity >= 2) {
    std::cout << "Calculating kamada kawai layout for graph with "
              << graph.num_vertices() << " vertices: \n";
  }

  // distances matrix
  std::vector<std::vector<double>> distances = graph.floyd_warschall();
  double biggest_distance = 0;
  for (int i = 0; i < graph.num_vertices(); i++) {
    int m = *std::max_element(distances[i].begin(), distances[i].end());
    if (m > biggest_distance) {
      biggest_distance = m;
    }
  }

  // "radius"=1
  double edge_target_length = 1.0;

  // length matrix, with edge length times shortest path
  std::vector<std::vector<double>> length_mat(
      graph.num_vertices(), std::vector<double>(graph.num_vertices()));
  for (int i = 0; i < graph.num_vertices(); i++) {
    for (int j = i; j < graph.num_vertices(); j++) {
      if (i == j) {
        length_mat[i][j] = 0.0;
      } else {
        length_mat[i][j] = distances[i][j] * edge_target_length;
        length_mat[j][i] = length_mat[i][j];
      }
    }
  }

  // spring matrix, with spring constants times shortest path
  std::vector<std::vector<double>> spring_mat(
      graph.num_vertices(), std::vector<double>(graph.num_vertices()));
  for (int i = 0; i < graph.num_vertices(); i++) {
    for (int j = i; j < graph.num_vertices(); j++) {
      if (i == j) {
        spring_mat[i][j] = 0.0;
      } else {
        spring_mat[i][j] = edge_strength / (distances[i][j] * distances[i][j]);
        spring_mat[j][i] = spring_mat[i][j];
      }
    }
  }

  // energy matrix, with energy between each pair of nodes. these are complex
  // numbers for energy from x direction and y direction
  std::vector<std::vector<std::complex<double>>> energy_mat(
      graph.num_vertices(),
      std::vector<std::complex<double>>(graph.num_vertices()));
  std::vector<std::complex<double>> energy_mat_sums(graph.num_vertices(),
                                                    std::complex<double>(0, 0));
  for (int i = 0; i < graph.num_vertices(); i++) {
    for (int j = i; j < graph.num_vertices(); j++) {
      if (i == j) {
        energy_mat[i][j] = 0.0;
      } else {
        std::complex<double> diff =
            _particles[i].get_pos() - _particles[j].get_pos();
        double E_dx = spring_mat[i][j] * (diff.real()) *
                      (1 - length_mat[i][j] / (std::abs(diff)));
        double E_dy = spring_mat[i][j] * (diff.imag()) *
                      (1 - length_mat[i][j] / (std::abs(diff)));
        energy_mat[i][j] = std::complex<double>(E_dx, E_dy);
        energy_mat[j][i] = std::complex<double>(-E_dx, -E_dy);
        energy_mat_sums[i] += energy_mat[i][j];
        energy_mat_sums[j] += energy_mat[j][i];
      }
    }
  }
  std::complex<double> biggest_energy =
      (std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
  int biggest_energy_id = 0;
  int iterations = 0;
  double verbose_step = 0.01;
  while (thresh < abs(biggest_energy) && iterations < maxIter) {
    if (double(iterations) / double(maxIter) > verbose_step && verbosity >= 1) {
      std::cout << "\b\b\b\b" << verbose_step * 100 << "%";
      verbose_step += 0.01;
    }

    iterations++;

    // find node with highest energy
    biggest_energy = (0, 0);
    for (int i = 0; i < graph.num_vertices(); i++) {
      if (abs(energy_mat_sums[i]) > abs(biggest_energy)) {
        biggest_energy = energy_mat_sums[i];
        biggest_energy_id = i;
      }
    }
    // move node to better energy position
    double E_dx = biggest_energy.real();
    double E_dy = biggest_energy.imag();
    double E_d2x = 0;
    double E_d2y = 0;
    double E_dxdy = 0;

    for (int i = 0; i < graph.num_vertices(); i++) {
      if (i != biggest_energy_id) {
        std::complex<double> diff =
            _particles[biggest_energy_id].get_pos() - _particles[i].get_pos();
        double denom = 0.001 + pow(std::abs(diff), 3);

        double del_E_d2x = spring_mat[biggest_energy_id][i] *
                           (1 - length_mat[biggest_energy_id][i] *
                                    pow(diff.imag(), 2) / denom);
        double del_E_d2y = spring_mat[biggest_energy_id][i] *
                           (1 - length_mat[biggest_energy_id][i] *
                                    pow(diff.real(), 2) / denom);
        double del_E_dxdy = spring_mat[biggest_energy_id][i] *
                            length_mat[biggest_energy_id][i] * diff.real() *
                            diff.imag() / denom;

        E_d2x += del_E_d2x;
        E_d2y += del_E_d2y;
        E_dxdy += del_E_dxdy;
      }
    }
    double denom = E_d2x * E_d2y - E_dxdy * E_dxdy;
    double dx = (E_dxdy * E_dy - E_d2y * E_dx) / denom;
    double dy = (E_dxdy * E_dx - E_d2x * E_dy) / denom;

    _particles[biggest_energy_id].add_pos(std::complex<double>(dx, dy));

    // update energy matrix
    for (int j = 0; j < graph.num_vertices(); j++) {
      if (biggest_energy_id == j) {
        energy_mat[biggest_energy_id][j] = 0.0;
      } else {
        energy_mat_sums[biggest_energy_id] -= energy_mat[biggest_energy_id][j];
        energy_mat_sums[j] -= energy_mat[j][biggest_energy_id];
        std::complex<double> diff =
            _particles[biggest_energy_id].get_pos() - _particles[j].get_pos();
        double E_dx = spring_mat[biggest_energy_id][j] * (diff.real()) *
                      (1 - length_mat[biggest_energy_id][j] / (std::abs(diff)));
        double E_dy = spring_mat[biggest_energy_id][j] * (diff.imag()) *
                      (1 - length_mat[biggest_energy_id][j] / (std::abs(diff)));
        energy_mat[biggest_energy_id][j] = std::complex<double>(E_dx, E_dy);
        energy_mat[j][biggest_energy_id] = std::complex<double>(-E_dx, -E_dy);
        energy_mat_sums[biggest_energy_id] += energy_mat[biggest_energy_id][j];
        energy_mat_sums[j] += energy_mat[j][biggest_energy_id];
      }
    }
  }
  if (verbosity >= 1) {
    std::cout << "\b\b\b\b100"
              << "%\n";
  }
}

void LayoutGraph::calc_fa2(Graph &graph, int iterations, int _terms,
                           int _thresh, double _edge_force, double _gravity,
                           double _jitter_tol, double _speed,
                           double _speed_efficiency, int verbosity) {

  std::vector<Particle> &_particles = graph.get_particles();

  double speed = _speed;
  double speed_efficiency = _speed_efficiency;
  std::vector<Compl> old_forces;
  if (verbosity >= 2) {
    std::cout << "Calculating fa2 layout for graph with "
              << graph.num_vertices() << " vertices: \n";
  }
  if (_thresh = -1)
    _thresh = std::max(1, int(std::log10(graph.num_vertices())));
  double verbose_step = 0.01;
  std::vector<int> degrees(graph.num_vertices(), 0);
  for (auto e : graph.get_edges()) {
    degrees[e.source] += 1;
    degrees[e.target] += 1;
  }
  Forces _ffm = Forces(_particles);
  for (int i = 0; i < iterations; ++i) {
    if (double(i) / double(iterations) > verbose_step && verbosity >= 1) {
      std::cout << "\b\b\b\b" << verbose_step * 100 << "%";
      verbose_step += 0.01;
    }
    // calc forces
    _ffm.reset(_particles);
    _ffm.calc_ffm_forces(_terms, _thresh);
    _ffm.calc_gravity_forces(_gravity);
    _ffm.calc_edge_forces(_edge_force, graph.get_edges(), degrees);
    std::vector<Compl> forces = _ffm.get_forces();

    if (i == 0)
      old_forces = forces;
    // calc speed
    double total_swing = 0;    // measurement of "erratic" movement
    double total_traction = 0; // measurement of "sensible" movement
    for (int j = 0; j < _particles.size(); j++) {
      if (is_normal(forces[j] - old_forces[j])) {
        total_swing +=
            _particles[j].get_charge() * std::abs(forces[j] - old_forces[j]);
        total_traction += .5 * _particles[j].get_charge() *
                          std::abs(forces[j] + old_forces[j]);
      }
    }
    double estimatedOptimalJitterTolerance = .05 * std::sqrt(_particles.size());
    double minJT = std::sqrt(estimatedOptimalJitterTolerance);
    double maxJT = 10.0;
    double jt =
        _jitter_tol *
        std::max(
            minJT,
            std::min(maxJT, estimatedOptimalJitterTolerance * total_traction /
                                (_particles.size() * _particles.size())));

    double minSpeedEfficiency = 0.05;

    if (total_traction > 0.0 && total_swing / total_traction > 2.0) {
      if (speed_efficiency > minSpeedEfficiency)
        speed_efficiency *= .5;
      jt = std::max(jt, _jitter_tol);
    }
    double target_speed;
    if (total_swing == 0) {
      target_speed = std::numeric_limits<double>::max();
    } else {
      target_speed = jt * speed_efficiency * total_traction / total_swing;
    }
    if (total_swing > jt * total_traction) {
      if (speed_efficiency > minSpeedEfficiency) {
        speed_efficiency *= .7;
      }
    } else if (speed < 1000) {
      speed_efficiency *= 1.3;
    }
    // std::cout<<"speed"<<speed<<total_swing<<"\n";
    speed = speed + std::min(target_speed - speed, 0.5 * speed);

    // apply forces
    for (int j = 0; j < _particles.size(); j++) {
      if (is_normal(speed /
                    (Compl(1, 0) +
                     std::sqrt(_particles[j].get_charge() *
                               std::abs(forces[j] - old_forces[j]) * speed)) *
                    forces[j])) {
        _particles[j].add_pos(
            speed /
            (Compl(1, 0) +
             std::sqrt(_particles[j].get_charge() *
                       std::abs(forces[j] - old_forces[j]) * speed)) *
            forces[j]);
      } else {
        // std::cout<<"\nIS NOT NORMAL "<<speed<<",
        // "<<_particles[j].get_charge()<<", "<<forces[j]<<"\n";
      }
    }
    old_forces = forces;
  }
  if (verbosity >= 1) {
    std::cout << "\b\b\b\b100"
              << "%\n";
  }
}

std::vector<std::vector<Compl>>
LayoutGraph::calc_components_layout(std::vector<std::vector<Compl>> positions,
                                    int iterations) {

  Graph G;
  std::vector<int> component_sizes;

  int id = 0;
  std::vector<double> areas;
  std::vector<Compl> means;
  for (auto l : positions) {
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::min();
    double max_y = std::numeric_limits<double>::min();
    Compl mean = 0;
    for (auto p : l) {
      if (p.real() < min_x) {
        min_x = p.real();
      } else if (p.real() > max_x) {
        max_x = p.real();
      }
      if (p.imag() < min_y) {
        min_y = p.imag();
      } else if (p.imag() > max_y) {
        max_y = p.imag();
      }
      mean += p;
    }
    means.push_back(mean / Compl(l.size(), 0));
    areas.push_back(std::max((max_x - min_x), (max_y - min_y)));
    // or 1 as charge?
    G.add_particle(Particle(std::complex<double>((double)rand() / RAND_MAX,
                                                 (double)rand() / RAND_MAX),
                            1, id),
                   ParticleData());
    id++;
    component_sizes.push_back(l.size());
  }
  double max_area = 0;
  int max_area_id = 0;
  for (int i = 1; i < areas.size(); i++) {
    Edge ee;
    ee.source = max_area_id;
    ee.target = i;
    ee.length = 10 + (areas[max_area_id] + areas[i]) / 1.5;
    G.add_edge(ee);
    if (max_area < areas[i]) {
      max_area_id = i;
      max_area = areas[i];
    }
  }

  calc_kamada_kawai(G, 10, 0.0001, iterations, 2);
  auto _particles = G.get_particles();

  int index = 0;
  std::vector<std::vector<Compl>> result;
  for (int counter = 0; counter < component_sizes.size(); counter++) {
    std::vector<Compl> comp_result;
    for (int sub = 0; sub < component_sizes[counter]; sub++) {
      comp_result.push_back(positions[counter][sub] - means[counter] +
                            _particles[counter].get_pos());
      index++;
    }
    result.push_back(comp_result);
  }
  return result;
}

std::vector<std::vector<Compl>> LayoutGraph::calc_many_components_layout(
    std::vector<std::vector<Compl>> positions, int iterations,
    int comp_per_calc) {
  if (positions.size() < comp_per_calc) {
    return calc_components_layout(positions, iterations);
  } else {
    int counter = 0;
    std::vector<std::vector<std::vector<Compl>>> split_positions;
    std::vector<std::vector<Compl>> split;
    for (auto p : positions) {
      if (p.size() > 500) {
        split_positions.push_back(std::vector<std::vector<Compl>>{p});
        counter--;
      } else {
        split.push_back(p);
      }
      if (counter < comp_per_calc - 1) {
        counter++;
      } else {
        split_positions.push_back(split);
        split.clear();
        counter = 0;
      }
    }
    if (counter != 0) {
      split_positions.push_back(split);
    }
    std::vector<std::vector<Compl>> layouts;
    std::vector<int> sizes;
    for (auto split : split_positions) {
      std::vector<std::vector<Compl>> comp_layout =
          calc_components_layout(split, iterations);
      std::vector<Compl> flat_comp_layout;
      for (auto layout : comp_layout) {
        sizes.push_back(layout.size());
        for (auto position : layout) {
          flat_comp_layout.push_back(position);
        }
      }
      layouts.push_back(flat_comp_layout);
    }
    std::vector<std::vector<Compl>> res =
        calc_many_components_layout(layouts, iterations, comp_per_calc);
    std::vector<Compl> flat_res;
    for (auto rr : res) {
      for (auto r : rr) {
        flat_res.push_back(r);
      }
    }
    std::vector<std::vector<Compl>> result;
    std::vector<Compl> comp_result;
    counter = 0;
    int size_counter = 0;
    for (auto r : flat_res) {
      comp_result.push_back(r);
      if (counter < sizes[size_counter] - 1) {
        counter++;
      } else {
        size_counter += 1;
        counter = 0;
        result.push_back(comp_result);
        comp_result.clear();
      }
    }
    return result;
  }
}
