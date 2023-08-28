// ffm.h
#ifndef CLUSTER_H
#define CLUSTER_H

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

class ClusterNode {
public:
  ClusterNode(int representant, std::vector<int> neighbors,
              std::vector<double> neighbors_distances,
              std::complex<double> repr_pos, double charge) {
    _representant = representant;
    _neighbors = neighbors;
    _neighbors_distances = neighbors_distances;
    _node_ids = {representant};
    _pos = repr_pos;
    _charge = charge;
  };
  std::vector<int> get_neighours() { return _neighbors; }
  std::vector<double> get_neighbors_distances() { return _neighbors_distances; }
  void add_neighbour(int neighbour_id, int distance) {
    if (neighbour_id != _representant) {
      if (std::find(_neighbors.begin(), _neighbors.end(), neighbour_id) ==
          _neighbors.end()) {
        _neighbors.push_back(neighbour_id);
        _neighbors_distances.push_back(distance);
      }
    }
  }
  void remove_neighbour(int id) {
    for (int i = 0; i <= _neighbors.size(); i++) {
      if (_neighbors[i] == id) {
        _neighbors.erase(_neighbors.begin() + i);
        _neighbors_distances.erase(_neighbors_distances.begin() + i);
      }
    }
  }

  int get_representant() { return _representant; }
  std::vector<int> get_node_ids() { return _node_ids; }
  void add_node_id(int node_id) { _node_ids.push_back(node_id); }
  bool is_dead() { return _is_dead; }
  void die() { _is_dead = true; }
  void add_to_cluster(ClusterNode other, int other_distance) {
    for (int n : other.get_node_ids()) {
      _node_ids.push_back(n);
    }
    for (int n = 0; n < other.get_neighours().size(); n++) {
      int id = other.get_neighours()[n];
      if (id != _representant) {
        if (std::find(_neighbors.begin(), _neighbors.end(), id) ==
            _neighbors.end()) {
          _neighbors.push_back(id);
          _neighbors_distances.push_back(other_distance +
                                         other.get_neighbors_distances()[id]);
        }
      }
    }
  }
  std::complex<double> get_pos() { return _pos; }
  double get_charge() { return _charge; }

private:
  bool _is_dead = false;
  int _representant;
  std::vector<int> _node_ids;
  std::vector<int> _neighbors;
  std::vector<double> _neighbors_distances;
  std::complex<double> _pos;
  double _charge;
};

class ClusterPositioner {
public:
  // type 0 leaf cluster, 1 bridge, 2 central
  ClusterPositioner(int type, int id, int neigh1, double dist1) {
    _type = type;
    _id = id;
    _neigh1 = neigh1;
    _distance1 = dist1;
  }
  ClusterPositioner(int type, int id, int neigh1, int neigh2, double dist1,
                    double dist2) {
    _type = type;
    _id = id;
    _neigh1 = neigh1;
    _neigh2 = neigh2;
    _distance1 = dist1;
    _distance2 = dist2;
  }
  int get_id() { return _id; }
  Compl get_position(std::unordered_map<int, Compl> &pos_map,
                     std::vector<int> neighbors = {}) {
    if (_type == 0) {
      // clustered because leaf node
      Compl neighpos = pos_map[_neigh1];
      double angle = (double)rand() / RAND_MAX * 360;
      return Compl(sin(angle), cos(angle)) * Compl(_distance1, 0) + neighpos;

    } else if (_type == 1) {
      // clustered because bridge node
      Compl neigh1pos = pos_map[_neigh1];
      Compl neigh2pos = pos_map[_neigh2];
      return neigh1pos + (neigh2pos - neigh1pos) *
                             Compl(_distance1 / (_distance1 + _distance2), 0);

    } else if (_type == 2) {
      // clustered because next to central node
      Compl neighpos = pos_map[_neigh1];
      Compl center = 0;
      for (auto n : neighbors) {
        center += pos_map[n];
      }
      if (center == Compl(0, 0) || center == pos_map[_id]) {
        Compl neighpos = pos_map[_neigh1];
        double angle = (double)rand() / RAND_MAX * 360;
        return Compl(sin(angle), cos(angle)) * Compl(_distance1, 0) + neighpos;
      } else {
        return center;
      }
    }
    return Compl(0, 0);
  }

private:
  int _id;
  int _type;
  int _neigh1;
  int _neigh2;
  double _distance1;
  double _distance2;
};

class ClusterGraph {

public:
  ClusterGraph(Graph graph) {
    std::vector<std::vector<int>> neighbors(graph.num_vertices());
    for (auto e : graph.get_edges()) {
      neighbors[e.source].push_back(e.target);
      neighbors[e.target].push_back(e.source);
    }
    int i = 0;
    for (auto particle : graph.get_particles()) {
      std::vector<double> distances(neighbors[i].size(), 1.0);
      cluster_nodes.push_back(ClusterNode(particle.get_id(), neighbors[i],
                                          distances, particle.get_pos(),
                                          particle.get_charge()));
      i++;
    }
  };

  void cluster_leaves();
  void cluster_bridges();
  void cluster_central_node(int min_neigh);
  Graph get_graph();
  std::vector<Compl>
  get_unclustered_layout(std::vector<Compl> clustered_layout);

private:
  std::vector<ClusterNode> cluster_nodes;
  std::vector<ClusterPositioner> cluster_positions;
};

#endif
