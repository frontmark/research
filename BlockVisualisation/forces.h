// ffm.h
#ifndef FFM_H
#define FFM_H

#include "quadtree.h"

#include <vector>
#include <tuple>
#include <complex>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <bitset>


using Compl = std::complex<double>;

struct hash_pair {
    size_t operator()(std::pair<Compl, double> a)
    {
        return std::hash<double>{}(a.first.imag()) ^ std::hash<double>{}(a.first.real()) ^ std::hash<double>{}(a.second);          
    }
};

class Forces
{

public:
    Forces(){};
    Forces(const std::vector<std::pair<Compl, double>>& particles);
    Forces(const std::vector<Particle>& particles);
    
    void calc_ffm_forces(int terms=3, int thresh=10, int const node_id=0, double strength=.5);

    void calc_gravity_forces(double gravity){
        // centers around 0
        for (Particle p: _particles){
            double dist = std::abs(p.get_pos());
            if(dist>0.0001){
                _forces[p.get_id()]-=gravity*p.get_pos()*p.get_charge()/std::sqrt(dist);
            }
        }
    };

    void calc_gravity_forces2(double gravity){
        // centers around 0
        for (Particle p: _particles){
            double dist = std::abs(p.get_pos());
            if(dist>0.0001){
                _forces[p.get_id()]-=gravity*p.get_pos()*p.get_charge();
            }
        }
    };

    void calc_edge_forces(double edge_factor, const std::vector<Edge>& edges, std::vector<int> degree){
        for (auto edge: edges){
            double distance = std::max(std::abs(_particles[edge.source].get_pos()-_particles[edge.target].get_pos()), 0.0001);
            double opt_edge_length = 1;
            Compl force = double(degree[edge.source]+1)*double(degree[edge.target]+1)*Compl(edge_factor,0)*(_particles[edge.source].get_pos()-_particles[edge.target].get_pos())*(1.0-distance)/distance;
            _forces[edge.source] += force;
            _forces[edge.target] -= force;
        }
    };

    void calc_edge_forces_2(double edge_factor, const std::vector<Edge>& edges){
        for (auto edge: edges){
            double distance = std::max(std::abs(_particles[edge.source].get_pos()-_particles[edge.target].get_pos()), 0.0001);
            double opt_edge_length = 1;
            Compl force = Compl(edge_factor,0)*(_particles[edge.source].get_pos()-_particles[edge.target].get_pos())*(1.0-distance)/distance;
            _forces[edge.source] += force;
            _forces[edge.target] -= force;
        }
    };

    void update(std::vector<Particle>& new_parts){

        _particles = new_parts;
        _forces = std::vector<Compl> (_particles.size(), Compl(0,0));
        _tree = QuadTree(_particles);
    }

    std::vector<Compl> get_forces(){return _forces;};

private:
    void precomputations();
    void split(int const id);
    void calc_outgoing();
    void calc_incoming();
    
    void calc_outgoing_coef(int const id);
    void calc_incoming_coef(int const id);

    int _terms = 0;
    int _thresh = 0;
    std::vector<Compl> _forces;
    std::vector<Particle> _particles;
    QuadTree _tree;

    std::unordered_map<std::bitset<8>, Compl> _binom_map;
};



#endif