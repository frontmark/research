
#ifndef QUADTREE_H
#define QUADTREE_H

#include "graph.h"

#include <vector>
#include <tuple>
#include <complex>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <bitset>
#include <iostream>




class Node
{
public:
    Node(int id, int parent, std::vector<Particle> particles, double x0, double y0, double width, double height, int x_coord, int y_coord, int depth): _id(id),_parent(parent), _x0(x0),_y0(y0), _width(width), _height(height), _particles(particles), _is_leaf(true), _incoming_coeff(), _outgoing_coeff(), _x_coord(x_coord), _y_coord(y_coord), _depth(depth) {}

    void split(Node& node1, Node& node2, Node& node3, Node& node4, int id1, int id2, int id3,  int id4);

    double const get_x0()  {return _x0;}
    double const get_y0()  {return _y0;}
    double const get_width()  {return _width;}
    double const get_height()  {return _height;}

    int get_id () {return _id;}

    int const  get_parent(){return _parent;}

    std::vector<int> const get_children() const {
        return {_child_NW, _child_NE, _child_SW, _child_SE};
    }

    const std::vector<int> & get_close_neighbours()const {
        return close_nodes;
    }

    void set_close_neighbours(const std::vector<int>& close){
        close_nodes = close;
    }

    const std::vector<Particle> & get_particles()const  {
        return _particles;
    }

    int const size()const {
        return _particles.size();
    }

    std::complex<double> const get_center()const {
        return std::complex<double>(_x0+_width/2, _y0+_height/2);
    }

    bool const is_leaf()const {
        return _is_leaf;
    }

    void not_leaf_anymore(){
        _is_leaf=false;
    }

    const std::vector< std::complex<double>> & get_outgoing()const {
        return _outgoing_coeff;
    }

    const std::vector< std::complex<double>> &  get_incoming()const {
        return _incoming_coeff;
    }

    
    unsigned int const&  get_x_coord()const {
        return _x_coord;
    }

    unsigned int const& get_y_coord()const {
        return _y_coord;
    }

    int const& get_depth()const {
        return _depth;
    }

    void set_outgoing(std::vector< std::complex<double>> const out){
        _outgoing_coeff = out;
    }

    void set_incoming(std::vector< std::complex<double>> const inc){
        _incoming_coeff = inc;
    }

    void butput (){
        std::cout <<std::endl<<std::endl<<"===============butput============"<< std::endl;
        std::cout << "id "<<_id<< std::endl;
        std::cout << "x0 "<<_x0<< std::endl;
        std::cout << "y0 "<<_y0<< std::endl;
        std::cout << "xcode "<<_x_coord<< std::endl;
        std::cout << "ycode "<<_y_coord<< std::endl;
        std::cout << "width "<<_width<< std::endl;
        std::cout << "height "<<_height<< std::endl;
        std::cout << "depth "<<_depth<< std::endl;
        std::cout << "parent "<<_parent<< std::endl;
        std::cout << "particles "<<_particles.size()<< std::endl;
        for (auto i: _particles){std::cout << i.get_id() << ' ';}
        std::cout << std::endl<< "outgoing_coeff "<< std::endl;
        for (auto i: _outgoing_coeff){std::cout << i << ' ';}
        std::cout << std::endl<< "incoming coeff "<< std::endl;
        for (auto i: _incoming_coeff){std::cout << i << ' ';}
        std::cout << std::endl<< "close neigh "<< std::endl;
        for (auto i: close_nodes){std::cout << i << ' ';}
        std::cout << std::endl<< "isleaf "<<_is_leaf<< std::endl;
        std::cout << "ne "<<_child_NE<< std::endl;
        std::cout << "se "<<_child_SE<< std::endl;
        std::cout << "nw "<<_child_NW<< std::endl;
        std::cout << "sw "<<_child_SW<< std::endl;
    }


private:

    int const _id = -1;

    std::vector<Particle> const _particles;

    
    std::vector<std::complex<double>> _outgoing_coeff;
    std::vector<std::complex<double>> _incoming_coeff;

    unsigned int const _x_coord = 0;
    unsigned int const _y_coord = 0;
    int const _depth = 0;

    std::vector<int> close_nodes = {};
    
    int _child_NW=-1; int _child_NE=-1; int _child_SW=-1; int _child_SE=-1;
    
    int const _parent;
    int _is_leaf;

    double const _x0; double const _y0;
    double const _width; double const _height;
};


class QuadTree
{
public:
    // default constructor 
    QuadTree(){};
    QuadTree(std::vector<Particle> particles)
    {
        double _x0 =    (*std::min_element(particles.begin(),particles.end(), [&](Particle a, Particle b){return a.get_pos().real()<b.get_pos().real();})).get_pos().real();
        double _y0 =    (*std::min_element(particles.begin(),particles.end(), [&](Particle a, Particle b){return a.get_pos().imag()<b.get_pos().imag();})).get_pos().imag();
        double _width = (*std::max_element(particles.begin(),particles.end(), [&](Particle a, Particle b){return a.get_pos().real()<b.get_pos().real();})).get_pos().real()-_x0;
        double _height = (*std::max_element(particles.begin(),particles.end(), [&](Particle a, Particle b){return a.get_pos().imag()<b.get_pos().imag();})).get_pos().imag()-_y0;
        node_map.emplace_back(Node(0, 0, particles, _x0, _y0, _width, _height, 0, 0, 0));
    }

    void split(int const node);
    void calc_close();

    void butput (){for(auto p:node_map){p.butput();}};

    std::vector< Particle> const get_close_particles(int const id);
    std::vector<int> const get_active(int const id);

    Node& get_node(int id){
        return node_map[id];
    }

private:
    std::vector<Node> node_map={};

};

#endif