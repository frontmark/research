// A2DD.cpp

#include "graph.h"
#include "quadtree.h"

#include <vector>
#include <tuple>
#include <complex>
#include <unordered_map>
#include <algorithm>
#include <iostream>

using Compl = std::complex<double>;

int bit_ceil2(unsigned int a){
    int c = 0;
    while (a!= 0){
        a = a>>1;
        c+=1;
    }
    return c;
};


// Represents a node in the quadtree.
void Node::split(Node& node1, Node& node2, Node& node3, Node& node4, int id1, int id2, int id3,  int id4)
{
    _child_NW = id1; _child_NE = id2; _child_SW = id3; _child_SE = id4;
}



std::vector< Particle> const QuadTree::get_close_particles(int const id){
    Node& node = get_node(id);
    std::vector<Particle> vec;
    for (auto n_id: node.get_close_neighbours()){
        //if (n_id == -1) continue;
        Node& cl_node = get_node(n_id);
        for (auto particle: get_node(n_id).get_particles()){
            vec.push_back(particle);
        }
    }
    return vec;
}

std::vector<int> const QuadTree::get_active(int const id){
    //top level has no active neighbours...
    if (id<4){return {};}

    std::vector<int> vec;
    Node& node  = get_node(id);
    auto  my_close = node.get_close_neighbours();
    auto parent_close = get_node(node.get_parent()).get_close_neighbours();
    for (int p_close: parent_close){
        Node& close_Node = node_map[p_close];
        if(close_Node.is_leaf()){
            if(std::find(my_close.begin(), my_close.end(), p_close) == my_close.end()) {
                vec.push_back(p_close);
            }
        } else{
            for (int c: close_Node.get_children()){
                if(std::find(my_close.begin(), my_close.end(), c) == my_close.end()) {
                    vec.push_back(c);
                }
            } 
        }
    }
    return vec;
}



void QuadTree::split(int id)
{
    
    Node& node = get_node(id);

    double const w2 = node.get_width()/2;
    double const h2 = node.get_height()/2;
    double const x0 = node.get_x0();
    double const y0 = node.get_y0();

    int x_coord = node.get_x_coord()<<1;
    int y_coord = node.get_y_coord()<<1;
    int const& depth = node.get_depth()+1;
    
    std::vector<Particle> particles1;
    std::vector<Particle> particles2;
    std::vector<Particle> particles3;
    std::vector<Particle> particles4;
    
    for (auto p: node.get_particles()){
        if (p.get_pos().real()< x0+w2){
            if (p.get_pos().imag()>= y0+h2){
                particles3.push_back(p);
                continue;
            } else{
                particles1.push_back(p);
                continue;
            }
        }
        else {
            if (p.get_pos().imag()>= y0+h2){
                particles4.push_back(p);
                continue;
            } else{
                particles2.push_back(p);
                continue;
            }
        }
    }

    int id1 = node_map.size();
    int id2 = id1+1;
    int id3 = id2+1;
    int id4 = id3+1;

    int y_coord1 = y_coord|1;
    int x_coord1 = x_coord|1;

    node_map.emplace_back(Node (id1, id, particles1, x0, y0, w2, h2, x_coord, y_coord, depth));
    node_map.emplace_back(Node (id2, id, particles2, x0+w2, y0, w2, h2, x_coord1, y_coord, depth));
    node_map.emplace_back(Node (id3, id, particles3, x0, y0+h2, w2, h2, x_coord, y_coord1, depth));
    node_map.emplace_back(Node (id4, id, particles4, x0+w2, y0+h2, w2, h2, y_coord1, y_coord1, depth));

    Node& node1 = node_map[id1];
    Node& node2 = node_map[id2];
    Node& node3 = node_map[id3];
    Node& node4 = node_map[id4];
    Node& nodeP = get_node(id);

    nodeP.split(node1, node2, node3, node4, id1, id2, id3, id4);
    nodeP.not_leaf_anymore();
    
}


//TODO probably COMPLETTTTEEEE TRASH :/
void QuadTree::calc_close()
{
    int size = node_map.size();
    for (int i = 0; i < size; i++){
        Node& p = node_map[i];
        int const& depth = p.get_depth();
        unsigned int const& x_coord = p.get_x_coord();
        unsigned int const& y_coord = p.get_y_coord();
        std::vector<std::pair<int, int>> const shifts = {{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}};
        std::vector<int> close = {};
        for (std::pair<int, int> shift: shifts){
            if (((shift.first == -1 & x_coord>=1) || (shift.first == 1 & bit_ceil2(x_coord+1)<=depth) || shift.first == 0)&
                ((shift.second == -1 & y_coord>=1) || (shift.second == 1 & bit_ceil2(y_coord+1)<=depth) || shift.second == 0)){
                unsigned int close_x_coord = x_coord+shift.first;
                unsigned int close_y_coord = y_coord+shift.second;
                int counter = 1;
                int temp = 0;
                while (node_map[temp].is_leaf() == false & depth-counter>=0){
                    if (0 == ( (close_x_coord >> depth-counter) & 1)){
                        if (0 == ( (close_y_coord >> depth-counter) & 1)){
                            temp = node_map[temp].get_children()[0];
                        } else{
                            temp = node_map[temp].get_children()[2];
                        }
                    } else {
                        if (0 == ( (close_y_coord >> depth-counter) & 1)){
                            temp = node_map[temp].get_children()[1];
                        } else{
                            temp = node_map[temp].get_children()[3];
                        }
                    }
                    counter +=1;
                }
                if(std::find(close.begin(), close.end(), temp) == close.end()) {
                    close.push_back(temp);
                } 
            }
        }
        node_map[i].set_close_neighbours(close);
    }
}