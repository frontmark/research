#include "cluster.h"
#include "quadtree.h"

#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <bitset>
#include <chrono>
#include <fstream>
#include <string>




void ClusterGraph::cluster_leaves(){
    int id=0; 
    while ( id<=cluster_nodes.size()){

        ClusterNode & cl_node = cluster_nodes[id];
        if (cl_node.is_dead() ==false){
            if (cl_node.get_neighours().size()==1){
                cl_node.die();
                cluster_positions.push_back(ClusterPositioner(0, id, cl_node.get_neighours()[0], cl_node.get_neighbors_distances()[0]));
                ClusterNode neigh = cluster_nodes[cl_node.get_neighours()[0]];
                cl_node.die();
                neigh.add_neighbour(cl_node.get_representant(), cl_node.get_neighbors_distances()[0]);
            }
        }
        id++;
    }
}

void ClusterGraph::cluster_bridges(){
    int id=0; 
    while ( id<=cluster_nodes.size()){
        ClusterNode & cl_node = cluster_nodes[id];
        if (cl_node.is_dead() ==false){
            std::vector<int> neighbors = cl_node.get_neighours();
            if (neighbors.size()==2){
                cl_node.die();
                cluster_positions.push_back(ClusterPositioner(1, id, cl_node.get_neighours()[0], cl_node.get_neighours()[1], cl_node.get_neighbors_distances()[0], cl_node.get_neighbors_distances()[1]));
                ClusterNode & neigh1 = cluster_nodes[cl_node.get_neighours()[0]];
                ClusterNode & neigh2 = cluster_nodes[cl_node.get_neighours()[1]];
                neigh1.add_neighbour(neigh2.get_representant(), cl_node.get_neighbors_distances()[0]+cl_node.get_neighbors_distances()[1]);
                neigh2.add_neighbour(neigh1.get_representant(), cl_node.get_neighbors_distances()[0]+cl_node.get_neighbors_distances()[1]);
                neigh1.remove_neighbour(id);
                neigh2.remove_neighbour(id);
            }
        }
        id++;
    }
}

void ClusterGraph::cluster_central_node(int min_neigh){
    int id=0; 
    while ( id<=cluster_nodes.size()){
        ClusterNode & central_node = cluster_nodes[id];
        if (central_node.is_dead() ==false){
            std::vector<int> neighbors = central_node.get_neighours();
            if (neighbors.size()>=min_neigh){
                int neigh_number = 0;
                for (auto neigh_id: neighbors){    
                    ClusterNode & neigh = cluster_nodes[neigh_id];
                    neigh.die();
                    cluster_positions.push_back(ClusterPositioner(2, neigh_id, id, central_node.get_neighbors_distances()[neigh_number]));
                    central_node.add_to_cluster(neigh, central_node.get_neighbors_distances()[neigh_number]);
                    neigh_number++;
                }
            }
        }
        id++;
    }
}

Graph ClusterGraph::get_graph(){
    Graph g;
    int id=0;
    std::unordered_map<int, int> id_map;
    for (int n = 0; n < cluster_nodes.size(); n++){
        if (cluster_nodes[n].is_dead() == false){
            g.add_particle(Particle(cluster_nodes[n].get_pos(), cluster_nodes[n].get_charge(), id), ParticleData());
            id_map[cluster_nodes[n].get_representant()] = id;
            id++;
        }
    }
    for (auto node: cluster_nodes){
        if (node.is_dead()==false){
            for (int m=0; m<node.get_neighours().size(); m++){
                if (node.get_representant() <node.get_neighours()[m] & cluster_nodes[node.get_neighours()[m]].is_dead()==false){
                    Edge edge;
                    edge.source = id_map[node.get_representant()];
                    edge.target = id_map[node.get_neighours()[m]];
                    edge.length = node.get_neighbors_distances()[m];
                    g.add_edge(edge);
                }
            }
        }
    }
    return g;
}


std::vector<Compl> ClusterGraph::get_unclustered_layout(std::vector<Compl> clustered_layout){
    std::unordered_map<int, Compl> pos_map;

    Compl center = 0;
    for (auto a: clustered_layout){
        center += a;
    }
    center=center/Compl(clustered_layout.size(),0);

    int i = 0;
    for (int n = 0; n < cluster_nodes.size(); n++){
        if (cluster_nodes[n].is_dead() == false){
            pos_map[n]=clustered_layout[i];
            i++;
        }else{
           pos_map[n] = Compl(0,0); 
        }
    }

    for (std::vector<ClusterPositioner>::reverse_iterator i = cluster_positions.rbegin();  i != cluster_positions.rend(); ++i ) { 
            pos_map[i->get_id()] = i->get_position(pos_map, cluster_nodes[i->get_id()].get_neighours()); 

    }

    std::vector<Compl> res;

    for (i=0;i<cluster_nodes.size(); i++){
        res.push_back(pos_map[i]);
    }

    return res;
} 
