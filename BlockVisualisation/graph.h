
#ifndef GRAPH_H
#define GRAPH_H


#include <vector>
#include <tuple>
#include <complex>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <bitset>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>


struct ParticleData {
    bool type; //true for address, false for transactions
    int timestamp = 0;
    std::string label;
    std::string address = "NONE";
    std::string txid = "NONE";
};

struct Edge {
    int source;
    int target;
    std::string type;
    double value;
    double length = 1.0;
};


class Particle
{
public:
    Particle(std::complex<double> pos, double charge, int id): _pos(pos), _charge(charge), _id(id){}

    double const get_charge() const {return _charge;}
    const std::complex<double>& get_pos() const {return _pos;}
    int const get_id () const {return _id;}
    void set_id (int id) {_id=id;}
    bool operator==(const Particle& rhs)const  {return _id == rhs.get_id(); }

    void add_pos (std::complex< double> _add) {_pos += _add;}
    void mul_pos (std::complex< double> _add) {_pos *= _add;}

private:
    std::complex< double> _pos;
    double _charge;
    int _id;
};

class Graph
{
public:

    Graph(){};

    Graph(std::vector<Particle> particles, std::vector<ParticleData> particle_data, std::vector<Edge> edges){
        _particles=particles;
        _particle_data = particle_data;
        _edges=edges;
    }

    Graph(std::string file_type, std::string file_path){
        std::cout<<"Reading file "<<file_path<<".\n";
        if (file_type=="standard"){
            // reads file. first line is int with number of nodes and then each line has two ints the endpoints of an edge
            std::ifstream infile(file_path);
            std::string str; 
            std::getline(infile, str);
            int id_number = std::stoi(str);
            for (int i = 0; i<id_number;i++){
                _particles.emplace_back(Particle(std::complex<double>((double)rand() / RAND_MAX,(double)rand() / RAND_MAX), 1, i));
                _particle_data.push_back(ParticleData());
            }
            while (std::getline(infile, str)){
                int slash_index = str.find(" ");
                Edge edge;
                edge.source = std::stoi(str.substr(0, slash_index));
                edge.target = std::stoi(str.substr(slash_index + 1, str.length()));
                _edges.push_back(edge);
            }
        } else if (file_type == "gml"){
             // reads gml file into grpah format
            std::ifstream infile(file_path);
            std::string str; 
            std::getline(infile, str);
            int reading_node_edge = 0;
            ParticleData particle_data;
            Edge edge;
            int id_counter = 0;
            std::map<int, int> id_map;
            
            while (std::getline(infile, str)){
                if (str == "  node ["){
                    reading_node_edge=1;
                    id_counter += 1;
                } else if (str == "  edge ["){
                    reading_node_edge=2;
                }
                if (reading_node_edge==1){
                    if (str.substr(0, 9)=="    label"){
                        particle_data.label = str.substr(10, str.size());
                    } else if (str.substr(0, 7)=="    id "){
                        id_map[std::stoi(str.substr(7, str.size()))] = _particles.size();
                    }else if (str.substr(0, 9)=="    txid "){
                        particle_data.txid = str.substr(9, str.size());
                        particle_data.type = false;
                    } else if (str.substr(0, 12)=="    address "){
                        particle_data.address = str.substr(12, str.size());
                        particle_data.type = true;
                    } else if (str.substr(0, 14)=="    timestamp "){
                        particle_data.timestamp = std::stoi(str.substr(14, str.size()));
                    } else if (str == "  ]"){
                        reading_node_edge=0;
                        _particle_data.push_back(particle_data);
                        _particles.emplace_back(Particle(std::complex<double>((double)rand() / RAND_MAX,(double)rand() / RAND_MAX), 1, _particles.size()));
                        particle_data = ParticleData();
                    }
                } else if (reading_node_edge==2){
                    
                    if (str.substr(0, 11)=="    source "){
                        edge.source = std::stoi(str.substr(10, str.size()));
                    } else if (str.substr(0, 11)=="    target "){
                        edge.target = std::stoi(str.substr(10, str.size()));
                    } else if (str.substr(0, 9)=="    type "){
                        edge.type = str.substr(9, str.size());
                    } else if (str.substr(0, 10)=="    value "){
                        edge.value = std::stod(str.substr(10, str.size()));
                    } else if (str == "  ]"){
                        try{
                            reading_node_edge=0;
                            _edges.push_back(edge);
                            edge = Edge();
                        }catch (const std::out_of_range& ex) {}
                    }
                }
            }
        } else {
            std::cout << " Graph init recieved bad file format, valid formats are: standard, gml";
        }
    }

    int num_edges(){
        return _edges.size();
    }

    int num_vertices(){
        return _particles.size();
    }

    void add_particle(Particle p, ParticleData pd){
        _particles.push_back(p);
        _particle_data.push_back(pd);
    }

    void add_edge(Edge e){
        _edges.push_back(e);
    }

    std::vector<std::vector<double>> floyd_warschall(){

        std::vector<std::vector<double>> distances(_particles.size(), std::vector<double>(_particles.size(), std::numeric_limits<double>::max()));

        // init distances matrix with zeros for nodes to itself and edge length from nodes to neighbours
        for (int id = 0; id < _particles.size(); id++) {
            distances[id][id] = 0;
        }
        for (auto e: _edges){
            distances[e.source][e.target] = e.length;
            distances[e.target][e.source] = e.length;
        }

        // floyd warshall
        for (int k = 0; k < _particles.size(); k++) {
            for (int i = 0; i < _particles.size(); i++) {
                for (int j = 0; j < _particles.size(); j++) {
                    if (distances[i][k] > 0 && distances[k][j] > std::numeric_limits<double>::max() - distances[i][k]) {
                        distances[i][j] = std::min(distances[i][j], std::numeric_limits<double>::max());
                    } else{
                        distances[i][j] = std::min(distances[i][j], distances[i][k] + distances[k][j]);
                    }
                }
            }
        }

        return distances;

    }

    std::tuple<std::vector<Graph>, std::unordered_map<int, int>, std::vector<int>>  get_components(){

        std::vector<std::vector<int>> neighbours(_particles.size());

        std::vector<int> parent_vec;
        for (int i=0;i<_particles.size(); i++){
            parent_vec.push_back(-1);
        }

        for(auto e: _edges){
            neighbours[e.source].push_back(e.target);
            neighbours[e.target].push_back(e.source);
        }

        class get
        {
        public:
            static void get_neighbours(int node, int val, std::vector<std::vector<int>>& neighbours, std::vector<int>& parent_vec) { 
                if (parent_vec[node]== -1){
                    parent_vec[node] = val;
                    for (auto n: neighbours[node]){
                        get_neighbours(n, val, neighbours, parent_vec);
                    }
                }
            }
        };

        for (int i=0;i<_particles.size(); i++){
            get::get_neighbours(i, i, neighbours, parent_vec);
        }

        int num_components = std::set<double>( parent_vec.begin(), parent_vec.end() ).size();
        //for (auto i: parent_vec) std::cout << i << ' ';
        std::unordered_map<int, Graph> _components;
        std::unordered_map<int, int> _new_id;
        std::unordered_map<int, int> graph_number;
        int counter=0;
        for (int i: std::set<double>( parent_vec.begin(), parent_vec.end() )){
            _components[i] = Graph();
            graph_number[i]=counter;
            counter++;
        }
        for (int i=0;i<parent_vec.size();i++){
            int new_id = _components[parent_vec[i]].num_vertices();
            _components[parent_vec[i]].add_particle(Particle(_particles[i].get_pos(), _particles[i].get_charge(),new_id), _particle_data[i]);
            _new_id[i] = new_id;
        }
        for (auto e: _edges){
            int component = parent_vec[e.source];
            e.source = _new_id[e.source];
            e.target = _new_id[e.target];
            _components[component].add_edge(e);
        }

        auto value_selector = [](auto pair){return pair.second;};
        //std::vector<Graph> components(_components.size());
        //std::transform(_components.begin(), _components.end(), components.begin(), value_selector);
        std::vector<Graph> components;
        std::set<int> unique_parent_vec;
        for (auto p: parent_vec){
            if(unique_parent_vec.find(p) == unique_parent_vec.end()){
                components.push_back(_components[p]);
                unique_parent_vec.insert(p);
            }
        } 

        std::vector<int> component_nr;
        for (int i=0;i<parent_vec.size();i++){
            component_nr.push_back(graph_number[parent_vec[i]]);
        }

        return std::tuple<std::vector<Graph>, std::unordered_map<int, int>, std::vector<int>> (components, _new_id, component_nr); 
    }

    std::set<int> get_timestampset(){
        std::set<int> timestamps;
        for (auto p:_particle_data){
            if (p.timestamp != 0){
                timestamps.insert(p.timestamp);
            }
        }
        return timestamps;
    }

    Graph timestamp_subgraph(int timestamp){
        
        std::vector<bool> part_mask(_particles.size(), false); 
        std::vector<bool> edge_mask(_edges.size(),false); 
        for (int i=0;i<_edges.size(); i++){
            if (_particle_data[_edges[i].source].timestamp==timestamp || 
                    _particle_data[_edges[i].target].timestamp==timestamp){
                edge_mask[i]=true;
                part_mask[_edges[i].source] = true;
                part_mask[_edges[i].target] = true;
            }
        }
        Graph subgraph = Graph();
        std::unordered_map<int, int> _new_id;
        for (int i=0;i<_particles.size();i++){
            if (part_mask[i]){
                int new_id = subgraph.num_vertices();
                subgraph.add_particle(Particle(_particles[i].get_pos(), _particles[i].get_charge(),new_id), _particle_data[i]);
                _new_id[i] = new_id;
            }
        }
        for (int i=0;i<_edges.size();i++){
            if(edge_mask[i]){
                _edges[i].source = _new_id[_edges[i].source];
                _edges[i].target = _new_id[_edges[i].target];
                subgraph.add_edge(_edges[i]);
            }
        }
        return subgraph;
    }


    std::vector<Particle>& get_particles(){ return _particles;}
    std::vector<ParticleData>& get_particle_data(){ return _particle_data;}
    std::vector<Edge>& get_edges(){ return _edges;}

private:
    std::vector<Particle> _particles = {};
    std::vector<ParticleData> _particle_data = {};
    std::vector<Edge> _edges;
};


#endif