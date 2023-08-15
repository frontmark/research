#include "block_visual.h"
#include "layout.h"
#include "cluster.h"

#include <map>



BlockVisual::BlockVisual(std::string file_settings = ""){
    _settings = Settings();
    if (file_settings != ""){
        std::ifstream setfile(file_settings.c_str());
        std::string str;
        std::getline(setfile, str);std::getline(setfile, str);std::getline(setfile, str);std::getline(setfile, str);
        _file_path_in = str;
        std::getline(setfile, str);std::getline(setfile, str);
        _file_path_out = str;
        std::getline(setfile, str);std::getline(setfile, str);std::getline(setfile, str);
        _settings.componentwise = std::stoi(str);
        std::getline(setfile, str); std::getline(setfile, str);
        _settings.iterations = std::stoi(str);
        std::getline(setfile, str); std::getline(setfile, str);
        _settings.edge_force = std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.gravity = std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.jitter_tol = std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.speed = std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.speed_efficiency = std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.layout_verbosity = std::stoi(str);
        std::getline(setfile, str);std::getline(setfile, str);std::getline(setfile, str);
        _settings.terms = std::stoi(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.thresh = std::stoi(str);
        std::getline(setfile, str);std::getline(setfile, str);std::getline(setfile, str);
        _settings.kk_edge_strength =  std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.kk_thresh = std::stod(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.kk_iterations = std::stoi(str);
        std::getline(setfile, str);std::getline(setfile, str);
        _settings.kk_verbosity = std::stoi(str);
        std::getline(setfile, str);std::getline(setfile, str);std::getline(setfile, str);
        _settings.max_kamada_kawai_nodes = std::stoi(str);
    }
}


void write_particles(std::vector<Particle> _particles, int component, std::string filename){

    if (component == -1){
        std::ofstream myfile;
        myfile.open((filename+".txt").c_str(), std::ios::out|std::ios::trunc);
        for (Particle p: _particles){myfile << p.get_pos().real() << "|" << p.get_pos().imag() << std::endl;}
        myfile.close();
    } else {
        std::ofstream myfile;
        myfile.open((filename+".txt").c_str(), std::ios::out|std::ios::app);
        int p_index=0;
        for (Particle p: _particles){myfile <<component << "|" << p.get_pos().real() << "|" << p.get_pos().imag() << std::endl;p_index++;}
        myfile.close();
    }

};

void write_positions(std::vector<Compl> _positions, std::string filename){


    std::ofstream myfile;
    myfile.open((filename+".txt").c_str());
    for (Compl p: _positions){myfile<< p.real() << "|" << p.imag() << std::endl;}
    myfile.close();

};


void BlockVisual::calc_visual(){
    std::string type = "standard";
    if (_file_path_in.find("gml") != std::string::npos) {
        type="gml";
    }
    Graph big_graph = Graph(type, _file_path_in);


    LayoutGraph l;

    std::tuple<std::vector<Graph>, std::unordered_map<int, int>, std::vector<int>> components = big_graph.get_components();
    std::vector<Graph> graphs = std::get<0>(components);
    std::unordered_map<int, int> new_id = std::get<1>(components);
    std::vector<int> component_vec = std::get<2>(components);

    std::vector<std::vector<Compl>> component_positions;

    for (auto graph: graphs){

        ClusterGraph CL = ClusterGraph(graph);
        Graph g = CL.get_graph();
        if (g.num_vertices()<_settings.max_kamada_kawai_nodes){
            //l.calc_kamada_kawai(g, 10000, 0.0001, 200, 2);
            l.calc_kamada_kawai(g, _settings.kk_edge_strength, _settings.kk_thresh, _settings.kk_iterations,_settings.kk_verbosity);
            std::vector<Compl> clustered_layout;
            for (auto p: g.get_particles()){
                clustered_layout.push_back(p.get_pos());
            }
            component_positions.push_back(clustered_layout);
        }else{
            l.calc_fa2(g, _settings.iterations, _settings.terms, _settings.thresh,_settings.edge_force, _settings.gravity, _settings.jitter_tol, _settings.speed, _settings.speed_efficiency, _settings.layout_verbosity);

            //CL.cluster_leaves();
            //CL.cluster_bridges();
            //CL.cluster_central_node(4);

            std::vector<Compl> clustered_layout;
            double mean_edge = 0;
            for (auto e: g.get_edges()){
                mean_edge+=std::abs(g.get_particles()[e.source].get_pos()-g.get_particles()[e.target].get_pos());
            }
            mean_edge /= g.num_edges();
            for (auto p: g.get_particles()){
                clustered_layout.push_back(p.get_pos()/mean_edge);
            }
            component_positions.push_back(clustered_layout);
            //std::vector<Compl> decluster_positions = CL.get_unclustered_layout(clustered_layout);
        }



    }

    component_positions = l.calc_many_components_layout(component_positions, component_positions.size()*15, 250);
    std::vector<Compl> positions;
    for (int id=0;id<big_graph.num_vertices();id++){
        positions.push_back(component_positions[component_vec[id]][new_id[id]]);
    }
    write_positions(positions, _file_path_out);
}


void BlockVisual::fa2_visual(){
        std::string type = "standard";
        if (_file_path_in.find("gml") != std::string::npos) {
            type="gml";
        }
        Graph big_graph = Graph(type, _file_path_in);


        LayoutGraph l;
        l.calc_fa2(big_graph, _settings.iterations, _settings.terms, _settings.thresh,_settings.edge_force, _settings.gravity, _settings.jitter_tol, _settings.speed, _settings.speed_efficiency, _settings.layout_verbosity);
        std::vector<Compl> positions;
        for (auto p: big_graph.get_particles()){
            positions.push_back(p.get_pos());
        }
        write_positions(positions, _file_path_out);
    };
