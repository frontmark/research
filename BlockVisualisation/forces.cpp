#include "forces.h"
#include "quadtree.h"

#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <bitset>

using Compl = std::complex<double>;



std::uint64_t int_to_key(int a,int b){
    if (a<b){    
        std::uint64_t key = a<<5;
        key |= b;
        return key;
    }
    else{
        std::uint64_t key = b<<5;
        key |= a;
        return key;
    }
}

Compl poly_der_val(const std::vector< Compl>& poly, const Compl& value){
    // Horner to calc derivative of polynomial at value
    int poly_size = poly.size();
    Compl res = poly[poly_size-1]*Compl(poly_size-1,0);
    for (int i = poly_size-2; i>=1;--i){
        res = res * value + poly[i]*Compl(i, 0);
    }
    return res;
}

Compl calc_pair_force(const Particle& p1, const Particle& p2){
    // Calc forces between two particles
    Compl dist = p1.get_pos()-p2.get_pos();
    return std::conj(Compl(p1.get_charge()* p2.get_charge(),0)/dist);
    
}

Forces::Forces(const std::vector< std::pair< Compl, double>>& particles){
    _binom_map={};
    _particles={};
    int i = 0;
    for (auto p: particles){
        _particles.emplace_back(Particle(p.first, p.second, i));
        i=i+1;
    }
    _forces = std::vector<std::complex<double>> (_particles.size(), Compl(0,0));
    _tree = QuadTree(_particles);
}

Forces::Forces(const std::vector<Particle>& particles){
    _binom_map={};
    _particles=particles;
    _forces = std::vector<Compl> (_particles.size(), Compl(0,0));
    _tree = QuadTree(_particles);

}

void Forces::precomputations(){
    
    for (int a=0; a<16; ++a){
        for (int b=0; b<=a; ++b){
            if (a==b){
                _binom_map[int_to_key(a,b)]=Compl(1,0);
            } else if (b==0){
                _binom_map[int_to_key(a,b)]=Compl(1,0);
            } else{
                _binom_map[int_to_key(a,b)]=_binom_map[int_to_key(a-1,b-1)]+_binom_map[int_to_key(a-1,b)];
            }
        }
    }
}


void Forces::split(int const id){
    //split muss vor get_node da die pointer sich aendern koennen in split
    _tree.split(id);
    Node& node = _tree.get_node(id);
    std::vector<int> children  = node.get_children();

    for (int c: children){
        if (_tree.get_node(c).size()>_thresh){
            split(c);
        } else{
            calc_outgoing_coef(c);
        }
    }
    calc_outgoing_coef(id);
}


void Forces::calc_ffm_forces(int terms, int thresh, int const node_id, double strength){
    
    if (node_id==0){
        _terms=terms;
        _thresh=thresh;
        precomputations();
        split(0);
        _tree.calc_close();
    }
    calc_incoming_coef(node_id);
    Node& node = _tree.get_node(node_id);
    //node.butput();
    if (node.is_leaf()){
        
        const std::vector< Compl>  & inc_coef = node.get_incoming();
        std::vector<Particle> const close_particles = _tree.get_close_particles(node_id);
        for (const Particle& p: node.get_particles()){
            //TODO save factor two since we calculate a force between a pair twice...
            Compl shift = p.get_pos() - node.get_center();
            Compl& particle_force = _forces[p.get_id()];
            particle_force = std::conj(poly_der_val(inc_coef, shift));
            //for (auto i: inc_coef){std::cout << i << ' ';}
            for (auto other_p: close_particles){
                particle_force -= calc_pair_force(other_p, p);
            }
            for (auto other_p: node.get_particles()){
                if (p==other_p) continue;
                particle_force -= calc_pair_force(other_p, p);  
            }
            particle_force*=Compl(strength,0);
        }
    } else {
            for (int child_id: node.get_children()){
                calc_ffm_forces(terms, thresh, child_id);
            }
    }
}


void Forces::calc_outgoing_coef(int const id){
    Node& node=_tree.get_node(id);
    Compl center = node.get_center();
    std::vector<Compl> out_coeff(_terms+1, Compl(0,0));
    if (node.is_leaf()){
        std::vector<Particle> particles = node.get_particles();
        Compl charge_sum = Compl(0,0);
        for(auto p: particles){
            charge_sum += p.get_charge();
        }
        out_coeff[0]=charge_sum;
        for (int k=1; k<_terms+1; ++k){
            Compl particle_am = Compl(0,0);
            for(auto p: particles){
                particle_am += -p.get_charge()*pow(p.get_pos()-center,k)/Compl(k,0);
            }
            out_coeff[k]=particle_am;
        }
    } else {
        for (int child_id: node.get_children()){
            Node& child = _tree.get_node(child_id);
            Compl shift = child.get_center()-center;
            std::vector<Compl> const& child_coeff = child.get_outgoing();
            out_coeff[0] += child_coeff[0];
            for (int l=1; l<_terms+1; ++l){
                Compl res = 0;
                for (int k=1; k<l; ++k){
                     res += child_coeff[k]*pow(shift, l-k)*_binom_map[int_to_key(l-1, k-1)]-child_coeff[0]*pow(shift,l)/Compl(l,0);
                }
                out_coeff[l] += res;
            }
        }
    }
    _tree.get_node(id).set_outgoing(out_coeff);
}


void Forces::calc_incoming_coef(int const id){
    
    Node& node=_tree.get_node(id);
    const Compl& center = node.get_center();
    std::vector<Compl> in_coeff(_terms+1, Compl(0,0));
    Compl& in_coeff_at_0 = in_coeff[0];
    
    std::vector<Compl> temp_k_shift_facs(_terms+1, Compl(0,0));
    for (int act_id: _tree.get_active(id)){

        Node& act_node = _tree.get_node(act_id);
        Compl shift = act_node.get_center()-center;
        const std::vector< Compl> & act_coeff = act_node.get_outgoing();

        for (int k=1; k<_terms+1; ++k){ 
            temp_k_shift_facs[k] = act_coeff[k]/pow(shift,k)*(k & 1==1? Compl(-1, 0): Compl(1, 0));
        }

        for (int l=1; l<_terms+1; ++l){
            Compl powlshift = pow(shift,l);
            in_coeff_at_0 += act_coeff[l]/powlshift*(l & 1==1? Compl(-1, 0): Compl(1, 0));
            Compl& in_coeff_at_l = in_coeff[l];
            for (int k=1; k<_terms+1; ++k){
                in_coeff_at_l += temp_k_shift_facs[k]/powlshift*_binom_map[int_to_key(l+k-1, k-1)];
            }
            in_coeff_at_l -= act_coeff[0]/(powlshift*Compl(l,0));
        }
        in_coeff_at_0 += act_coeff[0]*std::log(-shift);
    }
    int parent_id = node.get_parent();
    if (parent_id!=node.get_id()){
        Node& parent =_tree.get_node(parent_id);
        Compl shift = parent.get_center()-center;
        const std::vector< Compl>& parent_coeff = parent.get_incoming();
        for (int l=0; l<_terms+1; ++l){
            Compl res = 0;
            for (int k=1; k<_terms+1; ++k){
                res += parent_coeff[k]*_binom_map[int_to_key(k,l)]*pow(shift, k-l)*(k-l & 1==1? Compl(-1, 0): Compl(1, 0));
            }
            in_coeff[l] += res;
        }
    }
    node.set_incoming(in_coeff);
}
