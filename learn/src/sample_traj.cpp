#include <iostream>
#include <bmla>
#include <fstream>
#include <vector>
#include <map>

#include <sstream>

#include "dbm_ic.hpp"

using namespace dblz;
using namespace std;

// ***************
// MARK: - Read
// ***************

std::vector<double> read_ixn_traj(string fname) {
    std::vector<double> ret;
    
    std::ifstream f;
    f.open(fname);
    if (!f.is_open()) {
        std::cerr << ">>> read_ic <<< could not find file: " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string timepoint="";
    std::string val="";
    std::string line;
    std::istringstream iss;
    if (f.is_open()) { // make sure we found it
        while (getline(f,line)) {
            if (line == "") { continue; };
            iss = std::istringstream(line);
            iss >> timepoint;
            iss >> val;
            ret.push_back(atof(val.c_str()));
            timepoint=""; val="";
        };
    };
    f.close();
    
    return ret;
};


int main() {

    // ***************
    // MARK: - DBM
    // ***************
    
    auto dbm = DBM();

    // ***************
    // MARK: - Setup optimizer
    // ***************
    
    // No markov chains
    int no_markov_chains = 100;
    
    std::string dir = "../data/sample_traj/";
    
    dbm.latt->set_no_markov_chains(MCType::AWAKE, no_markov_chains);
    dbm.latt->set_no_markov_chains(MCType::ASLEEP, no_markov_chains);

    // ***************
    // MARK: - Read
    // ***************
    
    int timepoint_read = 9800;
    auto traj_hX = read_ixn_traj("../data/learn_params/ixn_params/hX_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_hY = read_ixn_traj("../data/learn_params/ixn_params/hY_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wXX1 = read_ixn_traj("../data/learn_params/ixn_params/wXX1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wYY1 = read_ixn_traj("../data/learn_params/ixn_params/wYY1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bX1 = read_ixn_traj("../data/learn_params/ixn_params/bX1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bY1 = read_ixn_traj("../data/learn_params/ixn_params/bY1_"+pad_str(timepoint_read, 5)+".txt");
    
    // ***************
    // MARK: - Go through all times
    // ***************
    
    int layer = 0;
    int no_cd_steps = 10;
    MCType type = MCType::ASLEEP;
    int no_timesteps=500;
    
    for (auto timepoint=0; timepoint<no_timesteps; timepoint++) {
        
        std::cout << timepoint << " / " << no_timesteps << std::endl;
        
        // Set vals
        dbm.bias.at(0).at("hX")->set_val(traj_hX.at(timepoint));
        dbm.bias.at(0).at("hY")->set_val(traj_hY.at(timepoint));
        dbm.w_ixns.at(0).at("wXX1")->set_val(traj_wXX1.at(timepoint));
        dbm.w_ixns.at(0).at("wYY1")->set_val(traj_wYY1.at(timepoint));
        dbm.bias.at(1).at("bX1")->set_val(traj_bX1.at(timepoint));
        dbm.bias.at(1).at("bY1")->set_val(traj_bY1.at(timepoint));

        // Random initial
        layer = 0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->set_random_all_units_in_layer(MCType::ASLEEP, i_chain, layer, true);
        };
        
        // Up
        layer = 1;
        dbm.latt->activate_layer_calculate_from_below(type, layer);
        dbm.latt->activate_layer_convert_to_probs(type, layer, true);
        dbm.latt->activate_layer_committ(type, layer);
        
        // Sample
        for (auto i_step=0; i_step<no_cd_steps-1; i_step++) {
            
            // Down
            layer = 0;
            dbm.latt->activate_layer_calculate_from_above(type, layer);
            dbm.latt->activate_layer_convert_to_probs(type, layer, true);
            dbm.latt->activate_layer_committ(type, layer);
            
            // Up
            layer = 1;
            dbm.latt->activate_layer_calculate_from_below(type, layer);
            dbm.latt->activate_layer_convert_to_probs(type, layer, true);
            dbm.latt->activate_layer_committ(type, layer);
        };
        
        // Down
        layer = 0;
        dbm.latt->activate_layer_calculate_from_above(type, layer);
        dbm.latt->activate_layer_convert_to_probs(type, layer, true);
        dbm.latt->activate_layer_committ(type, layer);
        // Write
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(type, i_chain, layer, dir+pad_str(timepoint, 4)+"/layer_0_"+pad_str(i_chain, 4)+".txt", true);
        }
        
        // Up with prob
        layer = 1;
        dbm.latt->activate_layer_calculate_from_below(type, layer);
        dbm.latt->activate_layer_convert_to_probs(type, layer, false);
        dbm.latt->activate_layer_committ(type, layer);
        // Write
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(type, i_chain, layer, dir+pad_str(timepoint, 4)+"/layer_1_"+pad_str(i_chain, 4)+".txt", false);
        };

    };
    
	return 0;
};

