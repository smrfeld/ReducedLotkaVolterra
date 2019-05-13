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

// ***************
// MARK: - Write moments
// ***************

void write_moments(string fname, double moment_hX, double moment_hY, double moment_wXX1, double moment_wYY1, double moment_bX1, double moment_bY1) {
    
    std::ofstream f;
    f.open (fname);
    if (!f.is_open()) { // make sure we found it
        std::cerr << ">>> write_moments <<< Error: could not open file: " << fname << " for writing" << std::endl;
        exit(EXIT_FAILURE);
    };
    
    f << "hX " << moment_hX << "\n";
    f << "hY " << moment_hY << "\n";
    f << "wXX1 " << moment_wXX1 << "\n";
    f << "wYY1 " << moment_wYY1 << "\n";
    f << "bX1 " << moment_bX1 << "\n";
    f << "bY1 " << moment_bY1 << "\n";
    
    f.close();
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
    
    std::string dir_read = "learn_params_rbm_centered";
    
    int timepoint_read = 9800;
    auto traj_hX = read_ixn_traj("../data/"+dir_read+"/ixn_params/hX_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_hY = read_ixn_traj("../data/"+dir_read+"/ixn_params/hY_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wXX1 = read_ixn_traj("../data/"+dir_read+"/ixn_params/wXX1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wYY1 = read_ixn_traj("../data/"+dir_read+"/ixn_params/wYY1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bX1 = read_ixn_traj("../data/"+dir_read+"/ixn_params/bX1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bY1 = read_ixn_traj("../data/"+dir_read+"/ixn_params/bY1_"+pad_str(timepoint_read, 5)+".txt");
    
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
            dbm.latt->write_layer_to_file(type, i_chain, layer, dir+pad_str(timepoint, 4)+"/lattices/layer_0_"+pad_str(i_chain, 4)+".txt", true);
        }
        
        // Up with prob
        layer = 1;
        dbm.latt->activate_layer_calculate_from_below(type, layer);
        dbm.latt->activate_layer_convert_to_probs(type, layer, false);
        dbm.latt->activate_layer_committ(type, layer);
        // Write
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(type, i_chain, layer, dir+pad_str(timepoint, 4)+"/lattices/layer_1_"+pad_str(i_chain, 4)+".txt", false);
        };

        // ***************
        // MARK: - Reap moments
        // ***************
        
        double moment_hX = dbm.latt->reap_moment(type, 0, dbm.species_map.at(0).at("X"));
        double moment_hY = dbm.latt->reap_moment(type, 0, dbm.species_map.at(0).at("Y"));
        double moment_bX1 = dbm.latt->reap_moment(type, 1, dbm.species_map.at(1).at("X1"));
        double moment_bY1 = dbm.latt->reap_moment(type, 1, dbm.species_map.at(1).at("Y1"));
        double moment_wXX1 = dbm.latt->reap_moment(type, 0, dbm.species_map.at(0).at("X"), 1, dbm.species_map.at(1).at("X1"));
        double moment_wYY1 = dbm.latt->reap_moment(type, 0, dbm.species_map.at(0).at("Y"), 1, dbm.species_map.at(1).at("Y1"));
        
        write_moments(dir+pad_str(timepoint,4)+"/moments/moments.txt", moment_hX, moment_hY, moment_wXX1, moment_wYY1, moment_bX1, moment_bY1);
        
        // ***************
        // MARK: - Covariances
        // ***************
        
        double tmp2, tmp3;
        double term1, term2, term3;
        
        // hX with everything else
        
        term1=0.0; term2=0.0; term3=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"));
            tmp3 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"));
            term1 += tmp2 * tmp3;
            term2 += tmp2;
            term3 += tmp3;
        };
        double cov_hX_hX = (term1 / no_markov_chains) - (term2 / no_markov_chains) * (term3 / no_markov_chains);
        
        term1=0.0; term2=0.0; term3=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"));
            tmp3 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("Y"));
            term1 += tmp2 * tmp3;
            term2 += tmp2;
            term3 += tmp3;
        };
        double cov_hX_hY = (term1 / no_markov_chains) - (term2 / no_markov_chains) * (term3 / no_markov_chains);

        term1=0.0; term2=0.0; term3=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"));
            tmp3 = dbm.latt->reap_moment_sample(type, i_chain, 1, dbm.species_map.at(1).at("X1"));
            term1 += tmp2 * tmp3;
            term2 += tmp2;
            term3 += tmp3;
        };
        double cov_hX_bX1 = (term1 / no_markov_chains) - (term2 / no_markov_chains) * (term3 / no_markov_chains);

        term1=0.0; term2=0.0; term3=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"));
            tmp3 = dbm.latt->reap_moment_sample(type, i_chain, 1, dbm.species_map.at(1).at("Y1"));
            term1 += tmp2 * tmp3;
            term2 += tmp2;
            term3 += tmp3;
        };
        double cov_hX_bY1 = (term1 / no_markov_chains) - (term2 / no_markov_chains) * (term3 / no_markov_chains);

        term1=0.0; term2=0.0; term3=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"));
            tmp3 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("X"), 1, dbm.species_map.at(1).at("X1"));
            term1 += tmp2 * tmp3;
            term2 += tmp2;
            term3 += tmp3;
        };
        double cov_hX_wXX1 = (term1 / no_markov_chains) - (term2 / no_markov_chains) * (term3 / no_markov_chains);

        term1=0.0; term2=0.0; term3=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("Y"));
            tmp3 = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_map.at(0).at("Y"), 1, dbm.species_map.at(1).at("Y1"));
            term1 += tmp2 * tmp3;
            term2 += tmp2;
            term3 += tmp3;
        };
        double cov_hX_wYY1 = (term1 / no_markov_chains) - (term2 / no_markov_chains) * (term3 / no_markov_chains);

        write_moments(dir+pad_str(timepoint,4)+"/covs/hX.txt", cov_hX_hX, cov_hX_hY, cov_hX_wXX1, cov_hX_wYY1, cov_hX_bX1, cov_hX_bY1);
    };
    
	return 0;
};

