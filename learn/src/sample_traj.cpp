#include <iostream>
#include <bmla>
#include <fstream>
#include <vector>
#include <map>

#include <sstream>

#include "dbm_static.hpp"

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

void write_moments(string fname, double moment_hX, double moment_hY, double moment_wXX1, double moment_wYY1, double moment_bX1, double moment_bY1, double moment_wX1X2, double moment_wY1Y2, double moment_bX2, double moment_bY2) {
    
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
    f << "wX1X2 " << moment_wX1X2 << "\n";
    f << "wY1Y2 " << moment_wY1Y2 << "\n";
    f << "bX2 " << moment_bX2 << "\n";
    f << "bY2 " << moment_bY2;

    f.close();
};

// ***************
// MARK: - Write moments
// ***************

void write_moments_other(string fname, double moment_nn_XX, double moment_nn_YY, double moment_nn_XY, double moment_nnn_XX, double moment_nnn_YY, double moment_nnn_XY) {
    
    std::ofstream f;
    f.open (fname);
    if (!f.is_open()) { // make sure we found it
        std::cerr << ">>> write_moments_other <<< Error: could not open file: " << fname << " for writing" << std::endl;
        exit(EXIT_FAILURE);
    };
    
    f << "NN XX " << moment_nn_XX << "\n";
    f << "NN YY " << moment_nn_YY << "\n";
    f << "NN XY " << moment_nn_XY << "\n";
    f << "NNN XX " << moment_nnn_XX << "\n";
    f << "NNN YY " << moment_nnn_YY << "\n";
    f << "NNN XY " << moment_nnn_XY;
    
    f.close();
};

// ***************
// MARK: - Get covariances
// ***************

map<string,double> get_covariances(MCType type, std::map<int,double> tmp2, double term2, int no_markov_chains, shared_ptr<Lattice> latt, Sptr species_X, Sptr species_Y, Sptr species_X1, Sptr species_Y1, Sptr species_X2, Sptr species_Y2) {

    map<string,double> ret;
    double term1,term3, tmp3;
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 0, species_X);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["hX"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 0, species_Y);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["hY"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 1, species_X1);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["bX1"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 1, species_Y1);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["bY1"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 0, species_X, 1, species_X1);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["wXX1"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 0, species_Y, 1, species_Y1);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["wYY1"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 2, species_X2);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["bX2"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 2, species_Y2);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["bY2"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 1, species_X1, 2, species_X2);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["wX1X2"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    term1=0.0; term3=0.0;
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        tmp3 = latt->reap_moment_sample(type, i_chain, 1, species_Y1, 2, species_Y2);
        term1 += tmp2.at(i_chain) * tmp3;
        term3 += tmp3;
    };
    ret["wY1Y2"] = (term1 / no_markov_chains) - (term2) * (term3 / no_markov_chains);
    
    return ret;
};

int main() {

    // ***************
    // MARK: - DBM
    // ***************
    
    auto dbm = DBM(4,4);

    // ***************
    // MARK: - Setup optimizer
    // ***************
    
    // No markov chains
    int no_markov_chains = 100;
    
    dbm.latt->set_no_markov_chains(MCType::AWAKE, no_markov_chains);
    dbm.latt->set_no_markov_chains(MCType::ASLEEP, no_markov_chains);

    // ***************
    // MARK: - Read
    // ***************
    
    std::string dir_read = "learn_centered";
    std::string dir = "../data/"+dir_read+"/sample_traj_lp/";
    std::string dir_read_ixns = "ixn_params_lp";
    
    int timepoint_read = 980;
    auto traj_hX = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/hX_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_hY = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/hY_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wXX1 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/wXX1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wYY1 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/wYY1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bX1 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/bX1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bY1 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/bY1_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wX1X2 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/wX1X2_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_wY1Y2 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/wY1Y2_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bX2 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/bX2_"+pad_str(timepoint_read, 5)+".txt");
    auto traj_bY2 = read_ixn_traj("../data/"+dir_read+"/"+dir_read_ixns+"/bY2_"+pad_str(timepoint_read, 5)+".txt");

    // ***************
    // MARK: - Go through all times
    // ***************
    
    int no_layers = 3;
    int no_sampling_steps = 10;
    MCType type = MCType::ASLEEP;
    int no_timesteps=500;
    
    for (auto timepoint=370; timepoint<=370; timepoint++) {
    //for (auto timepoint=0; timepoint<no_timesteps; timepoint++) {

        std::cout << timepoint << " / " << no_timesteps << std::endl;
        
        // Set vals
        dbm.hX->set_val(traj_hX.at(timepoint));
        dbm.hY->set_val(traj_hY.at(timepoint));
        dbm.wXX1->set_val(traj_wXX1.at(timepoint));
        dbm.wYY1->set_val(traj_wYY1.at(timepoint));
        dbm.bX1->set_val(traj_bX1.at(timepoint));
        dbm.bY1->set_val(traj_bY1.at(timepoint));
        dbm.wX1X2->set_val(traj_wX1X2.at(timepoint));
        dbm.wY1Y2->set_val(traj_wY1Y2.at(timepoint));
        dbm.bX2->set_val(traj_bX2.at(timepoint));
        dbm.bY2->set_val(traj_bY2.at(timepoint));

        // Random initial
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            for (auto layer=0; layer<no_layers; layer++) {
                dbm.latt->set_random_all_units_in_layer(MCType::ASLEEP, i_chain, layer, true);
            };
        };
        
        // Sample
        for (auto i_step=0; i_step<no_sampling_steps-1; i_step++) {
            
            // Do odd layers
            for (auto layer=1; layer<no_layers; layer+=2) {
                if (layer != no_layers-1) {
                    dbm.latt->activate_layer_calculate_from_both(MCType::ASLEEP, layer);
                } else {
                    dbm.latt->activate_layer_calculate_from_below(MCType::ASLEEP, layer);
                };
                dbm.latt->activate_layer_convert_to_probs(MCType::ASLEEP, layer, true);
                dbm.latt->activate_layer_committ(MCType::ASLEEP, layer);
            };
            
            // Do even layers
            for (auto layer=0; layer<no_layers; layer+=2) {
                if (layer == 0) {
                    dbm.latt->activate_layer_calculate_from_above(MCType::ASLEEP, layer);
                } else if (layer == no_layers-1) {
                    dbm.latt->activate_layer_calculate_from_below(MCType::ASLEEP, layer);
                } else {
                    dbm.latt->activate_layer_calculate_from_both(MCType::ASLEEP, layer);
                };
                dbm.latt->activate_layer_convert_to_probs(MCType::ASLEEP, layer, true);
                dbm.latt->activate_layer_committ(MCType::ASLEEP, layer);
            };
        };
        
        // Final: in parallel, use probs for hidden layers, binary for visible
        for (auto layer=0; layer<no_layers; layer++) {
            if (layer == 0) {
                dbm.latt->activate_layer_calculate_from_above(MCType::ASLEEP, layer);
            } else if (layer == no_layers-1) {
                dbm.latt->activate_layer_calculate_from_below(MCType::ASLEEP, layer);
            } else {
                dbm.latt->activate_layer_calculate_from_both(MCType::ASLEEP, layer);
            };
            if (layer == 0) {
                dbm.latt->activate_layer_convert_to_probs(MCType::ASLEEP, layer, false);
            } else {
                dbm.latt->activate_layer_convert_to_probs(MCType::ASLEEP, layer, false);
            };
        };
        for (auto layer=0; layer<no_layers; layer++) {
            dbm.latt->activate_layer_committ(MCType::ASLEEP, layer);
        };
        
        // Write
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(type, i_chain, 0, dir+pad_str(timepoint, 4)+"/lattices/layer_0_"+pad_str(i_chain, 4)+".txt", false);
            dbm.latt->write_layer_to_file(type, i_chain, 1, dir+pad_str(timepoint, 4)+"/lattices/layer_1_"+pad_str(i_chain, 4)+".txt", false);
            dbm.latt->write_layer_to_file(type, i_chain, 2, dir+pad_str(timepoint, 4)+"/lattices/layer_2_"+pad_str(i_chain, 4)+".txt", false);
        };

        // ***************
        // MARK: - MOMENT: NUMBER OF X
        // ***************
        
        // ***************
        // MARK: - Reap moments
        // ***************
        
        double moment_hX = dbm.latt->reap_moment(type, 0, dbm.species_X);
        double moment_hY = dbm.latt->reap_moment(type, 0, dbm.species_Y);
        double moment_bX1 = dbm.latt->reap_moment(type, 1, dbm.species_X1);
        double moment_bY1 = dbm.latt->reap_moment(type, 1, dbm.species_Y1);
        double moment_wXX1 = dbm.latt->reap_moment(type, 0, dbm.species_X, 1, dbm.species_X1);
        double moment_wYY1 = dbm.latt->reap_moment(type, 0, dbm.species_Y, 1, dbm.species_Y1);
        double moment_bX2 = dbm.latt->reap_moment(type, 2, dbm.species_X2);
        double moment_bY2 = dbm.latt->reap_moment(type, 2, dbm.species_Y2);
        double moment_wX1X2 = dbm.latt->reap_moment(type, 1, dbm.species_X1, 2, dbm.species_X2);
        double moment_wY1Y2 = dbm.latt->reap_moment(type, 1, dbm.species_Y1, 2, dbm.species_Y2);

        write_moments(dir+pad_str(timepoint,4)+"/moments/moments.txt", moment_hX, moment_hY, moment_wXX1, moment_wYY1, moment_bX1, moment_bY1, moment_wX1X2, moment_wY1Y2, moment_bX2, moment_bY2);
        
        // ***************
        // MARK: - Covariances
        // ***************
        
        // Sample values
        std::map<int,double> tmp2;
        double term2=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2[i_chain] = dbm.latt->reap_moment_sample(type, i_chain, 0, dbm.species_X);
            term2 += tmp2.at(i_chain);
        };
        term2 /= no_markov_chains;
                
        auto cov_map = get_covariances(type, tmp2, term2, no_markov_chains, dbm.latt, dbm.species_X, dbm.species_Y, dbm.species_X1, dbm.species_Y1, dbm.species_X2, dbm.species_Y2);
        
        write_moments(dir+pad_str(timepoint,4)+"/covs/hX.txt", cov_map.at("hX"), cov_map.at("hY"), cov_map.at("wXX1"), cov_map.at("wYY1"), cov_map.at("bX1"), cov_map.at("bY1"), cov_map.at("wX1X2"), cov_map.at("wY1Y2"), cov_map.at("bX2"), cov_map.at("bY2"));
        
        // ***************
        // MARK: - MOMENT: NNs and NNNs
        // ***************
        
        // ***************
        // MARK: - Reap moments
        // ***************
        
        double moment_nn_XX = dbm.latt->reap_moment_nn(type, dbm.species_X, dbm.species_X);
        double moment_nn_YY = dbm.latt->reap_moment_nn(type, dbm.species_Y, dbm.species_Y);
        double moment_nn_XY = dbm.latt->reap_moment_nn(type, dbm.species_X, dbm.species_Y);
        double moment_nnn_XX = dbm.latt->reap_moment_nnn(type, dbm.species_X, dbm.species_X);
        double moment_nnn_YY = dbm.latt->reap_moment_nnn(type, dbm.species_Y, dbm.species_Y);
        double moment_nnn_XY = dbm.latt->reap_moment_nnn(type, dbm.species_X, dbm.species_Y);
        
        write_moments_other(dir+pad_str(timepoint,4)+"/moments/moments_other.txt", moment_nn_XX, moment_nn_YY, moment_nn_XY, moment_nnn_XX, moment_nnn_YY, moment_nnn_XY);

        // ***************
        // MARK: - Covariances
        // ***************
        
        // moment_nn_XX with everything else
        
        tmp2.clear();
        term2=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2[i_chain] = dbm.latt->reap_moment_nn_sample(type, i_chain, dbm.species_X, dbm.species_X);
            term2 += tmp2.at(i_chain);
        };
        term2 /= no_markov_chains;
        
        cov_map = get_covariances(type, tmp2, term2, no_markov_chains, dbm.latt, dbm.species_X, dbm.species_Y, dbm.species_X1, dbm.species_Y1, dbm.species_X2, dbm.species_Y2);
        
        write_moments(dir+pad_str(timepoint,4)+"/covs/nn_XX.txt", cov_map.at("hX"), cov_map.at("hY"), cov_map.at("wXX1"), cov_map.at("wYY1"), cov_map.at("bX1"), cov_map.at("bY1"), cov_map.at("wX1X2"), cov_map.at("wY1Y2"), cov_map.at("bX2"), cov_map.at("bY2"));
        
        // moment_nnn_XY with everything else
        
        tmp2.clear();
        term2=0.0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            tmp2[i_chain] = dbm.latt->reap_moment_nnn_sample(type, i_chain, dbm.species_X, dbm.species_Y);
            term2 += tmp2.at(i_chain);
        };
        term2 /= no_markov_chains;
        
        cov_map = get_covariances(type, tmp2, term2, no_markov_chains, dbm.latt, dbm.species_X, dbm.species_Y, dbm.species_X1, dbm.species_Y1, dbm.species_X2, dbm.species_Y2);
        
        write_moments(dir+pad_str(timepoint,4)+"/covs/nnn_XY.txt", cov_map.at("hX"), cov_map.at("hY"), cov_map.at("wXX1"), cov_map.at("wYY1"), cov_map.at("bX1"), cov_map.at("bY1"), cov_map.at("wX1X2"), cov_map.at("wY1Y2"), cov_map.at("bX2"), cov_map.at("bY2"));
    };
    
	return 0;
};

