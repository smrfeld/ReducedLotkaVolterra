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
// MARK: - Read means
// ***************

// Read init conds
map<string,double> read_means(string fname) {
    map<string,double> ret;
    
    std::ifstream f;
    f.open(fname);
    if (!f.is_open()) {
        std::cerr << ">>> read_means <<< could not find file: " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string nu="";
    std::string val="";
    std::string line;
    std::istringstream iss;
    if (f.is_open()) { // make sure we found it
        while (getline(f,line)) {
            if (line == "") { continue; };
            iss = std::istringstream(line);
            iss >> nu;
            iss >> val;
            ret[nu] = (atof(val.c_str())) / 1000.0;
            nu=""; val="";
        };
    };
    f.close();
    
    return ret;
};

// ***************
// MARK: - Read init conds
// ***************

// Read init conds
map<string,double> read_ic(string fname) {
    map<string,double> ret;
    
    std::ifstream f;
    f.open(fname);
    if (!f.is_open()) {
        std::cerr << ">>> read_ic <<< could not find file: " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string nu="";
    std::string val="";
    std::string line;
    std::istringstream iss;
    if (f.is_open()) { // make sure we found it
        while (getline(f,line)) {
            if (line == "") { continue; };
            iss = std::istringstream(line);
            iss >> nu;
            iss >> val;
            ret[nu] = (atof(val.c_str()));
            nu=""; val="";
        };
    };
    f.close();
    
    return ret;
};


int main() {

    // ***************
    // MARK: - DBM
    // ***************
    
    DBM dbm;

    // No markov chains
    int no_markov_chains = 10;
    
    // Set no markov chains
    dbm.latt->set_no_markov_chains(MCType::AWAKE, no_markov_chains);
    dbm.latt->set_no_markov_chains(MCType::ASLEEP, no_markov_chains);

    // ***************
    // MARK: - Go over timepoints
    // ***************
    
    for (auto timepoint=160; timepoint<=165; timepoint++) {
        
        // ***************
        // MARK: - Set values
        // ***************

        switch(timepoint) {
            case 157:   dbm.bias.at(0).at("hX")->set_val(-2.06578);
                        dbm.bias.at(0).at("hY")->set_val(-1.71815);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.28006);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.212472);
                        dbm.bias.at(1).at("bX1")->set_val(-2.69083);
                        dbm.bias.at(1).at("bY1")->set_val(-2.61613);
                        break;
            case 159:   dbm.bias.at(0).at("hX")->set_val(-2.15172);
                        dbm.bias.at(0).at("hY")->set_val(-1.65791);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.32847);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.267789);
                        dbm.bias.at(1).at("bX1")->set_val(-2.7325);
                        dbm.bias.at(1).at("bY1")->set_val(-2.60391);
                        break;
            case 160:   dbm.bias.at(0).at("hX")->set_val(-2.17778);
                        dbm.bias.at(0).at("hY")->set_val(-1.64267);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.25453);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.27609);
                        dbm.bias.at(1).at("bX1")->set_val(-2.755);
                        dbm.bias.at(1).at("bY1")->set_val(-2.60096);
                        break;
            case 161:   dbm.bias.at(0).at("hX")->set_val(-2.16443);
                        dbm.bias.at(0).at("hY")->set_val(-1.62611);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.37291);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.289772);
                        dbm.bias.at(1).at("bX1")->set_val(-2.73574);
                        dbm.bias.at(1).at("bY1")->set_val(-2.60009);
                        break;
            case 162:   dbm.bias.at(0).at("hX")->set_val(-2.31274);
                        dbm.bias.at(0).at("hY")->set_val(-1.51821);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.18779);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.449188);
                        dbm.bias.at(1).at("bX1")->set_val(-2.82509);
                        dbm.bias.at(1).at("bY1")->set_val(-2.56256);
                        break;
            case 163:   dbm.bias.at(0).at("hX")->set_val(-2.29682);
                        dbm.bias.at(0).at("hY")->set_val(-1.51671);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.28987);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.461686);
                        dbm.bias.at(1).at("bX1")->set_val(-2.80688);
                        dbm.bias.at(1).at("bY1")->set_val(-2.56156);
                        break;
            case 164:   dbm.bias.at(0).at("hX")->set_val(-2.29461);
                        dbm.bias.at(0).at("hY")->set_val(-1.5038);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.41886);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.500612);
                        dbm.bias.at(1).at("bX1")->set_val(-2.79542);
                        dbm.bias.at(1).at("bY1")->set_val(-2.55449);
                        break;
            case 165:   dbm.bias.at(0).at("hX")->set_val(-2.35213);
                        dbm.bias.at(0).at("hY")->set_val(-1.48572);
                        dbm.w_ixns.at(0).at("wXX1")->set_val(1.26528);
                        dbm.w_ixns.at(0).at("wYY1")->set_val(0.532876);
                        dbm.bias.at(1).at("bX1")->set_val(-2.84515);
                        dbm.bias.at(1).at("bY1")->set_val(-2.54355);
        };
        
        // ***************
        // MARK: - CD sampling func
        // ***************
        
        /*
        std::vector<FName> fnames;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            fnames.push_back(FName("../../stoch_sims/data/lattice_v"+pad_str(i_chain+4,3)+"/lattice/"+pad_str(timepoint,4)+".txt", true));
        };
        int no_cd_steps = 50;
        dbm.latt->wake_sleep_loop_rbm_cd(0, no_cd_steps, fnames, OptionsWakeSleep_RBM_CD());
        
        // Write
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(MCType::ASLEEP, i_chain, 0, "../data/cd_sampling/"+pad_str(timepoint,4)+"_layer_0_"+pad_str(i_chain,3)+"_"+pad_str(no_cd_steps,2)+".txt", true);
        };
         */
        
        // ***************
        // MARK: - Read initial lattice
        // ***************
        
        /*
        int layer = 0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->read_layer_from_file(MCType::ASLEEP, i_chain, layer, "../../stoch_sims/data/lattice_v"+pad_str(i_chain+4,3)+"/lattice/"+pad_str(timepoint,4)+".txt", true);
        };
         */
        
        // ***************
        // MARK: - Sample
        // ***************

        // Random initial
        int layer = 0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->set_random_all_units_in_layer(MCType::ASLEEP, i_chain, layer, true);
        };
        
        MCType type = MCType::ASLEEP;
        
        // Write init
        layer = 0;
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(type, i_chain, layer, "../data/cd_sampling/"+pad_str(timepoint,4)+"_layer_0_"+pad_str(i_chain,3)+"_00.txt", true);
        };
        
        int no_steps = 50;
        for (auto i_step=0; i_step<no_steps; i_step++) {
            
            // Up
            layer = 1;
            dbm.latt->activate_layer_calculate_from_below(type, layer);
            dbm.latt->activate_layer_convert_to_probs(type, layer, true);
            dbm.latt->activate_layer_committ(type, layer);
            // Write
            for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
                dbm.latt->write_layer_to_file(type, i_chain, layer, "../data/cd_sampling/"+pad_str(timepoint,4)+"_layer_1_"+pad_str(i_chain,3)+"_"+pad_str(i_step,2)+".txt", true);
            };

            // Down
            layer = 0;
            dbm.latt->activate_layer_calculate_from_above(type, layer);
            dbm.latt->activate_layer_convert_to_probs(type, layer, true);
            dbm.latt->activate_layer_committ(type, layer);
            // Write
            for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
                dbm.latt->write_layer_to_file(type, i_chain, layer, "../data/cd_sampling/"+pad_str(timepoint,4)+"_layer_0_"+pad_str(i_chain,3)+"_"+pad_str(i_step+1,2)+".txt", true);
            };

        };
        
        // Up with prob
        layer = 1;
        dbm.latt->activate_layer_calculate_from_below(type, layer);
        dbm.latt->activate_layer_convert_to_probs(type, layer, false);
        dbm.latt->activate_layer_committ(type, layer);
        // Write
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            dbm.latt->write_layer_to_file(type, i_chain, layer, "../data/cd_sampling/"+pad_str(timepoint,4)+"_layer_1_"+pad_str(i_chain,3)+"_"+pad_str(no_steps,2)+".txt", false);
        };
    };
    
	return 0;
};

