#include <iostream>
#include <bmla>
#include <fstream>
#include <vector>
#include <map>

#include <sstream>

#include "dbm_static.hpp"

using namespace dblz;
using namespace std;

int main() {

    // ***************
    // MARK: - DBM
    // ***************
    
    DBM dbm;
    
    // ***************
    // MARK: - Set values
    // ***************
    
    dbm.hX->set_val(-2.612);
    dbm.hY->set_val(-4.654);
    dbm.bX1->set_val(-3.057);
    dbm.bY1->set_val(-3.11355);
    dbm.wXX1->set_val(1.382);
    dbm.wYY1->set_val(1.890);
    dbm.bY2->set_val(-2.635);
    dbm.bX2->set_val(-2.630);
    dbm.wX1X2->set_val(-0.151);
    dbm.wY1Y2->set_val(0.000);
    
    int no_markov_chains = 10;
    
    // Set no markov chains
    dbm.latt->set_no_markov_chains(MCType::AWAKE, no_markov_chains);
    dbm.latt->set_no_markov_chains(MCType::ASLEEP, no_markov_chains);
    
    // ***************
    // MARK: - AWAKE PHASE
    // ***************
    
    // Read in the batch
    int timepoint = 300;
    for (int i_chain=0; i_chain<no_markov_chains; i_chain++)
    {
        dbm.latt->read_layer_from_file(MCType::AWAKE, i_chain, 0, "../../stoch_sims/data/lattice_v" + pad_str(i_chain+5,3) + "/lattice/"+pad_str(timepoint,4)+".txt", true);
    };
    
    // Option (2): upward pass with 2x weights (DBM) to activate probabilitsic units
    // (faster to converge!!!)
    int no_layers = 3;
    for (auto layer=1; layer<no_layers; layer++) {
        dbm.latt->activate_layer_calculate_from_below(MCType::AWAKE, layer, 2.0);
        dbm.latt->activate_layer_convert_to_probs(MCType::AWAKE, layer, true); // binary upward bass
        dbm.latt->activate_layer_committ(MCType::AWAKE, layer);
    };
    
    // Variational inference
    int no_steps_awake = 10;
    // Sample vis, hidden
    for (int i_sampling_step=0; i_sampling_step<no_steps_awake; i_sampling_step++)
    {
        // Do odd layers
        for (auto layer=1; layer<no_layers; layer+=2) {
            if (layer != no_layers-1) {
                dbm.latt->activate_layer_calculate_from_both(MCType::AWAKE, layer);
            } else {
                dbm.latt->activate_layer_calculate_from_below(MCType::AWAKE, layer);
            };
            dbm.latt->activate_layer_convert_to_probs(MCType::AWAKE, layer, true);
            dbm.latt->activate_layer_committ(MCType::AWAKE, layer);
        };
        
        // Do even layers
        for (auto layer=2; layer<no_layers; layer+=2) {
            if (layer == no_layers-1) {
                dbm.latt->activate_layer_calculate_from_below(MCType::AWAKE, layer);
            } else {
                dbm.latt->activate_layer_calculate_from_both(MCType::AWAKE, layer);
            };
            dbm.latt->activate_layer_convert_to_probs(MCType::AWAKE, layer, true);
            dbm.latt->activate_layer_committ(MCType::AWAKE, layer);
        };
    };
    
    // Final: in parallel, use probs for hidden layers, binary for visible
    /*
    for (auto layer=1; layer<no_layers; layer++) {
        if (layer == no_layers-1) {
            dbm.latt->activate_layer_calculate_from_below(MCType::AWAKE, layer);
        } else {
            dbm.latt->activate_layer_calculate_from_both(MCType::AWAKE, layer);
        };
        dbm.latt->activate_layer_convert_to_probs(MCType::AWAKE, layer, false);
        dbm.latt->activate_layer_committ(MCType::AWAKE, layer);
    };
     */

    // Write
    for (int i_chain=0; i_chain<no_markov_chains; i_chain++) {
        dbm.latt->write_layer_to_file(MCType::AWAKE, i_chain, 0, "../data/diagnose_dbm/awake_chain_"+pad_str(i_chain, 2)+"_layer_0.txt", true);
        dbm.latt->write_layer_to_file(MCType::AWAKE, i_chain, 1, "../data/diagnose_dbm/awake_chain_"+pad_str(i_chain, 2)+"_layer_1.txt", true);
        dbm.latt->write_layer_to_file(MCType::AWAKE, i_chain, 2, "../data/diagnose_dbm/awake_chain_"+pad_str(i_chain, 2)+"_layer_2.txt", true);
    };
    
    // ***************
    // MARK: - ASLEEP PHASE
    // ***************
    
    // Reset chains (visible layer)
    for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
        for (auto layer=0; layer<no_layers; layer++) {
            dbm.latt->set_random_all_units_in_layer(MCType::ASLEEP, i_chain, layer, true);
        };
    };
    
    // Run CD sampling
    
    // Sample vis, hidden
    int no_steps_asleep = 10;
    for (int i_sampling_step=0; i_sampling_step<no_steps_asleep; i_sampling_step++)
    {
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
    /*
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
        dbm.latt->activate_layer_committ(MCType::ASLEEP, layer);
    };
     */

    // Write
    for (int i_chain=0; i_chain<no_markov_chains; i_chain++) {
        dbm.latt->write_layer_to_file(MCType::ASLEEP, i_chain, 0, "../data/diagnose_dbm/asleep_chain_"+pad_str(i_chain, 2)+"_layer_0.txt", true);
        dbm.latt->write_layer_to_file(MCType::ASLEEP, i_chain, 1, "../data/diagnose_dbm/asleep_chain_"+pad_str(i_chain, 2)+"_layer_1.txt", true);
        dbm.latt->write_layer_to_file(MCType::ASLEEP, i_chain, 2, "../data/diagnose_dbm/asleep_chain_"+pad_str(i_chain, 2)+"_layer_2.txt", true);
    };
    
	return 0;
};

