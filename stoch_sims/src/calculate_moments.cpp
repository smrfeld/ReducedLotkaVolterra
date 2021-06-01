#include <iostream>
#include <bmla>
#include <fstream>
#include <vector>
#include <map>

#include <sstream>

using namespace dblz;
using namespace std;

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

int main() {

    // ***************
    // MARK: - Lattice
    // ***************
    
    auto species_X = make_shared<Species>("X");
    auto species_Y = make_shared<Species>("Y");

    int no_dims = 2;
    int side_length = 40;
    auto latt = make_shared<Lattice>(no_dims, side_length, std::vector<Sptr>({species_X,species_Y}));

    // ***************
    // MARK: - Setup optimizer
    // ***************
    
    // No markov chains
    int no_markov_chains = 100;
    
    latt->set_no_markov_chains(MCType::AWAKE, no_markov_chains);
    latt->set_no_markov_chains(MCType::ASLEEP, no_markov_chains);

    // ***************
    // MARK: - Go through all times
    // ***************
    
    MCType type = MCType::ASLEEP;
    int no_timesteps=500;
    
    for (auto timepoint=0; timepoint<no_timesteps; timepoint++) {
        
        std::cout << timepoint << " / " << no_timesteps << std::endl;
        
        // Read
        for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
            latt->read_layer_from_file(type, i_chain, 0, "../data/lattice_v"+pad_str(i_chain+1,3)+"/lattice/"+pad_str(timepoint,4)+".txt", true);
        };
        
        // ***************
        // MARK: - Reap moments
        // ***************

        double moment_nn_XX = latt->reap_moment_nn(type, species_X, species_X);
        double moment_nn_YY = latt->reap_moment_nn(type, species_Y, species_Y);
        double moment_nn_XY = latt->reap_moment_nn(type, species_X, species_Y);
        double moment_nnn_XX = latt->reap_moment_nnn(type, species_X, species_X);
        double moment_nnn_YY = latt->reap_moment_nnn(type, species_Y, species_Y);
        double moment_nnn_XY = latt->reap_moment_nnn(type, species_X, species_Y);

        write_moments_other("moments_other/"+pad_str(timepoint,4)+".txt", moment_nn_XX, moment_nn_YY, moment_nn_XY, moment_nnn_XX, moment_nnn_YY, moment_nnn_XY);
    };
    
	return 0;
};

