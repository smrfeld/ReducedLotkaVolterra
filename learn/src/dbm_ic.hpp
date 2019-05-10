#include <iostream>
#include <bmla>
#include <fstream>
#include <vector>
#include <map>

#include <sstream>

using namespace dblz;
using namespace std;

// ***************
// MARK: - DBM
// ***************

struct DBM {
    
    /********************
     Members
     ********************/
    
    // Species
    map<int,map<string,shared_ptr<Species>>> species_map;
    map<int,vector<shared_ptr<Species>>> species_vec;
    
    // Ixn funcs
    map<int,map<string,shared_ptr<IxnParam>>> bias;
    map<int,map<string,shared_ptr<IxnParam>>> w_ixns;
    vector<shared_ptr<IxnParam>> ixns;
    
    // Scale and shift parameters
    map<int,Iptr> beta;
    map<int,Iptr> gamma;
    
    // Lattice
    shared_ptr<Lattice> latt;
    
    /********************
     Constructor
     ********************/
    
    DBM() {
        
        // ***************
        // MARK: - Species
        // ***************
        
        cout << "--- Making species ---" << endl;
        
        string sp_name;
        
        // 0
        sp_name = "X";
        species_map[0][sp_name] = make_shared<Species>(sp_name);
        species_vec[0].push_back(species_map[0][sp_name]);
        sp_name = "Y";
        species_map[0][sp_name] = make_shared<Species>(sp_name);
        species_vec[0].push_back(species_map[0][sp_name]);
        
        sp_name = "X1";
        species_map[1][sp_name] = make_shared<Species>(sp_name);
        species_vec[1].push_back(species_map[1][sp_name]);
        sp_name = "Y1";
        species_map[1][sp_name] = make_shared<Species>(sp_name);
        species_vec[1].push_back(species_map[1][sp_name]);
        
        cout << "--- [Finished] Making species ---" << endl;
        cout << endl;
        
        // ***************
        // MARK: - Ixn funcs
        // ***************
        
        cout << "--- Making ixn func ---" << endl;
        
        std::string ixn_name;
        
        // Layer 0
        // Ws
        ixn_name = "wXX1";
        w_ixns[0][ixn_name] = make_shared<IxnParam>(ixn_name,IxnParamType::W, 0.0, 0.0);
        ixns.push_back(w_ixns.at(0).at(ixn_name));
        
        ixn_name = "wYY1";
        w_ixns[0][ixn_name] = make_shared<IxnParam>(ixn_name,IxnParamType::W, 0.0, 0.0);
        ixns.push_back(w_ixns.at(0).at(ixn_name));
        
        // Bias visibles
        ixn_name = "hX";
        bias[0][ixn_name] = make_shared<IxnParam>(ixn_name,IxnParamType::H, 0.0, 0.0);
        ixns.push_back(bias.at(0).at(ixn_name));
        
        ixn_name = "hY";
        bias[0][ixn_name] = make_shared<IxnParam>(ixn_name,IxnParamType::H, 0.0, 0.0);
        ixns.push_back(bias.at(0).at(ixn_name));
        
        // Bias hidden
        ixn_name = "bX1";
        bias[1][ixn_name] = make_shared<IxnParam>(ixn_name,IxnParamType::B, 0.0, 0.0);
        ixns.push_back(bias.at(1).at(ixn_name));
        
        ixn_name = "bY1";
        bias[1][ixn_name] = make_shared<IxnParam>(ixn_name,IxnParamType::B, 0.0, 0.0);
        ixns.push_back(bias.at(1).at(ixn_name));
        
        cout << "--- [Finished] Making ixn func ---" << endl;
        cout << endl;
        
        // ***************
        // MARK: - Lattice
        // ***************
        
        cout << "--- Making lattice ---" << endl;
        
        int no_dims = 2;
        int side_length = 40;
        latt = make_shared<Lattice>(no_dims,side_length,species_vec.at(0));
        
        cout << " > begin visible" << endl;
        
        // Add visible bias
        latt->set_bias_of_layer(0, species_map.at(0).at("X"), bias.at(0).at("hX"));
        latt->set_bias_of_layer(0, species_map.at(0).at("Y"), bias.at(0).at("hY"));
        
        cout << " > end visible" << endl;
        
        cout << " > begin hidden" << endl;
        
        // Add layer
        latt->add_layer(1, side_length, species_vec[1]);
        
        // Connectivity
        for (auto i=1; i<=side_length; i++) {
            for (auto j=1; j<=side_length; j++) {
                
                // Displacements
                for (auto i2=0; i2<=1; i2++) {
                    for (auto j2=0; j2<=1; j2++) {
                        auto i3 = i+i2;
                        if (i3 > side_length) {
                            i3 = i3-side_length;
                        };
                        if (i3 < 1) {
                            i3 = side_length + i3;
                        };
                        
                        auto j3 = j+j2;
                        if (j3 > side_length) {
                            j3 = j3-side_length;
                        };
                        if (j3 < 1) {
                            j3 = side_length + j3;
                        };
                        
                        latt->add_conn(1-1, i, j, 1, i3, j3);
                    };
                };
            };
        };
        
        // Biases
        latt->set_bias_of_layer(1, species_map.at(1).at("X1"), bias.at(1).at("bX1"));
        latt->set_bias_of_layer(1, species_map.at(1).at("Y1"), bias.at(1).at("bY1"));
        
        // Ixns
        latt->set_ixn_between_layers(0, species_map.at(0).at("X"), 1, species_map.at(1).at("X1"),w_ixns.at(0).at("wXX1"));
        latt->set_ixn_between_layers(0, species_map.at(0).at("Y"), 1, species_map.at(1).at("Y1"),w_ixns.at(0).at("wYY1"));
        
        cout << " > end hidden" << endl;
        
        cout << "--- [Finished] Making lattice ---" << endl;
        cout << endl;
        
    };
};

// ***************
// MARK: - Write function
// ***************

void write(std::map<string,double> aves, string fname) {
    std::ofstream f;
    
    // Open
    f.open(fname);
    
    // Make sure we found it
    if (!f.is_open()) {
        std::cerr << ">>> Error: write <<< could not write to file: " << fname << std::endl;
        exit(EXIT_FAILURE);
    };
    
    // Go through all
    int i=0;
    for (auto pr: aves) {
        f << pr.first << " " << pr.second;
        if (i != aves.size()-1) {
            f << "\n";
        };
        i++;
    };
    
    // Close
    f.close();
};
