#include <iostream>
#include <dblz>
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
    
    // Centers
    map<int,map<string,CTptr>> center_map;
    map<int,vector<CTptr>> center_vec;

    // Species
    map<int,map<string,Sptr>> species_map;
    map<int,vector<Sptr>> species_vec;
    
    // Ixn funcs
    map<int,map<string,ITptr>> bias;
    map<int,map<string,ITptr>> w_ixns;
    vector<ITptr> ixns;
    
    // Scale and shift parameters
    map<int,Iptr> beta;
    map<int,Iptr> gamma;
    
    // Lattice
    shared_ptr<LatticeTrajCenteredHom> latt;
    
    // Domain
    std::vector<Domain1DParam*> dom;
    
    /********************
     Constructor
     ********************/
    
    DBM(double spacing, map<string, double> init_conds_ixns, std::map<string,double> lrs) {
        
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
        // MARK: - Centers
        // ***************
        
        cout << "--- Making centers ---" << endl;

        double init_center_val = 100.0/1600.0;
        
        sp_name = "X";
        center_map[0][sp_name] = make_shared<CenterTraj>(0,species_map.at(0).at(sp_name),init_center_val);
        center_vec[0].push_back(center_map.at(0).at(sp_name));
        sp_name = "Y";
        center_map[0][sp_name] = make_shared<CenterTraj>(0,species_map.at(0).at(sp_name),init_center_val);
        center_vec[0].push_back(center_map.at(0).at(sp_name));
        
        sp_name = "X1";
        center_map[1][sp_name] = make_shared<CenterTraj>(1,species_map.at(1).at(sp_name),init_center_val);
        center_vec[1].push_back(center_map.at(1).at(sp_name));
        sp_name = "Y1";
        center_map[1][sp_name] = make_shared<CenterTraj>(1,species_map.at(1).at(sp_name),init_center_val);
        center_vec[1].push_back(center_map.at(1).at(sp_name));

        cout << "--- [Finished] Making centers ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Ixn funcs
        // ***************

        cout << "--- Making ixn func ---" << endl;
        
        std::string ixn_name;
        
        // Layer 0
        // Ws
        ixn_name = "wXX1";
        w_ixns[0][ixn_name] = make_shared<IxnParamTraj>(ixn_name,IxnParamType::W, init_conds_ixns.at(ixn_name));
        ixns.push_back(w_ixns.at(0).at(ixn_name));
        
        ixn_name = "wYY1";
        w_ixns[0][ixn_name] = make_shared<IxnParamTraj>(ixn_name,IxnParamType::W, init_conds_ixns.at(ixn_name));
        ixns.push_back(w_ixns.at(0).at(ixn_name));

        // Bias visibles
        ixn_name = "hX";
        bias[0][ixn_name] = make_shared<IxnParamTraj>(ixn_name,IxnParamType::H, init_conds_ixns.at(ixn_name));
        ixns.push_back(bias.at(0).at(ixn_name));

        ixn_name = "hY";
        bias[0][ixn_name] = make_shared<IxnParamTraj>(ixn_name,IxnParamType::H, init_conds_ixns.at(ixn_name));
        ixns.push_back(bias.at(0).at(ixn_name));
        
        // Bias hidden
        ixn_name = "bX1";
        bias[1][ixn_name] = make_shared<IxnParamTraj>(ixn_name,IxnParamType::B, init_conds_ixns.at(ixn_name));
        ixns.push_back(bias.at(1).at(ixn_name));
        
        ixn_name = "bY1";
        bias[1][ixn_name] = make_shared<IxnParamTraj>(ixn_name,IxnParamType::B, init_conds_ixns.at(ixn_name));
        ixns.push_back(bias.at(1).at(ixn_name));

        cout << "--- [Finished] Making ixn func ---" << endl;
        cout << endl;
        
        // ***************
        // MARK: - RHS for ixn funcs
        // ***************

        cout << "--- Making diff eq rhs ---" << endl;

        double spacing_bias_v = 0.1;
        double spacing_bias_h = 0.1;
        double spacing_weight = 0.1;
        Domain1DParam *dom_hX = new Domain1DParam(bias.at(0).at("hX"), spacing_bias_v, init_conds_ixns.at("hX") - 0.5*spacing_bias_v);
        Domain1DParam *dom_hY = new Domain1DParam(bias.at(0).at("hY"), spacing_bias_v, init_conds_ixns.at("hY") - 0.5*spacing_bias_v);
        Domain1DParam *dom_wXX1 = new Domain1DParam(w_ixns.at(0).at("wXX1"), spacing_weight, init_conds_ixns.at("wXX1") - 0.5*spacing_weight);
        Domain1DParam *dom_wYY1 = new Domain1DParam(w_ixns.at(0).at("wYY1"), spacing_weight, init_conds_ixns.at("wYY1") - 0.5*spacing_weight);
        Domain1DParam *dom_bX1 = new Domain1DParam(bias.at(1).at("bX1"), spacing_bias_h, init_conds_ixns.at("bX1") - 0.5*spacing_bias_h);
        Domain1DParam *dom_bY1 = new Domain1DParam(bias.at(1).at("bY1"), spacing_bias_h, init_conds_ixns.at("bY1") - 0.5*spacing_bias_h);
        dom = std::vector<Domain1DParam*>({dom_hX,dom_hY,dom_wXX1,dom_wYY1,dom_bX1,dom_bY1});
        
        // biases
        for (auto pr1: bias) {
            for (auto pr2: pr1.second) {
                auto rhs = make_shared<DiffEqRHS>("rhs " + pr2.second->get_name(), pr2.second, dom, lrs.at(pr2.first));
                pr2.second->set_diff_eq_rhs(rhs);
            };
        };
        
        // weights
        auto rhs_wXX1 = make_shared<DiffEqRHSCenteredHomWeight>("rhs wXX1", w_ixns.at(0).at("wXX1"), dom, lrs.at("wXX1"), 4, bias.at(0).at("hX"), bias.at(1).at("bX1"), center_map.at(0).at("X"), center_map.at(1).at("X1"));
        w_ixns.at(0).at("wXX1")->set_diff_eq_rhs(rhs_wXX1);

        auto rhs_wYY1 = make_shared<DiffEqRHSCenteredHomWeight>("rhs wYY1", w_ixns.at(0).at("wYY1"), dom, lrs.at("wYY1"), 4, bias.at(0).at("hY"), bias.at(1).at("bY1"), center_map.at(0).at("Y"), center_map.at(1).at("Y1"));
        w_ixns.at(0).at("wYY1")->set_diff_eq_rhs(rhs_wYY1);

        cout << "--- [Finished] Making diff eq rhs ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Adjoint
        // ***************
        
        cout << "--- Making adjoint ---" << endl;
        
        // All biases
        std::map<int,std::map<Sptr,ITptr>> all_biases;
        all_biases[0][species_map.at(0).at("X")] = bias.at(0).at("hX");
        all_biases[0][species_map.at(0).at("Y")] = bias.at(0).at("hY");
        all_biases[1][species_map.at(1).at("X1")] = bias.at(1).at("bX1");
        all_biases[1][species_map.at(1).at("Y1")] = bias.at(1).at("bY1");
        // All weights
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> all_weights;
        all_weights[0][species_map.at(0).at("X")][1][species_map.at(1).at("X1")] = w_ixns.at(0).at("wXX1");
        all_weights[1][species_map.at(1).at("X1")][0][species_map.at(0).at("X")] = w_ixns.at(0).at("wXX1");
        all_weights[0][species_map.at(0).at("Y")][1][species_map.at(1).at("Y1")] = w_ixns.at(0).at("wYY1");
        all_weights[1][species_map.at(1).at("Y1")][0][species_map.at(0).at("Y")] = w_ixns.at(0).at("wYY1");
        // Centers
        std::map<int,std::map<Sptr,CTptr>> all_center_trajs;
        all_center_trajs[0][species_map.at(0).at("X")] = center_map.at(0).at("X");
        all_center_trajs[0][species_map.at(0).at("Y")] = center_map.at(0).at("Y");
        all_center_trajs[1][species_map.at(1).at("X1")] = center_map.at(1).at("X1");
        all_center_trajs[1][species_map.at(1).at("Y1")] = center_map.at(1).at("Y1");
        // Conn mults
        std::map<int, std::map<int,int>> conn_mults;
        conn_mults[0][1] = 4;
        conn_mults[1][0] = 4;
        
        // bias
        auto deriv_term_hX = make_shared<AdjointParamsCenteredHomDerivTerm>(bias.at(0).at("hX"), all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_hX = make_shared<AdjointParamsCenteredHomBias>("adjoint "+bias.at(0).at("hX")->get_name(),bias.at(0).at("hX"),deriv_term_hX);
        bias.at(0).at("hX")->set_adjoint(adjoint_hX);

        auto deriv_term_hY = make_shared<AdjointParamsCenteredHomDerivTerm>(bias.at(0).at("hY"), all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_hY = make_shared<AdjointParamsCenteredHomBias>("adjoint "+bias.at(0).at("hY")->get_name(),bias.at(0).at("hY"),deriv_term_hY);
        bias.at(0).at("hY")->set_adjoint(adjoint_hY);

        auto deriv_term_bX1 = make_shared<AdjointParamsCenteredHomDerivTerm>(bias.at(1).at("bX1"), all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_bX1 = make_shared<AdjointParamsCenteredHomBias>("adjoint "+bias.at(1).at("bX1")->get_name(),bias.at(1).at("bX1"),deriv_term_bX1);
        bias.at(1).at("bX1")->set_adjoint(adjoint_bX1);

        auto deriv_term_bY1 = make_shared<AdjointParamsCenteredHomDerivTerm>(bias.at(1).at("bY1"), all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_bY1 = make_shared<AdjointParamsCenteredHomBias>("adjoint "+bias.at(1).at("bY1")->get_name(),bias.at(1).at("bY1"),deriv_term_bY1);
        bias.at(1).at("bY1")->set_adjoint(adjoint_bY1);

        // weights
        auto deriv_term_wXX1 = make_shared<AdjointParamsCenteredHomDerivTerm>(w_ixns.at(0).at("wXX1"), all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_wXX1 = make_shared<AdjointParamsCenteredHomWeight>("adjoint "+w_ixns.at(0).at("wXX1")->get_name(),w_ixns.at(0).at("wXX1"),deriv_term_wXX1,deriv_term_hX,deriv_term_bX1,4,center_map.at(0).at("X"),center_map.at(1).at("X1"),adjoint_hX,adjoint_bX1,center_vec.at(0),center_vec.at(1));
        w_ixns.at(0).at("wXX1")->set_adjoint(adjoint_wXX1);

        auto deriv_term_wYY1 = make_shared<AdjointParamsCenteredHomDerivTerm>(w_ixns.at(0).at("wYY1"), all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_wYY1 = make_shared<AdjointParamsCenteredHomWeight>("adjoint "+w_ixns.at(0).at("wYY1")->get_name(),w_ixns.at(0).at("wYY1"),deriv_term_wYY1,deriv_term_hY,deriv_term_bY1,4,center_map.at(0).at("Y"),center_map.at(1).at("Y1"),adjoint_hY,adjoint_bY1,center_vec.at(0),center_vec.at(1));
        w_ixns.at(0).at("wYY1")->set_adjoint(adjoint_wYY1);
        
        cout << "--- [Finished] Making adjoint ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Lattice
        // ***************
        
        cout << "--- Making lattice ---" << endl;
        
        int no_dims = 2;
        int side_length = 40;
        latt = make_shared<LatticeTrajCenteredHom>(no_dims, side_length, species_vec.at(0), center_vec.at(0));

        cout << " > begin visible" << endl;
        
        // Add visible bias
        latt->set_bias_of_layer(0, species_map.at(0).at("X"), bias.at(0).at("hX"));
        latt->set_bias_of_layer(0, species_map.at(0).at("Y"), bias.at(0).at("hY"));

        cout << " > end visible" << endl;
        
        cout << " > begin hidden" << endl;
    
        // Add layer
        int conn_mult = 4;
        latt->add_layer(1, side_length, species_vec.at(1), center_vec.at(1), conn_mult);
        
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
