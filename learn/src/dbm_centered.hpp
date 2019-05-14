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
    CTptr center_X, center_Y, center_X1, center_Y1, center_X2, center_Y2;
    vector<CTptr> centers;

    // Species
    Sptr species_X, species_Y, species_X1, species_Y1, species_X2, species_Y2;
    
    // Ixn funcs
    ITptr hX, hY, bX1, bY1, bX2, bY2, wXX1, wYY1, wX1X2, wY1Y2;
    vector<ITptr> ixns;
    
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
        
        species_X = make_shared<Species>("X");
        species_Y = make_shared<Species>("Y");
        species_X1 = make_shared<Species>("X1");
        species_Y1 = make_shared<Species>("Y1");
        species_X2 = make_shared<Species>("X2");
        species_Y2 = make_shared<Species>("Y2");

        cout << "--- [Finished] Making species ---" << endl;
        cout << endl;
        
        // ***************
        // MARK: - Centers
        // ***************
        
        cout << "--- Making centers ---" << endl;

        double init_center_val = 100.0/1600.0;
        
        center_X = make_shared<CenterTraj>(0,species_X,init_center_val);
        center_Y = make_shared<CenterTraj>(0,species_Y,init_center_val);
        center_X1 = make_shared<CenterTraj>(1,species_X1,init_center_val);
        center_Y1 = make_shared<CenterTraj>(1,species_Y1,init_center_val);
        center_X2 = make_shared<CenterTraj>(2,species_X2,init_center_val);
        center_Y2 = make_shared<CenterTraj>(2,species_Y2,init_center_val);
        centers = std::vector<CTptr>({center_X,center_Y,center_X1,center_Y1,center_X2,center_Y2});
        
        cout << "--- [Finished] Making centers ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Ixn funcs
        // ***************

        cout << "--- Making ixn func ---" << endl;
        
        hX = make_shared<IxnParamTraj>("hX", IxnParamType::H, init_conds_ixns.at("hX"));
        hY = make_shared<IxnParamTraj>("hY", IxnParamType::H, init_conds_ixns.at("hY"));
        bX1 = make_shared<IxnParamTraj>("bX1", IxnParamType::B, init_conds_ixns.at("bX1"));
        bY1 = make_shared<IxnParamTraj>("bY1", IxnParamType::B, init_conds_ixns.at("bY1"));
        bX2 = make_shared<IxnParamTraj>("bX2", IxnParamType::B, init_conds_ixns.at("bX2"));
        bY2 = make_shared<IxnParamTraj>("bY2", IxnParamType::B, init_conds_ixns.at("bY2"));
        wXX1 = make_shared<IxnParamTraj>("wXX1", IxnParamType::W, init_conds_ixns.at("wXX1"));
        wYY1 = make_shared<IxnParamTraj>("wYY1", IxnParamType::W, init_conds_ixns.at("wYY1"));
        wX1X2 = make_shared<IxnParamTraj>("wX1X2", IxnParamType::X, init_conds_ixns.at("wX1X2"));
        wY1Y2 = make_shared<IxnParamTraj>("wY1Y2", IxnParamType::X, init_conds_ixns.at("wY1Y2"));
        ixns = std::vector<ITptr>({hX,hY,bX1,bY1,bX2,bY2,wXX1,wYY1,wX1X2,wY1Y2});
        
        cout << "--- [Finished] Making ixn func ---" << endl;
        cout << endl;
        
        // ***************
        // MARK: - RHS for ixn funcs
        // ***************

        cout << "--- Making diff eq rhs ---" << endl;

        double spacing_bias_v = spacing;
        double spacing_bias_h = spacing;
        double spacing_weight = spacing;
        Domain1DParam *dom_hX = new Domain1DParam(hX, spacing_bias_v, init_conds_ixns.at("hX") - 0.5*spacing_bias_v);
        Domain1DParam *dom_hY = new Domain1DParam(hY, spacing_bias_v, init_conds_ixns.at("hY") - 0.5*spacing_bias_v);
        Domain1DParam *dom_wXX1 = new Domain1DParam(wXX1, spacing_weight, init_conds_ixns.at("wXX1") - 0.5*spacing_weight);
        Domain1DParam *dom_wYY1 = new Domain1DParam(wYY1, spacing_weight, init_conds_ixns.at("wYY1") - 0.5*spacing_weight);
        //Domain1DParam *dom_bX1 = new Domain1DParam(bX1, spacing_bias_h, init_conds_ixns.at("bX1") - 0.5*spacing_bias_h);
        //Domain1DParam *dom_bY1 = new Domain1DParam(bY1, spacing_bias_h, init_conds_ixns.at("bY1") - 0.5*spacing_bias_h);
        Domain1DParam *dom_wX1X2 = new Domain1DParam(wX1X2, spacing_weight, init_conds_ixns.at("wX1X2") - 0.5*spacing_weight);
        Domain1DParam *dom_wY1Y2 = new Domain1DParam(wY1Y2, spacing_weight, init_conds_ixns.at("wY1Y2") - 0.5*spacing_weight);
        //Domain1DParam *dom_bX2 = new Domain1DParam(bX2, spacing_bias_h, init_conds_ixns.at("bX2") - 0.5*spacing_bias_h);
        //Domain1DParam *dom_bY2 = new Domain1DParam(bY2, spacing_bias_h, init_conds_ixns.at("bY2") - 0.5*spacing_bias_h);
        dom = std::vector<Domain1DParam*>({dom_hX,dom_hY,dom_wXX1,dom_wYY1,dom_wX1X2,dom_wY1Y2});
        // dom = std::vector<Domain1DParam*>({dom_wXX1,dom_wYY1,dom_wX1X2,dom_wY1Y2});
        
        // biases
        auto rhs_hX = make_shared<DiffEqRHS>("rhs hX", hX, dom, lrs.at("hX"));
        hX->set_diff_eq_rhs(rhs_hX);
        
        auto rhs_hY = make_shared<DiffEqRHS>("rhs hY", hY, dom, lrs.at("hY"));
        hY->set_diff_eq_rhs(rhs_hY);

        auto rhs_bX1 = make_shared<DiffEqRHS>("rhs bX1", bX1, dom, lrs.at("bX1"));
        bX1->set_diff_eq_rhs(rhs_bX1);

        auto rhs_bY1 = make_shared<DiffEqRHS>("rhs bY1", bY1, dom, lrs.at("bY1"));
        bY1->set_diff_eq_rhs(rhs_bY1);

        auto rhs_bX2 = make_shared<DiffEqRHS>("rhs bX2", bX2, dom, lrs.at("bX2"));
        bX2->set_diff_eq_rhs(rhs_bX2);

        auto rhs_bY2 = make_shared<DiffEqRHS>("rhs bY2", bY2, dom, lrs.at("bY2"));
        bY2->set_diff_eq_rhs(rhs_bY2);
        
        // weights
        auto rhs_wXX1 = make_shared<DiffEqRHSCenteredHomWeight>("rhs wXX1", wXX1, dom, lrs.at("wXX1"), 4, hX, bX1, center_X, center_X1);
        wXX1->set_diff_eq_rhs(rhs_wXX1);

        auto rhs_wYY1 = make_shared<DiffEqRHSCenteredHomWeight>("rhs wYY1", wYY1, dom, lrs.at("wYY1"), 4, hY, bY1, center_Y, center_Y1);
        wYY1->set_diff_eq_rhs(rhs_wYY1);

        auto rhs_wX1X2 = make_shared<DiffEqRHSCenteredHomWeight>("rhs wX1X2", wX1X2, dom, lrs.at("wX1X2"), 4, bX1, bX2, center_X1, center_X2);
        wX1X2->set_diff_eq_rhs(rhs_wX1X2);

        auto rhs_wY1Y2 = make_shared<DiffEqRHSCenteredHomWeight>("rhs wY1Y2", wY1Y2, dom, lrs.at("wY1Y2"), 4, bY1, bY2, center_Y1, center_Y2);
        wY1Y2->set_diff_eq_rhs(rhs_wY1Y2);

        cout << "--- [Finished] Making diff eq rhs ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Adjoint
        // ***************
        
        cout << "--- Making adjoint ---" << endl;
        
        // All biases
        std::map<int,std::map<Sptr,ITptr>> all_biases;
        all_biases[0][species_X] = hX;
        all_biases[0][species_Y] = hY;
        all_biases[1][species_X1] = bX1;
        all_biases[1][species_Y1] = bY1;
        all_biases[2][species_X2] = bX2;
        all_biases[2][species_Y2] = bY2;
        // All weights
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> all_weights;
        all_weights[0][species_X][1][species_X1] = wXX1;
        all_weights[1][species_X1][0][species_X] = wXX1;
        all_weights[0][species_Y][1][species_Y1] = wYY1;
        all_weights[1][species_Y1][0][species_Y] = wYY1;
        all_weights[1][species_X1][2][species_X2] = wX1X2;
        all_weights[2][species_X2][1][species_X1] = wX1X2;
        all_weights[1][species_Y1][2][species_Y2] = wY1Y2;
        all_weights[2][species_Y2][1][species_Y1] = wY1Y2;
        // Centers
        std::map<int,std::map<Sptr,CTptr>> all_center_trajs;
        all_center_trajs[0][species_X] = center_X;
        all_center_trajs[0][species_Y] = center_Y;
        all_center_trajs[1][species_X1] = center_X1;
        all_center_trajs[1][species_Y1] = center_Y1;
        all_center_trajs[2][species_X2] = center_X2;
        all_center_trajs[2][species_Y2] = center_Y2;
        // Conn mults
        std::map<int, std::map<int,int>> conn_mults;
        conn_mults[0][1] = 4;
        conn_mults[1][0] = 4;
        conn_mults[1][2] = 9;
        conn_mults[2][1] = 9;

        // bias
        auto deriv_term_hX = make_shared<AdjointParamsCenteredHomDerivTerm>(hX, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_hX = make_shared<AdjointParamsCenteredHomBias>("adjoint hX",hX,deriv_term_hX);
        hX->set_adjoint(adjoint_hX);

        auto deriv_term_hY = make_shared<AdjointParamsCenteredHomDerivTerm>(hY, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_hY = make_shared<AdjointParamsCenteredHomBias>("adjoint hY",hY,deriv_term_hY);
        hY->set_adjoint(adjoint_hY);

        auto deriv_term_bX1 = make_shared<AdjointParamsCenteredHomDerivTerm>(bX1, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_bX1 = make_shared<AdjointParamsCenteredHomBias>("adjoint bX1",bX1,deriv_term_bX1);
        bX1->set_adjoint(adjoint_bX1);

        auto deriv_term_bY1 = make_shared<AdjointParamsCenteredHomDerivTerm>(bY1, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_bY1 = make_shared<AdjointParamsCenteredHomBias>("adjoint bY1",bY1,deriv_term_bY1);
        bY1->set_adjoint(adjoint_bY1);

        auto deriv_term_bX2 = make_shared<AdjointParamsCenteredHomDerivTerm>(bX2, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_bX2 = make_shared<AdjointParamsCenteredHomBias>("adjoint bX2",bX2,deriv_term_bX2);
        bX2->set_adjoint(adjoint_bX2);
        
        auto deriv_term_bY2 = make_shared<AdjointParamsCenteredHomDerivTerm>(bY2, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_bY2 = make_shared<AdjointParamsCenteredHomBias>("adjoint bY2",bY2,deriv_term_bY2);
        bY2->set_adjoint(adjoint_bY2);

        // weights
        auto deriv_term_wXX1 = make_shared<AdjointParamsCenteredHomDerivTerm>(wXX1, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_wXX1 = make_shared<AdjointParamsCenteredHomWeight>("adjoint wXX1",wXX1,deriv_term_wXX1,deriv_term_hX,deriv_term_bX1,4,center_X,center_X1,adjoint_hX,adjoint_bX1,std::vector<CTptr>({center_X,center_Y}),std::vector<CTptr>({center_X1,center_Y1}));
        wXX1->set_adjoint(adjoint_wXX1);

        auto deriv_term_wYY1 = make_shared<AdjointParamsCenteredHomDerivTerm>(wYY1, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_wYY1 = make_shared<AdjointParamsCenteredHomWeight>("adjoint wYY1",wYY1,deriv_term_wYY1,deriv_term_hY,deriv_term_bY1,4,center_Y,center_Y1,adjoint_hY,adjoint_bY1,std::vector<CTptr>({center_X,center_Y}),std::vector<CTptr>({center_X1,center_Y1}));
        wYY1->set_adjoint(adjoint_wYY1);
        
        auto deriv_term_wX1X2 = make_shared<AdjointParamsCenteredHomDerivTerm>(wX1X2, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_wX1X2 = make_shared<AdjointParamsCenteredHomWeight>("adjoint wX1X2",wX1X2,deriv_term_wX1X2,deriv_term_bX1,deriv_term_bX2,4,center_X1,center_X2,adjoint_bX1,adjoint_bX2,std::vector<CTptr>({center_X1,center_Y1}),std::vector<CTptr>({center_X2,center_Y2}));
        wX1X2->set_adjoint(adjoint_wX1X2);
        
        auto deriv_term_wY1Y2 = make_shared<AdjointParamsCenteredHomDerivTerm>(wY1Y2, all_biases, all_weights, all_center_trajs, conn_mults);
        auto adjoint_wY1Y2 = make_shared<AdjointParamsCenteredHomWeight>("adjoint wY1Y2",wY1Y2,deriv_term_wY1Y2,deriv_term_bY1,deriv_term_bY2,4,center_Y1,center_Y2,adjoint_bY1,adjoint_bY2,std::vector<CTptr>({center_X1,center_Y1}),std::vector<CTptr>({center_X2,center_Y2}));
        wY1Y2->set_adjoint(adjoint_wY1Y2);
        
        cout << "--- [Finished] Making adjoint ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Lattice
        // ***************
        
        cout << "--- Making lattice ---" << endl;
        
        int no_dims = 2;
        int side_length = 40;
        latt = make_shared<LatticeTrajCenteredHom>(no_dims, side_length, std::vector<Sptr>({species_X,species_Y}), std::vector<CTptr>({center_X,center_Y}));

        cout << " > begin visible" << endl;
        
        // Add visible bias
        latt->set_bias_of_layer(0, species_X, hX);
        latt->set_bias_of_layer(0, species_Y, hY);

        cout << " > end visible" << endl;
        
        cout << " > begin hidden" << endl;
    
        // Add layer
        int conn_mult = 4;
        latt->add_layer(1, side_length, std::vector<Sptr>({species_X1,species_Y1}), std::vector<CTptr>({center_X1,center_Y1}), conn_mult);
        
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

                        latt->add_conn(0, i, j, 1, i3, j3);
                    };
                };
            };
        };
        
        // Biases
        latt->set_bias_of_layer(1, species_X1, bX1);
        latt->set_bias_of_layer(1, species_Y1, bY1);

        // Ixns
        latt->set_ixn_between_layers(0, species_X, 1, species_X1, wXX1);
        latt->set_ixn_between_layers(0, species_Y, 1, species_Y1, wYY1);
        
        // Layer 2
        latt->add_layer(2, side_length, std::vector<Sptr>({species_X2,species_Y2}), std::vector<CTptr>({center_X2,center_Y2}), conn_mult);
        
        // Connectivity
        for (auto i=1; i<=side_length; i++) {
            for (auto j=1; j<=side_length; j++) {
                
                // Displacements
                for (auto i2=0; i2<=2; i2++) {
                    for (auto j2=0; j2<=2; j2++) {
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
                        
                        latt->add_conn(1, i, j, 2, i3, j3);
                    };
                };
            };
        };
        
        // Biases
        latt->set_bias_of_layer(2, species_X2, bX2);
        latt->set_bias_of_layer(2, species_Y2, bY2);
        
        // Ixns
        latt->set_ixn_between_layers(1, species_X1, 2, species_X2, wX1X2);
        latt->set_ixn_between_layers(1, species_Y1, 2, species_Y2, wY1Y2);

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
