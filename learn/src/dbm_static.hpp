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
    Sptr species_X, species_Y, species_X1, species_Y1, species_X2, species_Y2;
    
    // Ixn funcs
    Iptr hX, hY, bX1, bY1, bX2, bY2, wXX1, wYY1, wX1X2, wY1Y2;
    vector<Iptr> ixns;
    
    // Lattice
    shared_ptr<Lattice> latt;
    
    /********************
     Constructor
     ********************/
    
    DBM(int conn_mult_01, int conn_mult_12) {
        
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
        // MARK: - Ixn funcs
        // ***************

        cout << "--- Making ixn func ---" << endl;
        
        hX = make_shared<IxnParam>("hX", IxnParamType::H, 0.0, 0.0);
        hY = make_shared<IxnParam>("hY", IxnParamType::H, 0.0, 0.0);
        bX1 = make_shared<IxnParam>("bX1", IxnParamType::B, 0.0, 0.0);
        bY1 = make_shared<IxnParam>("bY1", IxnParamType::B, 0.0, 0.0);
        bX2 = make_shared<IxnParam>("bX2", IxnParamType::B, 0.0, 0.0);
        bY2 = make_shared<IxnParam>("bY2", IxnParamType::B, 0.0, 0.0);
        wXX1 = make_shared<IxnParam>("wXX1", IxnParamType::W, 0.0, 0.0);
        wYY1 = make_shared<IxnParam>("wYY1", IxnParamType::W, 0.0, 0.0);
        wX1X2 = make_shared<IxnParam>("wX1X2", IxnParamType::X, 0.0, 0.0);
        wY1Y2 = make_shared<IxnParam>("wY1Y2", IxnParamType::X, 0.0, 0.0);
        ixns = std::vector<Iptr>({hX,hY,bX1,bY1,bX2,bY2,wXX1,wYY1,wX1X2,wY1Y2});
        
        cout << "--- [Finished] Making ixn func ---" << endl;
        cout << endl;

        // ***************
        // MARK: - Lattice
        // ***************
        
        cout << "--- Making lattice ---" << endl;
        
        int no_dims = 2;
        int side_length = 40;
        latt = make_shared<Lattice>(no_dims, side_length, std::vector<Sptr>({species_X,species_Y}));

        cout << " > begin visible" << endl;
        
        // Add visible bias
        latt->set_bias_of_layer(0, species_X, hX);
        latt->set_bias_of_layer(0, species_Y, hY);

        cout << " > end visible" << endl;
        
        cout << " > begin hidden" << endl;
    
        // Add layer
        latt->add_layer(1, side_length, std::vector<Sptr>({species_X1,species_Y1}));
        
        // Connectivity
        for (auto i=1; i<=side_length; i++) {
            for (auto j=1; j<=side_length; j++) {
                
                // Displacements
                for (auto i2=0; i2<=sqrt(conn_mult_01)-1; i2++) {
                    for (auto j2=0; j2<=sqrt(conn_mult_01)-1; j2++) {
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
        latt->add_layer(2, side_length, std::vector<Sptr>({species_X2,species_Y2}));
        
        // Connectivity
        for (auto i=1; i<=side_length; i++) {
            for (auto j=1; j<=side_length; j++) {
                
                // Displacements
                for (auto i2=0; i2<=sqrt(conn_mult_12)-1; i2++) {
                    for (auto j2=0; j2<=sqrt(conn_mult_12)-1; j2++) {
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
