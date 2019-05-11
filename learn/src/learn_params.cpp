#include <iostream>
#include <dblz>
#include <fstream>
#include <vector>
#include <map>

#include <sstream>

#include "dbm_params.hpp"

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
    // MARK: - Init conds ixns
    // ***************
    
    map<string,double> init_conds_ixns;
    
    init_conds_ixns["hX"] = -2.63; // -2.63
    init_conds_ixns["hY"] = -2.63;
    init_conds_ixns["bX1"] = -2.63;
    init_conds_ixns["bY1"] = -2.63;
    init_conds_ixns["wXX1"] = 0.0;
    init_conds_ixns["wYY1"] = 0.0;

    // ***************
    // MARK: - Learning rates
    // ***************
    
    map<string,double> lrs;

    double lr_bias = 2.50; // 25
    double lr_weight = 2.50; // 25
    lrs["hX"] = lr_bias;
    lrs["hY"] = lr_bias;
    lrs["bX1"] = lr_bias;
    lrs["bY1"] = lr_bias;
    lrs["wXX1"] = lr_weight;
    lrs["wYY1"] = lr_weight;

    // ***************
    // MARK: - DBM
    // ***************
    
    double spacing = 0.5;
    auto dbm = DBM(spacing, init_conds_ixns, lrs);

    // ***************
    // MARK: - Options
    // ***************
    
    OptionsSolveDynamic options;
    OptionsWakeSleep_RBM_CD options_wake_sleep;
    
    options.verbose_timing = true;
    options_wake_sleep.verbose_timing = false;
    
    options.verbose = false;
        
    // Solver
    options.solver = Solver::SGD;
    
    // ***************
    // MARK: - Setup optimizer
    // ***************
    
    int no_timesteps_ixn_params = 10;
    
    int timepoint_start_wake_sleep = 0;
    int no_timesteps_wake_sleep = 10;
    
    int timepoint_start_adjoint = 0;
    int no_timesteps_adjoint = 10;

    // No markov chains
    int no_markov_chains = 5; // 10
    
    // Opt solver class
    OptProblemDynamic opt;
    
    // Params
    int no_opt_steps = 9800; // 5000
    int no_cd_steps = 10;
    double dt = 0.01;
    int no_steps_move = 20; // 20

    std::string dir = "../data/learn_params/";

    // ***************
    // MARK: - Set up the initial optimizer with the timepoints
    // ***************
    
    for (auto ixn: dbm.ixns) {
        ixn->set_no_timesteps(no_timesteps_ixn_params); // NOTE: MUST do this before dbm.latt->set_no_timesteps!
    };
    
    dbm.latt->set_no_markov_chains(MCType::AWAKE, no_markov_chains);
    dbm.latt->set_no_markov_chains(MCType::ASLEEP, no_markov_chains);
    dbm.latt->set_no_timesteps(timepoint_start_wake_sleep, no_timesteps_wake_sleep);

    // ***************
    // MARK: - Filenames
    // ***************
    
    FNameTrajColl fnames;
    int timepoint_max = 500;
    for (auto i_batch=1; i_batch<=100; i_batch++) {
        // Add
        bool binary_read = true;
        FNameTraj fname_traj;
        for (auto timepoint=0; timepoint<=timepoint_max; timepoint++) {
            fname_traj.push_back(FName("../../stoch_sims/data/lattice_v" + pad_str(i_batch,3) + "/lattice/"+pad_str(timepoint,4)+".txt",binary_read));
        };
        
        fnames.add_fname_traj(fname_traj);
    };
    
    // ***************
    // MARK: - Solve
    // ***************
    
    for (auto opt_step=1; opt_step<=no_opt_steps; opt_step++) {

        cout << "------------------" << endl;
        cout << "Opt step: " << opt_step << " / " << no_opt_steps << endl;
        cout << "------------------" << endl;
        
        // Move forward?
        if (opt_step > no_steps_move && opt_step % no_steps_move == 1) {
            
            // Ixn params
            no_timesteps_ixn_params += 1;
            for (auto ixn: dbm.ixns) {
                ixn->set_no_timesteps(no_timesteps_ixn_params);
            };
            
            // Move adjoint forward
            timepoint_start_adjoint += 1;
            
            // Latt
            timepoint_start_wake_sleep = timepoint_start_adjoint;
            dbm.latt->set_no_timesteps(timepoint_start_wake_sleep, no_timesteps_wake_sleep);
        };
        
        // Solve
        opt.solve_one_step_rbm_cd_params(dbm.latt, opt_step, 0, no_timesteps_ixn_params, timepoint_start_wake_sleep, no_timesteps_wake_sleep, timepoint_start_adjoint, no_timesteps_adjoint, dt, no_cd_steps, fnames, options, options_wake_sleep);
        
        // Write
        if (opt_step == 1 || opt_step % 100 == 0 || opt_step == no_opt_steps) {
            int i_write = opt_step;
            if (opt_step == 1) {
                i_write = 0;
            };
            
            // Moments
            for (auto ixn: dbm.ixns) {
                ixn->write_moment_traj_to_file(0, no_timesteps_ixn_params, dir + "moments/" + ixn->get_name() + "_" + pad_str(i_write,5) + ".txt");
            };
            
            // Ixns
            for (auto ixn: dbm.ixns) {
                ixn->write_val_traj_to_file(0, no_timesteps_ixn_params, dir + "ixn_params/" + ixn->get_name() + "_" + pad_str(i_write,5) + ".txt");
            };
        };

	// Write the diff eqs (big!)
	if (opt_step == no_opt_steps) {
	    for (auto ixn: dbm.ixns) {
                ixn->get_diff_eq_rhs()->write_to_file(dir + "diff_eq_rhs/" + ixn->get_name() + "_" + pad_str(no_opt_steps,5) + ".txt");
            };
	};
    };
    
	return 0;
};

