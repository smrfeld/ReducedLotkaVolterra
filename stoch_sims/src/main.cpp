#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <sstream>

// Get the header
#include <lattgillespie>

using namespace lattg;

// Function to get the mean of a vector
double get_mean(std::vector<int> v) {
	double ctr = 0.0;
	for (auto i: v) { ctr += i; };
	return ctr / v.size();
};

// Function to get the std dev of a vector
double get_std(std::vector<int> v) {
	double mean = get_mean(v);
	double ctr = 0.0;
	for (auto i: v) { ctr += pow(i - mean,2); };
	return sqrt(ctr / (v.size()-1));
};

/********************
Zero pad a string
********************/

std::string pad_str(int i, int n_zeros) {
	std::stringstream fname;
	fname << std::setfill('0') << std::setw(n_zeros) << i;
	return fname.str();
};

int main() {

	/********************
	Params
	********************/

	// Seed random no
	srand (time(NULL));
	// srand (2);

	// box length 30
	// kr1 0.012 -> 0.096
	// kr2 0.022 -> 0.176
	// prob3 0.1 -> 0.8

	// box length 40
	// kr1 0.025
	// kr2 0.06
	// prob3 0.4

	// Rates
	double kr1 = 0.025;
	double kr2 = 0.06;
	double prob3 = 0.4;

	// Box length
	int box_length = 40;

	// Timestep
	double dt = 1.0;

	// Number of steps to run
	// int n_steps = 100;
	int n_steps = 500;

	// Run
	bool verbose = true;
	bool write_counts = true;
	bool write_nns = true;
	bool write_latt = true;
	int write_step = 1;
	bool periodic_bc = true;

	// 200 particles ~ log (0.2/(1-3*0.2)) = - log(2)
	std::map<std::string,int> pops;
	pops["X"] = 100;
	pops["Y"] = 100;

	// Number of samples to run
	int i_start = 1;
	int i_end = 100; // inclusive
	std::vector<int> idxs({6,22,70});

	for (auto write_version_no: idxs) {
	// for (int write_version_no=i_start; write_version_no<=i_end; write_version_no++) {

		std::cout << "------ " << write_version_no << " -------" << std::endl;

		/****************************************
		Make a simulation!
		****************************************/

		// Make a simulation
		Simulation sim(dt,box_length,2);

		/********************
		Add species
		********************/

		sim.add_species("X");
		sim.add_species("Y");

		/********************
		Add unimol rxns
		********************/

		// P1
		sim.add_uni_rxn("rxn1", kr1, "", "X");
		// P3
		sim.add_uni_rxn("rxn2", kr2, "Y");
		// P5
		sim.add_bi_rxn("rxn3", prob3, "X", "Y", "Y","Y");

		/********************
		Populate the lattice
		********************/

		sim.populate_lattice(pops);

		/********************
		Run
		********************/

		sim.run(n_steps,verbose,write_counts,write_nns,write_latt,write_step,write_version_no,"../data",periodic_bc);
	};

	return 0;
}
