#include "settings.h"
#ifndef NPARALLEL
    #include "computation/1_computation_parallel.h"
#else
    #include "computation/0_computation.h"
#endif
#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

/**
 * The entry point of the simulation.
 * 
 * Call with the parameter file as command line argument.
 * If no argument is given the user is asked to specify the parameter file.
 * The output files are written to ../out/*
 */
int main(int argc, char *argv[]) {
    // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
    string program = argv[0];
    string filename = "";
    if (argc == 1) {
        cout << "enter parameter file path: " << endl;
        cin >> filename;
    }
    else {
        // read in the first argument
        filename = argv[1];
    }
#ifndef NPARALLEL
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    cout << "ComputationParallel" << endl;
    auto computation = ComputationParallel();
#else
    auto computation = Computation();
#endif
    cout << "computation.initialize_" << endl;
    computation.initialize(filename);
    computation.runSimulation();
#ifndef NPARALLEL
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}

