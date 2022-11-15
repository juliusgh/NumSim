#include "settings.h"
#include "computation/computation.h"
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

    auto computation = Computation();
    computation.initialize(filename);
    computation.runSimulation();

    return EXIT_SUCCESS;
}

