#pragma once
#include "computation/0_computation.h"

class ComputationParallel : public Computation {
public:
    /**
     * Initialize the computation object
     *
     * Parse the settings from the parameter file that is given as the command line argument
     * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
     */
    void initialize(string filename);

    /**
     * Run the whole simulation until tend
     */
    void runSimulation();
};
