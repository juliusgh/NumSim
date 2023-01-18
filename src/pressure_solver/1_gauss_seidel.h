#pragma once

#include <cmath>
#include <iostream>
#include "discretization/1_discretization.h"
#include "pressure_solver/0_pressure_solver.h"

/**
 * Standard Gauss-Seidel solver for linear systems of equations..
 */

class GaussSeidel : public PressureSolver {
public:
    /**
     * constructor
     * @param discretization: pointer to distrectization implementation
     * @param epsilon: error tolerance
     * @param maximumNumberOfIterations: maximal number of iterations before ending the iteration
     */
    GaussSeidel(std::shared_ptr<Discretization> discretization,
                double epsilon,
                int maximumNumberOfIterations);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve();
};
