#pragma once

#include "pressure_solver/0_pressure_solver.h"

/**
 * Parallel red black solver for solving a linear system of equations.
 */

class RedBlack : public PressureSolver {
public:
    /**
    * constructor
    * @param discretization
    * @param epsilon
    * @param maximumNumberOfIterations
    * @param omega
    */
    RedBlack(std::shared_ptr <Discretization> discretization,
        double epsilon,
        int maximumNumberOfIterations,
        double omega);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve();

private:
    double omega_;
};
