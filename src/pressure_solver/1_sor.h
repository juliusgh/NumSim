#pragma once

#include "pressure_solver/0_pressure_solver.h"

/**
 * Successive over-relaxation solver for solving a linear system of equations.
 */

class SOR : public PressureSolver {
public:
    //! constructor
    SOR(std::shared_ptr <Discretization> discretization,
        double epsilon,
        int maximumNumberOfIterations,
        double omega);

    //! solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
    void solve();

private:
    double omega_;
};
