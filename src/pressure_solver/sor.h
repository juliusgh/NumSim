#pragma once

#include "pressure_solver/pressure_solver.h"

/** Successive over-relaxation solver.
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
