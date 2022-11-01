#pragma once

#include "discretization/1_discretization.h"
#include "pressure_solver/pressure_solver.h"

/** Successive over-relaxation solver.
 */

class SOR : PressureSolver
{
public:
    //! constructor
    using PressureSolver::PressureSolver;

    //! solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
    void solve();

private:
    double omega_;
};
