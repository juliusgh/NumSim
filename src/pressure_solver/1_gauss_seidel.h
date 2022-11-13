#pragma once

#include "discretization/1_discretization.h"
#include "pressure_solver/0_pressure_solver.h"

/** Standard Gauss-Seidel solver.
 */

class GaussSeidel : public PressureSolver {
public:
    //! constructor
    GaussSeidel(std::shared_ptr<Discretization> discretization,
                double epsilon,
                int maximumNumberOfIterations);

    //! solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
    void solve();
};
