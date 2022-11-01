#pragma once

#include "discretization/1_discretization.h"

/** Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 */

class PressureSolver
{
public:
    //! constructor
    PressureSolver(std::shared_ptr<Discretization> discretization,
                   double epsilon,
                   int maximumNumberOfIterations);
    
    //! solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
    virtual void solve() = 0;

protected:
    //! set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
    void setBoundaryValues();
    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int maximumNumberOfIterations_;
};
