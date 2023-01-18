#pragma once

#include "pressure_solver/1_pressure_solver_parallel.h"

/**
 * Implementation of a parallelisized version of the Conjugated Gradients (CG) solver.
 */
class ConjugateGradient : public PressureSolverParallel {
public:
    /**
    * constructor
    * @param discretization pointer to the implementation of the discretization
    * @param epsilon error tolerance below which we consider the solver to be converged
    * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
    * @param partitioning information about subdomain
    */
    ConjugateGradient(std::shared_ptr <Discretization> discretization,
        double epsilon,
        int maximumNumberOfIterations,
        std::shared_ptr<Partitioning> partitioning);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve() override;

protected:
    std::shared_ptr<Array2D> q_;
    /**
    *  Implementation of communication of pressure values between neighbouring subdomains
    */
    void qGhostLayer();

    /**
    * set boundary values at the bottom of the subdomain for the search direction
    */
    void setQBoundaryValuesBottom();

    /**
    * set boundary values at the top of the subdomain for the search direction
    */    
    void setQBoundaryValuesTop();

    /**
    * set boundary values at the left of the subdomain for the search direction
    */
    void setQBoundaryValuesLeft();

    /**
    * set boundary values at the right of the subdomain for the search direction
    */
    void setQBoundaryValuesRight();
};
