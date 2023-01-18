#pragma once

#include "pressure_solver/2_red_black.h"
#include "1_pressure_solver_parallel.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the SOR solver.
 */

class RedBlackSOR : public PressureSolverParallel {
public:
    /**
    * constructor
    * @param discretization pointer to the implementation of the discretization
    * @param epsilon error tolerance below which we consider the solver to be converged
    * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
    * @param partitioning information about subdomain
    */
    RedBlackSOR(std::shared_ptr <Discretization> discretization,
        double epsilon,
        int maximumNumberOfIterations,
        double omega,
        std::shared_ptr<Partitioning> partitioning);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve() override;
private:
    double omega_;
};
