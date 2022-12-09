#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include "mpi.h"
#include "pressure_solver/0_pressure_solver.h"
#include "partitioning/partitioning.h"
#include "1_pressure_solver_parallel.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the Gauss-Seidel solver.
 */

class RedBlack : public PressureSolverParallel {
public:
    /**
    * constructor
    * @param discretization pointer to the implementation of the discretization
    * @param epsilon error tolerance below which we consider the solver to be converged
    * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
    * @param partitioning information about subdomain
    */
    RedBlack(std::shared_ptr <Discretization> discretization,
        double epsilon,
        int maximumNumberOfIterations,
        std::shared_ptr<Partitioning> partitioning);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve();

protected:
    std::shared_ptr<Partitioning> partitioning_;
};
