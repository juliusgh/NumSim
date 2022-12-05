#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include "mpi.h"
#include "pressure_solver/0_pressure_solver.h"
#include "partitioning/partitioning.h"
#include "1_pressure_solver_parallel.h"

/**
 * Parallel red black solver for solving a linear system of equations.
 */

class RedBlack : public PressureSolverParallel {
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
        std::shared_ptr<Partitioning> partitioning);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve();

protected:
    std::shared_ptr<Partitioning> partitioning_;
};
