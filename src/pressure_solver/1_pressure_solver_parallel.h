#pragma once

#include "discretization/1_discretization.h"
#include "pressure_solver/0_pressure_solver.h"
#include <cmath>
#include <memory>
#include <cmath>
#include <iostream>
#include <vector>
#include "mpi.h"
#include "partitioning/partitioning.h"

/**
 * Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 */

class PressureSolverParallel : public PressureSolver
{
public:
    /**
    *  constructor
    * @param discretization: pointer to distrectization implementation
    * @param epsilon: error tolerance
    * @param maximumNumberOfIterations: maximal number of iterations before ending the iteration
    */
    PressureSolverParallel(std::shared_ptr<Discretization> discretization,
                   double epsilon,
                   int maximumNumberOfIterations,
                   std::shared_ptr<Partitioning> partitioning);
    
    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */

protected:
    virtual void computeResidualNorm();
    void pGhostLayer();
    std::shared_ptr<Partitioning> partitioning_;
};
