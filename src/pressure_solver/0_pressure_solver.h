#pragma once

#include "discretization/1_discretization.h"
#include <cmath>
#include <memory>

/**
 * Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 */

class PressureSolver
{
public:
    /**
    *  constructor
    * @param discretization: pointer to distrectization implementation
    * @param epsilon: error tolerance
    * @param maximumNumberOfIterations: maximal number of iterations before ending the iteration
    */
    PressureSolver(std::shared_ptr<Discretization> discretization,
                   double epsilon,
                   int maximumNumberOfIterations);
    
    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */

    virtual void solve() = 0;

    double residualNorm() const;

    int iterations() const;

protected:
    /**
     * set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
     */
    void setBoundaryValues();
    void setBoundaryValuesBottom();
    void setBoundaryValuesTop();
    void setBoundaryValuesLeft();
    void setBoundaryValuesRight();
    virtual void computeResidualNorm();
    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int maximumNumberOfIterations_;
    double residual_norm2_;
    int iterations_;
};
