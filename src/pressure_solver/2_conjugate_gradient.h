#pragma once

#include "pressure_solver/1_pressure_solver_parallel.h"

/**
 * Parallel red black solver for solving a linear system of equations.
 */

class ConjugateGradient : public PressureSolverParallel {
public:
    /**
    * constructor
    * @param discretization
    * @param epsilon
    * @param maximumNumberOfIterations
    * @param omega
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
    virtual void computeResidualNorm();
    std::shared_ptr<Array2D> residual_;
    std::shared_ptr<Array2D> q_;
    std::shared_ptr<Array2D> Aq_;
    void qGhostLayer();
    void setQBoundaryValuesBottom();
    void setQBoundaryValuesTop();
    void setQBoundaryValuesLeft();
    void setQBoundaryValuesRight();
};
