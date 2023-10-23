#pragma once

#include "pressure_solver/1_conjugate_gradient.h"

/**
 * Implementation of a parallelisized version of the Conjugated Gradients (CG) solver.
 */
class ConjugateGradientPreconditioned : public ConjugateGradient {
public:
    /**
    * constructor
    * @param discretization pointer to the implementation of the discretization
    * @param epsilon error tolerance below which we consider the solver to be converged
    * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
    * @param partitioning information about subdomain
    */
    ConjugateGradientPreconditioned(std::shared_ptr<Discretization> discretization,
                              double epsilon,
                              int maximumNumberOfIterations
                              );

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve() override;

protected:
    /**
    * get reference to field variable preconditioned variable
    */
    const FieldVariable &z() const;

    /**
     * evaluate field variable preconditioned variable in an element (i,j)
    */
    double z(int i, int j) const;

    /**
     * evaluate field variable preconditioned variable in an element (i,j)
    */
    double &z(int i, int j);


private:
    FieldVariable z_;

};
