#pragma once

#include "pressure_solver/0_pressure_solver.h"

/**
 * Implementation of a parallelisized version of the Conjugated Gradients (CG) solver.
 */
class ConjugateGradient : public PressureSolver {
public:
    /**
    * constructor
    * @param discretization pointer to the implementation of the discretization
    * @param epsilon error tolerance below which we consider the solver to be converged
    * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
    * @param partitioning information about subdomain
    */
    ConjugateGradient(std::shared_ptr<Discretization> discretization,
                              double epsilon,
                              int maximumNumberOfIterations
                              );

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve() override;

protected:

    /**
      * get reference to field variable p
     */
    const FieldVariable &s() const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double s(int i, int j) const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double &s(int i, int j);
    /**
  * get reference to field variable p
 */
    const FieldVariable &residual() const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double residual(int i, int j) const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double &residual(int i, int j);

    /**
  * get reference to field variable p
 */
    const FieldVariable &As() const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double As(int i, int j) const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double &As(int i, int j);

    /**
    *  Implementation of communication of pressure values between neighbouring subdomains
    */
    void applyBoundarySearchDirection();

private:
    FieldVariable s_;
    FieldVariable residual_;
    FieldVariable As_;

};
