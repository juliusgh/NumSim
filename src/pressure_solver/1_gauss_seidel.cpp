#include "pressure_solver/1_gauss_seidel.h"

/**
 * Standard Gauss-Seidel solver for linear systems of equations.
 * @param discretization
 * @param epsilon
 * @param maximumNumberOfIterations
 */

GaussSeidel::GaussSeidel(std::shared_ptr <Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{

}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void GaussSeidel::solve() {
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = pow(epsilon_, 2);
    int iteration = 0;
    do {
        iteration++;

        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
};