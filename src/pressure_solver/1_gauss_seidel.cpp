#include "pressure_solver/1_gauss_seidel.h"
#include <cmath>
#include <iostream>

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
    auto dx = discretization_->dx();
    auto dy = discretization_->dy();
    double dx2 = pow(dx, 2);
    double dy2 = pow(dy, 2);
    double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    int iteration = 0;
    do {
        iteration++;

        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / pow(dx, 2);
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / pow(dy, 2);
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    } while (iteration < maximumNumberOfIterations_ && residualNorm() > pow(epsilon_, 2));
    iterations_ = iteration;
};