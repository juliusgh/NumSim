#include "pressure_solver/gauss_seidel.h"
#include <cmath>
#include <iostream>

GaussSeidel::GaussSeidel(std::shared_ptr <Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{

}

void GaussSeidel::solve() {
    auto dx = discretization_->dx();
    auto dy = discretization_->dy();
    auto nCells = discretization_->nCells();
    int N = nCells[0] * nCells[1];
    int iteration = 0;
    double residual_norm2 = 0.0;
    do {
        iteration++;

        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++) {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / pow(dx, 2);
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / pow(dy, 2);
                discretization_->p(i, j) = pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
                        (px + py - discretization_->rhs(i, j));
            }
        }
        // stopping criterion
        residual_norm2 = 0.0;
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++) {
                double pxx = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / pow(dx, 2);
                double pyy = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / pow(dy, 2);
                residual_norm2 += pow(pxx + pyy - discretization_->rhs(i, j), 2);
            }
        }
    } while (iteration < maximumNumberOfIterations_ && residual_norm2 / N > pow(epsilon_, 2));
    setBoundaryValues();
};