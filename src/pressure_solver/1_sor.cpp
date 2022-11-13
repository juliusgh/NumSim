#include "pressure_solver/1_sor.h"
#include <cmath>
#include <iostream>

SOR::SOR(std::shared_ptr <Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         double omega) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations),
        omega_(omega)
{

}

void SOR::solve() {
    auto dx = discretization_->dx();
    auto dy = discretization_->dy();
    double dx2 = pow(dx, 2);
    double dy2 = pow(dy, 2);
    double k1 = 1 - omega_;
    double k2 = omega_ * (dx2 * dy2) / (2.0 * (dx2 + dy2));
    int iteration = 0;
    do {
        iteration++;

        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k1 * discretization_->p(i, j) + k2 * (px + py - discretization_->rhs(i, j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    } while (iteration < maximumNumberOfIterations_ && residualNorm() > pow(epsilon_, 2));
    iterations_ = iteration;
};