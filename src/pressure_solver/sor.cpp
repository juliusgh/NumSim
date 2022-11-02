#include "pressure_solver/sor.h"
#include <cmath>

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
    auto p = discretization_->p();
    auto rhs = discretization_->rhs();
    auto nCells = discretization_->nCells();
    int N = nCells[0] * nCells[1];
    int iteration = 0;
    double residual_norm2 = 0.0;
    while (iteration < maximumNumberOfIterations_ && residual_norm2 / N <= pow(epsilon_, 2)) {
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
                p(i, j) = (1 - omega_) * p(i, j) + omega_ * pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
                          ((p(i - 1, j) + p(i + 1, j)) / pow(dx, 2) +
                           (p(i, j - 1) + p(i, j + 1)) / pow(dy, 2) - rhs(i, j));
            }
        }
        // stopping criterion
        double residual_norm2 = 0.0;
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
                residual_norm2 += pow((p(i + 1, j) - 2 * p(i, j) + p(i - 1, j)) / pow(dx, 2) +
                        (p(i, j + 2) - 2 * p(i, j) + p(i, j - 1)) / pow(dx, 2) - rhs(i, j), 2);
            }
        }
    }
};