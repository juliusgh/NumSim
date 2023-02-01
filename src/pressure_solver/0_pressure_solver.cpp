#include "pressure_solver/0_pressure_solver.h"

/**
 *  Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 * @param discretization: pointer to distrectization implementation
 * @param epsilon: error tolerance
 * @param maximumNumberOfIterations: maximal number of iterations before ending the iteration
 */

PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization,
                               double epsilon,
                               int maximumNumberOfIterations) :
        discretization_(discretization),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations) {

}

/**
 * Set pressure solver boundary values
 */
void PressureSolver::setBoundaryValues() {
    discretization_->applyBoundaryPressure();
}

void PressureSolver::setBoundaryValuesBottom() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        // copy values to bottom boundary
        discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pInteriorJBegin());
    }
}

void PressureSolver::setBoundaryValuesTop() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        // copy values to top boundary
        discretization_->p(i, discretization_->pJEnd() - 1) = discretization_->p(i,
                                                                                 discretization_->pInteriorJEnd() - 1);
    }
}

void PressureSolver::setBoundaryValuesLeft() {
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
        // copy values to left boundary
        discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pInteriorIBegin(), j);
    }
}

void PressureSolver::setBoundaryValuesRight() {
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
        // copy values to right boundary
        discretization_->p(discretization_->pIEnd() - 1, j) = discretization_->p(discretization_->pInteriorIEnd() - 1,
                                                                                 j);
    }
}

/**
 * Compute squared Euclidean norm to measure pressure solver convergence performance
 * ||X||_2^2 = x^2 + y^2 for X = (x,y)^T
 */
void PressureSolver::computeResidualNorm() {
    double residual_norm2 = 0.0;
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const int N = discretization_->nCells()[0] * discretization_->nCells()[1];
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            if (discretization_->marker(i, j) != FLUID) {
                continue;
            }
            double pxx =
                    (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2;
            double pyy =
                    (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2;
            residual_norm2 += pow(pxx + pyy - discretization_->rhs(i, j), 2);
        }
    }
    residual_norm2_ = residual_norm2 / N;
};

double PressureSolver::residualNorm() const {
    return residual_norm2_;
}

int PressureSolver::iterations() const {
    return iterations_;
}
