#include "pressure_solver/pressure_solver.h"

PressureSolver::PressureSolver(std::shared_ptr <Discretization> discretization,
                               double epsilon,
                               int maximumNumberOfIterations) :
        discretization_(discretization),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations) {

}

void PressureSolver::setBoundaryValues() {
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++) {

        // copy values to left boundary
        discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pIBegin() + 1, j);

        // copy values to right boundary
        discretization_->p(discretization_->pIEnd() , j) = discretization_->p(discretization_->pIEnd() - 1, j);
    }
    for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++) {

        // copy values to bottom boundary
        discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pJBegin() + 1);

        // copy values to top boundary
        discretization_->p(i, discretization_->pJEnd()) = discretization_->p(i, discretization_->pJEnd() - 1);
    }
};

double PressureSolver::residualNorm() const {
    return residual_norm2_;
}

int PressureSolver::iterations() const {
    return iterations_;
}