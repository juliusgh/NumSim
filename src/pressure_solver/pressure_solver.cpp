#include "pressure_solver/pressure_solver.h"

PressureSolver::PressureSolver(std::shared_ptr <Discretization> discretization,
                               double epsilon,
                               int maximumNumberOfIterations) :
        discretization_(discretization),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations) {

}

void PressureSolver::setBoundaryValues() {
    for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {

        // copy values to left boundary
        discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pInteriorIBegin(), j);

        // copy values to right boundary
        discretization_->p(discretization_->pIEnd() - 1, j) = discretization_->p(discretization_->pInteriorIEnd() - 1, j);
    }
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {

        // copy values to bottom boundary
        discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pInteriorJBegin());

        // copy values to top boundary
        discretization_->p(i, discretization_->pJEnd() - 1) = discretization_->p(i, discretization_->pInteriorJEnd() - 1);
    }
};

double PressureSolver::residualNorm() const {
    return residual_norm2_;
}

int PressureSolver::iterations() const {
    return iterations_;
}