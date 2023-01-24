#include "pressure_solver/2_red_black.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the Gauss-Seidel solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
 * @param partitioning information about subdomain
 */

RedBlack::RedBlack(std::shared_ptr<Discretization> discretization,
                   double epsilon,
                   int maximumNumberOfIterations,
                   std::shared_ptr<Partitioning> partitioning) :
        PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning) {

}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void RedBlack::solve() {
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = pow(epsilon_, 2);
    int iteration = 0;
    do {
        iteration++;

        // black half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayer();

        // red half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pInteriorJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayer();
        computeResidualNorm();
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
}
