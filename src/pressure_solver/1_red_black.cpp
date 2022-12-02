#include "pressure_solver/1_red_black.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the Gauss-Seidel solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations
 */

RedBlack::RedBlack(std::shared_ptr<Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         std::shared_ptr<Partitioning> partitioning) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations),
partitioning_(partitioning)
{

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

        // red half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j += 2) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                std::cout << "RED: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayer();

        // black half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j += 2) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pInteriorJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                std::cout << "BLACK: i = " << i << ", j = " << j << std::endl;
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



/**
 * Calculation of the global residual norm using the squared Euclidean norm
 */
void RedBlack::computeResidualNorm() {
    double residual_norm2 = 0.0;
    const double dx2 = pow(discretization_->dx(),2);
    const double dy2 = pow(discretization_->dy(),2);
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            double pxx = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2;
            double pyy = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2;
            residual_norm2 += pow(pxx + pyy - discretization_->rhs(i, j), 2);
        }
    }

    double residual_norm2_global = partitioning_->globalSum(residual_norm2);
    residual_norm2_ = residual_norm2_global / N;
}