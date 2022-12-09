#include "pressure_solver/2_red_black_sor.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the SOR solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
 * @param partitioning information about subdomain
 */

RedBlackSOR::RedBlackSOR(std::shared_ptr<Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         double omega,
         std::shared_ptr<Partitioning> partitioning) :
PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning),
omega_(omega)
{

}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void RedBlackSOR::solve() {
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double k1 = 1 - omega_;
    const double k2 = omega_ * (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = pow(epsilon_, 2);
    int iteration = 0;

    const int pInteriorIBegin = discretization_->pInteriorIBegin();
    const int pInteriorIEnd = discretization_->pInteriorIEnd();
    const int pInteriorJBegin = discretization_->pInteriorJBegin();
    const int pInteriorJEnd = discretization_->pInteriorJEnd();
    int offset;
    if (((partitioning_->nodeOffset()[0] % 2) + (partitioning_->nodeOffset()[0] % 2)) % 2 == 0)
    {
        offset = 0;
    } else {
        offset = 1;
    }
    do {
        iteration++;

        // black half step
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            int iStart = pInteriorIBegin + (j + offset) % 2;
            for (int i = iStart; i < pInteriorIEnd; i += 2) {
                //std::cout << "BLACK: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k1 * discretization_->p(i, j) + k2 * (px + py - discretization_->rhs(i, j));
            }
        }
        pGhostLayer();

        // red half step
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            int iStart = pInteriorIBegin + (j + 1 + offset) % 2;
            for (int i = iStart; i < pInteriorIEnd; i += 2) {
                //std::cout << "RED: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k1 * discretization_->p(i, j) + k2 * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayer();
        computeResidualNorm();
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
}
