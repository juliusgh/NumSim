#include "pressure_solver/2_red_black_sor.h"

/**
 * Implementation of the red-black solver for SOR, a parallelisized version of the SOR solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations if not converged when this number is reached, solver is terminated nonetheless 
 */

RedBlackSOR::RedBlackSOR(std::shared_ptr<Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         double omega,
         std::shared_ptr<Partitioning> partitioning) :
RedBlack(discretization, epsilon, maximumNumberOfIterations, partitioning),
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
    do {
        iteration++;
        //getchar();
        /*if (partitioning_->ownRankNo() == 0) {
            discretization_->p().print();
        }
        if (getchar() == '1' && partitioning_->ownRankNo() == 1) {
            discretization_->p().print();
        }*/

        // black half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                //std::cout << "BLACK: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k1 * discretization_->p(i, j) + k2 * (px + py - discretization_->rhs(i, j));
            }
        }
        //std::cout << "RANK " << partitioning_->ownRankNo() << " : start pGhostLayer 1" << std::endl;

        // communication
        pGhostLayer();

        // red half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pInteriorJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                //std::cout << "RED: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k1 * discretization_->p(i, j) + k2 * (px + py - discretization_->rhs(i, j));
            }
        }

        //std::cout << "RANK " << partitioning_->ownRankNo() << " : start pGhostLayer 2" << std::endl;
        pGhostLayer();
        computeResidualNorm();
        //partitioning_->log("residual norm final:");
        //std::cout << "RANK " << partitioning_->ownRankNo() << " : " << residualNorm() << std::endl;
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
}
