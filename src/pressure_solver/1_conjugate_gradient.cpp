#include "pressure_solver/1_conjugate_gradient.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the Gauss-Seidel solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations
 */

ConjugateGradient::ConjugateGradient(std::shared_ptr<Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         std::shared_ptr<Partitioning> partitioning) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations),
partitioning_(partitioning),
residual_(discretization_->pSize)
{

}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void ConjugateGradient::solve() {
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = pow(epsilon_, 2);

    double alpha = 0;
    double beta = 0;
    Array2D q = Array2D(discretization_->pSize); // Search direction qₖ
    Array2D Aq = Array2D(discretization_->pSize); //precomputed matrix-vector product between the system matrix and the search direction Aqₖ

    int iteration = 0;

    // Calculate initial residuum r(i,j) = rhs(i,j) - (Ap)(i,j) = rhs(i,j) - (D2pDx2 + D2pDy2) and set inital search direction q(i,j) = r(i,j)
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            double D2pDx2 = (discretization_->p(i-1, j) - 2 * discretization_->p(i, j) + discretization_->p(i+1, j)) / dx2;
            double D2pDy2 = (discretization_->p(i, j-1) - 2 * discretization_->p(i, j) + discretization_->p(i, j+1)) / dy2;
            residual_(i, j) = discretization_->rhs(i,j) - (D2pDx2 + D2pDy2); 
            q(i, j) = residual_(i, j); 
        }
    }
    
    do {
        //Reset local quantities
        double local_residual_product = 0;      //local scalar product of the residual with itself rₖᵀ rₖ
        double local_newResidual_product = 0;   //local scalar product of the new residual with itself rₖ₊₁ᵀ rₖ₊₁
        double local_search_product = 0;        //local scalar product between the search direction with the precomputed matrix-vector product qₖᵀAqₖ
        iteration++;

        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                double D2qDx2 = (q(i-1, j) - 2 * q(i, j) + q(i+1, j)) / dx2;
                double D2qDy2 = (q(i, j-1) - 2 * q(i, j) + q(i, j+1)) / dy2;
                Aq(i, j) = (D2qDx2 + D2qDy2);

                local_search_product += q(i, j) * Aq(i, j);             // qₖᵀAqₖ
                local_residual_product += pow(residual_(i, j), 2);      // rₖᵀ rₖ
            }
        }
        double global_residual_product = partitioning_->globalSum(local_residual_product);
        double global_search_product = partitioning_->globalSum(local_search_product);

        alpha = global_residual_product / global_search_product;        // α = rₖᵀ rₖ / qₖᵀAqₖ

        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(i, j)= discretization_->p(i ,j) + alpha * q(i, j); // pₖ₊₁ = pₖ + α qₖ
                residual_(i, j) = residual_(i, j) - alpha * Aq(i ,j);   // rₖ₊₁ = rₖ - α Aqₖ
                local_newResidual_product += pow(residual(i, j), 2)     // rₖ₊₁ᵀ rₖ₊₁
            }
        }
        double global_newResidual_product = partitioning_->globalSum(local_newResidual_product);

        beta = global_newResidual_product / global_residual_product;

        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                q(i, j) = residual_(i, j) + beta * q(i, j);             // qₖ₊₁ = rₖ₊₁ + β qₖ
            }
        }
        computeResidualNorm();

        pGhostLayer();

    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
}

/**
 * Calculation of the global residual norm using the squared Euclidean norm
 */
void ConjugateGradient::computeResidualNorm() {
    double residual_norm2 = 0.0;
    const double dx2 = pow(discretization_->dx(),2);
    const double dy2 = pow(discretization_->dy(),2);
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            residual_norm2 += pow(residual_(i, j), 2);
        }
    }

    double residual_norm2_global = partitioning_->globalSum(residual_norm2);
    residual_norm2_ = residual_norm2_global / N;
}