#include "pressure_solver/1_conjugate_gradient.h"

/**
 * Implementation of a parallelisized version of the Conjugated Gradients (CG) solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
 * @param partitioning information about subdomain
 */

ConjugateGradient::ConjugateGradient(std::shared_ptr<Discretization> discretization,
                                                     double epsilon,
                                                     int maximumNumberOfIterations) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations) {
    q_ = std::make_unique<Array2D>(discretization_->pSize());// Search direction qₖ
}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void ConjugateGradient::solve() {
#ifndef NDEBUG
    std::cout << "Solving pressure using Conjugate Gradient method" << std::endl;
#endif
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double eps2 = pow(epsilon_, 2);
    const int pIBegin = discretization_->pIBegin();
    const int pJBegin = discretization_->pJBegin();
    const int pIEnd = discretization_->pIEnd();
    const int pJEnd = discretization_->pJEnd();
    const int pIIntBegin = discretization_->pInteriorIBegin();
    const int pJIntBegin = discretization_->pInteriorJBegin();
    const int pIIntEnd = discretization_->pInteriorIEnd();
    const int pJIntEnd = discretization_->pInteriorJEnd();
    const int N = discretization_->pSize()[0] * discretization_->pSize()[1];


    Array2D residual_ = Array2D(discretization_->pSize());
    Array2D Aq_ = Array2D(discretization_->pSize());

    int iteration = 0;
    // Initialization Loop
    double alpha = 0.0;
    for (int i = pIIntBegin; i < pIIntEnd; i++) {
        for (int j = pJIntBegin; j < pJIntEnd; j++) {
            double D2pDx2 =
                    (discretization_->p(i - 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / dx2;
            double D2pDy2 =
                    (discretization_->p(i, j - 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / dy2;
            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            residual_(i - pIBegin, j - pJBegin) = discretization_->rhs(i, j) - (D2pDx2 + D2pDy2);

            // set search direction to preconditioned defect q = z
            (*q_)(i - pIBegin, j - pJBegin) = residual_(i - pIBegin, j - pJBegin);


            alpha += pow(residual_(i - pIBegin, j - pJBegin), 2);
        }
    }
#ifndef NDEBUG
    std::cout << "Initial residual: " << sqrt(alpha) << std::endl;
#endif

    do {

        qGhostLayer();
        double lambda = 0.0;
        // Calculate auxillary variable Aq
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                double D2qDx2 = ((*q_)(i - 1, j) - 2 * (*q_)(i, j) + (*q_)(i + 1, j)) / dx2;
                double D2qDy2 = ((*q_)(i, j - 1) - 2 * (*q_)(i, j) + (*q_)(i, j + 1)) / dy2;
                Aq_(i, j) = D2qDx2 + D2qDy2;

                // qₖᵀAqₖ
                lambda += (*q_)(i, j) * Aq_(i, j);
            }
        }

        // λ = αₖ / qₖᵀAqₖ
        lambda = alpha / lambda;
        iteration++;

        double alphaold = alpha;
        alpha = 0.0;
        // Update variables in the search direction
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {

                // pₖ₊₁ = pₖ + λ qₖ
                discretization_->p(i + pIBegin, j + pIBegin) += lambda * (*q_)(i, j);

                // rₖ₊₁ = rₖ - λ Aqₖ
                residual_(i, j) -= lambda * Aq_(i, j);

                // αₖ₊₁ = rₖ₊₁ᵀ rₖ₊₁
                alpha += pow(residual_(i, j), 2);
            }
        }

        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                (*q_)(i, j) = residual_(i, j) + beta * (*q_)(i, j);             // qₖ₊₁ = rₖ₊₁ + β qₖ
            }
        }
        residual_norm2_ = alpha / N;

    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);

    discretization_->applyBoundaryPressure();

    iterations_ = iteration;
}

/**
 *  Implementation of communication of search directions q between neighbouring subdomains
 */
void ConjugateGradient::qGhostLayer() {
    setQBoundaryValuesTop();
    setQBoundaryValuesBottom();
    setQBoundaryValuesRight();
    setQBoundaryValuesLeft();

}

/**
 * set boundary values at the bottom of the subdomain for the search direction
*/
void ConjugateGradient::setQBoundaryValuesBottom() {
    for (int i = 0; i < discretization_->pIEnd() - discretization_->pIBegin(); i++) {
        // copy values to bottom boundary
        (*q_)(i, 0) = (*q_)(i, 1);
    }
}

/**
 * set boundary values at the top of the subdomain for the search direction
*/
void ConjugateGradient::setQBoundaryValuesTop() {
    for (int i = 0; i < discretization_->pIEnd() - discretization_->pIBegin(); i++) {
        // copy values to top boundary
        (*q_)(i, discretization_->pJEnd() - 1 - discretization_->pJBegin()) = (*q_)(i,
                                                                                    discretization_->pInteriorJEnd() -
                                                                                    1 - discretization_->pJBegin());
    }
}

/**
 * set boundary values at the left of the subdomain for the search direction
*/
void ConjugateGradient::setQBoundaryValuesLeft() {
    for (int j = 0; j < discretization_->pJEnd() - discretization_->pJBegin(); j++) {
        // copy values to left boundary
        (*q_)(0, j) = (*q_)(1, j);
    }
}

/**
 * set boundary values at the right of the subdomain for the search direction
*/
void ConjugateGradient::setQBoundaryValuesRight() {
    for (int j = 0; j < discretization_->pJEnd() - discretization_->pJBegin(); j++) {
        // copy values to right boundary
        (*q_)(discretization_->pIEnd() - 1 - discretization_->pIBegin(), j) = (*q_)(
                discretization_->pInteriorIEnd() - 1 - discretization_->pIBegin(), j);
    }
}