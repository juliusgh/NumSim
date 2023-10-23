#include "pressure_solver/2_conjugate_gradient_preconditioned.h"

/**
 * Implementation of a preconditioned version of the Conjugated Gradients (CG) solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
 * @param partitioning information about subdomain
 */

ConjugateGradientPreconditioned::ConjugateGradientPreconditioned(std::shared_ptr<Discretization> discretization,
                                                     double epsilon,
                                                     int maximumNumberOfIterations) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations),
        s_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth()),// Search direction sₖ
        residual_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth()),
        As_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth()),
        z_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth())
        {
}

void ConjugateGradientPreconditioned::solve() {
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double eps2 = pow(epsilon_, 2);
    const int pIIntBegin = discretization_->pInteriorIBegin();
    const int pJIntBegin = discretization_->pInteriorJBegin();
    const int pIIntEnd = discretization_->pInteriorIEnd();
    const int pJIntEnd = discretization_->pInteriorJEnd();
    const int N = discretization_->pSize()[0] * discretization_->pSize()[1];

    s_.setToZero();
    residual_.setToZero();
    As_.setToZero();
    z_.setToZero();

    int iteration = 0;

    // Initialization Loop
    double alpha = 0.0;
    double residual_norm = 0.0;
    for (int i = pIIntBegin; i < pIIntEnd; i++) {
        for (int j = pJIntBegin; j < pJIntEnd; j++) {
            if (!discretization_->isInnerFluid(i, j)) {
                continue;
            }

            double D2pDx2 =
                    (discretization_->p(i - 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / dx2;
            double D2pDy2 =
                    (discretization_->p(i, j - 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / dy2;

            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            residual(i, j) = discretization_->rhs(i, j) - (D2pDx2 + D2pDy2);

            // set search direction to preconditioned defect s = z
            //TODO: Preconditioning with s= P^{-1} residual
            s(i, j) = residual(i, j);

            alpha += residual(i, j) * s(i, j);

            residual_norm += pow(residual(i, j), 2);
        }
    }
    residual_norm2_ = residual_norm;


    do {
        applyBoundarySearchDirection();

        double lambda = 0.0;
        // Calculate auxiliary variable As
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (!discretization_->isInnerFluid(i, j)) {
                    continue;
                }
                double D2sDx2 = (s(i - 1, j) - 2 * s(i, j) + s(i + 1, j)) / dx2;
                double D2sDy2 = (s(i, j - 1) - 2 * s(i, j) + s(i, j + 1)) / dy2;
                As(i, j) = D2sDx2 + D2sDy2;

                // sₖᵀAsₖ
                lambda += s(i, j) * As(i, j);
            }
        }

        // λ = αₖ / sₖᵀAsₖ
        lambda = alpha / lambda;
        iteration++;

        // Update variables in the search direction
        double alphaold = alpha;
        residual_norm = 0.0;
        alpha = 0.0;
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (!discretization_->isInnerFluid(i, j)) {
                    continue;
                }
                // pₖ₊₁ = pₖ + λ sₖ
                discretization_->p(i, j) += lambda * s(i, j);

                // rₖ₊₁ = rₖ - λ Asₖ
                residual(i, j) -= lambda * As(i, j);

                residual_norm += pow(residual(i, j), 2);

                //TODO: Implement the preconditioning z = p^{-1} residual
                z(i, j) = residual(i, j);

                alpha += residual(i, j) * z (i, j);
            }
        }
        residual_norm2_ = residual_norm;

        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i, j) != FLUID) {
                    continue;
                }
                s(i, j) = z(i, j) + beta * s(i, j);             // sₖ₊₁ = zₖ₊₁ + β sₖ
            }
        }


    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);

    discretization_->applyBoundaryPressure();

    iterations_ = iteration;
};

const FieldVariable &ConjugateGradientPreconditoned::z() const {
    return z_;
};


double ConjugateGradientPreconditioned::z(int i, int j) const {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return z_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};


double &ConjugateGradientPreconditioned::z(int i, int j) {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return z_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};
