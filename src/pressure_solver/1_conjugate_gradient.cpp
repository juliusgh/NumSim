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
        PressureSolver(discretization, epsilon, maximumNumberOfIterations),
        s_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth()),// Search direction sₖ
        residual_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth()),
        As_(discretization_->pSize(), {discretization_->dx() / 2.0, discretization_->dy() / 2.0}, discretization_->meshWidth())
        {
}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void ConjugateGradient::solve() {
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

    int iteration = 0;

    // Initialization Loop
    double alpha = 0.0;
    for (int i = pIIntBegin; i < pIIntEnd; i++) {
        for (int j = pJIntBegin; j < pJIntEnd; j++) {
            if (discretization_->marker(i, j) != FLUID) {
                continue;
            }

            double D2pDx2 =
                    (discretization_->p(i - 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / dx2;
            double D2pDy2 =
                    (discretization_->p(i, j - 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / dy2;

            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            residual(i, j) = discretization_->rhs(i, j) - (D2pDx2 + D2pDy2);

            // set search direction to preconditioned defect s = z
            s(i, j) = residual(i, j);

            alpha += pow(residual(i, j), 2);
        }
    }

    do {
        applyBoundarySearchDirection();

        double lambda = 0.0;
        // Calculate auxiliary variable As
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i, j) != FLUID) {
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
        alpha = 0.0;
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i, j) != FLUID) {
                    continue;
                }
                // pₖ₊₁ = pₖ + λ sₖ
                discretization_->p(i, j) += lambda * s(i, j);

                // rₖ₊₁ = rₖ - λ Asₖ
                residual(i, j) -= lambda * As(i, j);

                alpha += pow(residual(i, j), 2);

            }
        }
        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i, j) != FLUID) {
                    continue;
                }
                s(i, j) = residual(i, j) + beta * s(i, j);             // sₖ₊₁ = rₖ₊₁ + β sₖ
            }
        }
        residual_norm2_ = alpha / N;

    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);

    discretization_->applyBoundaryPressure();

    iterations_ = iteration;
};

void ConjugateGradient::applyBoundarySearchDirection() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            switch (discretization_->marker(i, j)) {
                case OBSTACLE_LEFT:
                    s(i, j) = s(i - 1, j);
                    break;
                case OBSTACLE_RIGHT:
                    s(i, j) = s(i + 1, j);
                    break;
                case OBSTACLE_TOP:
                    s(i, j) = s(i,j + 1);
                    break;
                case OBSTACLE_BOTTOM:
                    s(i, j) = s(i, j - 1);
                    break;
                case OBSTACLE_LEFT_TOP:
                    s(i, j) = (s(i - 1, j) + s(i,j + 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    s(i, j) = (s(i + 1, j) + s(i,j + 1)) / 2.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    s(i, j) = (s(i - 1, j) + s(i,j - 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    s(i, j) = (s(i + 1, j) + s(i,j - 1)) / 2.0;
                    break;
                default:
                    break;
            }
        }
    }

    // set boundary values for p at bottom and top side (lower priority)
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        // set boundary values at bottom side
        switch (discretization_->marker(i, discretization_->pJBegin())) {
            case INFLOW:
            case NOSLIP:
                s(i, discretization_->pJBegin()) = s(i, discretization_->pInteriorJBegin());
                break;
            case OUTFLOW:
                s(i, discretization_->uJBegin()) = -s(i, discretization_->uInteriorJBegin());
                break;
            default:
                break;
        }
        // set boundary values for p at top side
        switch (discretization_->marker(i, discretization_->pJEnd() - 1)) {
            case NOSLIP:
            case INFLOW:
                s(i, discretization_->pJEnd() - 1) = s(i, discretization_->pInteriorJEnd() - 1);
                break;
            case OUTFLOW:
                s(i, discretization_->uJEnd() - 1) = -s(i, discretization_->uInteriorJEnd() - 1);
                break;
            default:
                break;
        }
    }

    // set boundary values for p at left and right side (higher priority)
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
        // set boundary values for p at left side
        switch (discretization_->marker(discretization_->pIBegin(), j)) {
            case NOSLIP:
            case INFLOW:
                s(discretization_->pIBegin(), j) = s(discretization_->pInteriorIBegin(), j);
                break;
            case OUTFLOW:
                s(discretization_->pIBegin(), j) = -s(discretization_->uInteriorIBegin(), j);
                break;
            default:
                break;
        }
        // set boundary values for p at right side
        switch (discretization_->marker(discretization_->pIEnd() - 1, j)) {
            case NOSLIP:
            case INFLOW:
                s(discretization_->pIEnd() - 1, j) = s(discretization_->pInteriorIEnd() - 1, j);
                break;
            case OUTFLOW:
                s(discretization_->pIEnd() - 1, j) = -s(discretization_->uInteriorIEnd() - 1, j);
                break;
            default:
                break;
        }
    }
};

const FieldVariable &ConjugateGradient::s() const {
    return s_;
};


double ConjugateGradient::s(int i, int j) const {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return s_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};


double &ConjugateGradient::s(int i, int j) {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return s_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};

const FieldVariable &ConjugateGradient::residual() const {
    return residual_;
};


double ConjugateGradient::residual(int i, int j) const {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return residual_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};


double &ConjugateGradient::residual(int i, int j) {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return residual_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};

const FieldVariable &ConjugateGradient::As() const {
    return As_;
};


double ConjugateGradient::As(int i, int j) const {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return As_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};


double &ConjugateGradient::As(int i, int j) {
#ifndef NDEBUG
    assert((discretization_->pIBegin() <= i) && (i <= discretization_->pIEnd()));
    assert((discretization_->pJBegin() <= j) && (j <= discretization_->pJEnd()));
#endif
    return As_(i - discretization_->pIBegin(), j - discretization_->pJBegin());
};
