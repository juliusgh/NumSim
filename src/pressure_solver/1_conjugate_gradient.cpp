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
        s_(discretization_->pSize(),{0,0},{0,0})// Search direction qₖ
        {
}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void ConjugateGradient::solve() {
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

    s_.setToZero();
    applyBoundarySearchDirection();
    Array2D residual_ = Array2D(discretization_->pSize());
    Array2D Aq_ = Array2D(discretization_->pSize());

    int iteration = 0;

    // Initialization Loop
    double alpha = 0.0;
    for (int i = pIIntBegin; i < pIIntEnd; i++) {
        for (int j = pJIntBegin; j < pJIntEnd; j++) {
            if (discretization_->marker(i,j) != FLUID) {
                continue;
            }
            double D2pDx2 =
                    (discretization_->p(i - 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / dx2;
            double D2pDy2 =
                    (discretization_->p(i, j - 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / dy2;
            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            residual_(i, j ) = discretization_->rhs(i, j) - (D2pDx2 + D2pDy2);

            // set search direction to preconditioned defect q = z
            s(i, j) = residual_(i , j);


            alpha += pow(residual_(i, j), 2);
        }
    }

    do {
        //std::cout << "p" << std::endl;
        //discretization_->p().print();
        //std::cout << "r" << std::endl;
        //residual_.print();
        //std::cout << "q" << std::endl;
        //q.print();

        double lambda = 0.0;
        // Calculate auxillary variable Aq
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i,j) != FLUID) {
                    continue;
                }
                double D2qDx2 = (s(i - 1, j) - 2 * s(i, j) + s(i + 1, j)) / dx2;
                double D2qDy2 = (s(i, j - 1) - 2 * s(i, j) + s(i, j + 1)) / dy2;
                Aq_(i, j) = D2qDx2 + D2qDy2;

                // qₖᵀAqₖ
                lambda += s(i, j) * Aq_(i, j);
            }
        }

        // λ = αₖ / qₖᵀAqₖ
        lambda = alpha / lambda;
        iteration++;

        double alphaold = alpha;
        alpha = 0.0;
        // Update variables in the search direction
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i,j) != FLUID) {
                    continue;
                }
                // pₖ₊₁ = pₖ + λ qₖ
                discretization_->p(i, j) += lambda * s(i, j);

                // rₖ₊₁ = rₖ - λ Aqₖ
                residual_(i, j) -= lambda * Aq_(i, j);

                // αₖ₊₁ = rₖ₊₁ᵀ rₖ₊₁
                alpha += pow(residual_(i, j), 2);
            }
        }

        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                if (discretization_->marker(i,j) != FLUID) {
                    continue;
                }
                s(i, j) = residual_(i, j) + beta * s(i, j);             // qₖ₊₁ = rₖ₊₁ + β qₖ
            }
        }
        residual_norm2_ = alpha / N;
        discretization_->applyBoundaryPressure();

    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);

    discretization_->applyBoundaryPressure();

    iterations_ = iteration;
}

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
