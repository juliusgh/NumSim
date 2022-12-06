#include "pressure_solver/2_conjugate_gradient.h"

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
PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning)
{
    residual_ = std::make_unique<Array2D>(discretization_->pSize());
    q_ = std::make_unique<Array2D>(discretization_->pSize()); // Search direction qₖ
    z_ = std::make_unique<Array2D>(discretization_->pSize()); // preconditioned search direction zₖ
    Aq_ = std::make_unique<Array2D>(discretization_->pSize()); //precomputed matrix-vector product between the system matrix and the search direction Aqₖ

}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void ConjugateGradient::solve() {
    const double dx2 = pow(discretization_->dx(), 2);
    const double dy2 = pow(discretization_->dy(), 2);
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = pow(epsilon_, 2);
    const int pIBegin = discretization_->pIBegin();
    const int pJBegin = discretization_->pJBegin();

    int iteration = 0;
    // Initialization Loop
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            double D2pDx2 = (discretization_->p(i-1, j) - 2 * discretization_->p(i, j) + discretization_->p(i+1, j)) / dx2;
            double D2pDy2 = (discretization_->p(i, j-1) - 2 * discretization_->p(i, j) + discretization_->p(i, j+1)) / dy2;
            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            (*residual_)(i - pIBegin, j - pJBegin) = discretization_->rhs(i,j) - (D2pDx2 + D2pDy2);
            
            // precondition initial defect: "z = M^{-1} r" for M = diag(A)
            (*z_)(i - pIBegin, j - pJBegin) = (*residual_)(i - pIBegin, j - pJBegin) * ((dx2 * dy2) / (-2 * dy2 - 2 * dx2));   
            
            // set search direction to preconditioned defect q = z
            (*q_)(i - pIBegin, j - pJBegin) = (*z_)(i - pIBegin, j - pJBegin);
        }
    }


        
    std::cout << "\n" <<std::endl;
    for (int i = discretization_->pIBegin() - pIBegin; i < discretization_->pIEnd() - pIBegin; i++) {
        for (int j = discretization_->pJBegin() - pJBegin; j < discretization_->pJEnd() - pJBegin; j++) {
            std::cout <<"("<< i <<","<< j<< ")" <<(*residual_)(i ,j) << "  ";
        }
        std::cout << " " << std::endl;
    }
    

    double local_alpha = 0.0;  
    // Calculate initial alpha value
    for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
        for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
            local_alpha += (*residual_)(i, j) * (*q_)(i, j); // α₀ = r₀ᵀ q₀
        }
    }
    double alpha = partitioning_->globalSum(local_alpha);

    do {
        // Calculate auxillary variable Aq
        for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
            for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
                double D2qDx2 = ((*q_)(i-1, j) - 2 * (*q_)(i, j) + (*q_)(i+1, j)) / dx2;
                double D2qDy2 = ((*q_)(i, j-1) - 2 * (*q_)(i, j) + (*q_)(i, j+1)) / dy2;
                (*Aq_)(i, j) = D2qDx2 + D2qDy2;
            }
        }

        double local_lambda = 0.0;        // λ = αₖ / qₖᵀAqₖ 
        for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
            for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
                local_lambda += (*q_)(i, j) * (*Aq_)(i, j);             // qₖᵀAqₖ
            }
        }
        double lambda = alpha / partitioning_->globalSum(local_lambda);
        iteration++;

        // Update variables in the search direction
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                
                // pₖ₊₁ = pₖ + α qₖ
                discretization_->p(i, j)= discretization_->p(i ,j) + lambda * (*q_)(i - pIBegin, j - pJBegin); 
                
                // rₖ₊₁ = rₖ - α Aqₖ
                (*residual_)(i - pIBegin, j - pJBegin) = (*residual_)(i - pIBegin, j - pJBegin) - lambda * (*Aq_)(i - pIBegin ,j - pJBegin);
            }
        }

        // Preconditioned search direction
        for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
            for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
                (*z_)(i, j) = (*residual_)(i, j) * ((dx2 * dy2) / (-2 * dy2 - 2 * dx2));
            }
        }

        double alphaold = alpha;
        local_alpha = 0.0;
        for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
            for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
                local_alpha += (*residual_)(i, j) * (*z_)(i, j); // αₖ₊₁ = rₖ₊₁ᵀ zₖ₊₁
            }
        }
        alpha = partitioning_->globalSum(local_alpha);

        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
            for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
                (*q_)(i, j) = (*z_)(i, j) + beta * (*q_)(i, j);             // qₖ₊₁ = zₖ₊₁ + β qₖ
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
    const int pIBegin = discretization_->pIBegin();
    const int pJBegin = discretization_->pJBegin();
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    for (int i = discretization_->pInteriorIBegin() - pIBegin; i < discretization_->pInteriorIEnd() - pIBegin; i++) {
        for (int j = discretization_->pInteriorJBegin() - pJBegin; j < discretization_->pInteriorJEnd() - pJBegin; j++) {
            residual_norm2 += pow((*residual_)(i, j), 2);
        }
    }

    double residual_norm2_global = partitioning_->globalSum(residual_norm2);
    residual_norm2_ = residual_norm2_global / N;
}

