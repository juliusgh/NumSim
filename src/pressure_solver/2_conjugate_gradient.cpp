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
    const double eps2 = pow(epsilon_, 2);
    const int pIBegin = discretization_->pIBegin();
    const int pJBegin = discretization_->pJBegin();
    const int pIIntBegin = discretization_->pInteriorIBegin();
    const int pJIntBegin = discretization_->pInteriorJBegin();
    const int pIIntEnd = discretization_->pInteriorIEnd();
    const int pJIntEnd = discretization_->pInteriorJEnd();

    double MInv = ((dx2 * dy2) / (-2 * dy2 - 2 * dx2));
    MInv = 1.0;

    pGhostLayer();

    int iteration = 0;
    // Initialization Loop
    for (int i = pIIntBegin; i < pIIntEnd; i++) {
        for (int j = pJIntBegin; j < pJIntEnd; j++) {
            double D2pDx2 = (discretization_->p(i-1, j) - 2 * discretization_->p(i, j) + discretization_->p(i+1, j)) / dx2;
            double D2pDy2 = (discretization_->p(i, j-1) - 2 * discretization_->p(i, j) + discretization_->p(i, j+1)) / dy2;
            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            (*residual_)(i - pIBegin, j - pJBegin) = discretization_->rhs(i,j) - (D2pDx2 + D2pDy2);

            //std::cout << "RANK " << partitioning_->ownRankNo() << " : res(" << i - pIBegin << )" << local_alpha << std::endl;

            // precondition initial defect: "z = M^{-1} r" for M = diag(A)
            (*z_)(i - pIBegin, j - pJBegin) = MInv * (*residual_)(i - pIBegin, j - pJBegin);   
            
            // set search direction to preconditioned defect q = z
            (*q_)(i - pIBegin, j - pJBegin) = (*z_)(i - pIBegin, j - pJBegin);
        }
    }

    double local_alpha = 0.0;  
    // Calculate initial alpha value
    for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
        for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
            local_alpha += (*residual_)(i, j) * (*z_)(i, j); // α₀ = r₀ᵀ z₀
        }
    }
    double alpha = partitioning_->globalSum(local_alpha);
    std::cout << "RANK " << partitioning_->ownRankNo() << " : local_alpha = " << local_alpha << std::endl;
    std::cout << "RANK " << partitioning_->ownRankNo() << " : alpha = " << alpha << std::endl;
    getchar();

    do {

        pGhostLayer();

        // Calculate auxillary variable Aq
        std::cout << "RANK " << partitioning_->ownRankNo() << " : dx2 = " << dx2 << ", dy2 =" << dy2 << std::endl;
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                double D2qDx2 = ((*q_)(i-1, j) - 2 * (*q_)(i, j) + (*q_)(i+1, j)) / dx2;
                double D2qDy2 = ((*q_)(i, j-1) - 2 * (*q_)(i, j) + (*q_)(i, j+1)) / dy2;
                (*Aq_)(i, j) = D2qDx2 + D2qDy2;
                std::cout << "RANK " << partitioning_->ownRankNo() << " : Aq(" << i << "," << j << ")=" << (*Aq_)(i, j) << std::endl;
            }
        }

        double local_lambda = 0.0;        // λ = αₖ / qₖᵀAqₖ 
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                local_lambda += (*q_)(i, j) * (*Aq_)(i, j);             // qₖᵀAqₖ
            }
        }
        double lambda = alpha / partitioning_->globalSum(local_lambda);
        std::cout << "RANK " << partitioning_->ownRankNo() << " : local_lambda = " << local_lambda << std::endl;
        std::cout << "RANK " << partitioning_->ownRankNo() << " : lambda = " << lambda << std::endl;
        getchar();
        iteration++;

        // Update variables in the search direction
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                
                // pₖ₊₁ = pₖ + λ qₖ
                discretization_->p(i, j)= discretization_->p(i ,j) + lambda * (*q_)(i - pIBegin, j - pJBegin); 
                
                // rₖ₊₁ = rₖ - λ Aqₖ
                (*residual_)(i - pIBegin, j - pJBegin) = (*residual_)(i - pIBegin, j - pJBegin) - lambda * (*Aq_)(i - pIBegin ,j - pJBegin);
            }
        }

        // Preconditioned search direction
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                (*z_)(i, j) =  MInv * (*residual_)(i, j);
            }
        }

        double alphaold = alpha;
        local_alpha = 0.0;
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                local_alpha += (*residual_)(i, j) * (*z_)(i, j); // αₖ₊₁ = rₖ₊₁ᵀ zₖ₊₁
            }
        }
        alpha = partitioning_->globalSum(local_alpha);
        std::cout << "RANK " << partitioning_->ownRankNo() << " : local_alpha = " << local_alpha << std::endl;
        std::cout << "RANK " << partitioning_->ownRankNo() << " : alpha = " << alpha << std::endl;
        getchar();

        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                (*q_)(i, j) = (*z_)(i, j) + beta * (*q_)(i, j);             // qₖ₊₁ = zₖ₊₁ + β qₖ
            }
        }
        computeResidualNorm();

    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);

    pGhostLayer();

    iterations_ = iteration;
}


/**
 * Calculation of the global residual norm using the squared Euclidean norm
 */
void ConjugateGradient::computeResidualNorm() {
    double residual_norm2 = 0.0;
    const int pIBegin = discretization_->pIBegin();
    const int pJBegin = discretization_->pJBegin();
    const int pIIntBegin = discretization_->pInteriorIBegin();
    const int pJIntBegin = discretization_->pInteriorJBegin();
    const int pIIntEnd = discretization_->pInteriorIEnd();
    const int pJIntEnd = discretization_->pInteriorJEnd();
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
        for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
            residual_norm2 += pow((*residual_)(i, j), 2);
        }
    }

    double residual_norm2_global = partitioning_->globalSum(residual_norm2);
    residual_norm2_ = residual_norm2_global / N;
}

