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
    q_ = std::make_unique<Array2D>(discretization_->pSize()); // Search direction qₖ
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
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];

    pGhostLayer();
    Array2D residual_ = Array2D(discretization_->pSize());
    Array2D Aq_ = Array2D(discretization_->pSize()); 

    int iteration = 0;
    // Initialization Loop
    for (int i = pIIntBegin; i < pIIntEnd; i++) {
        for (int j = pJIntBegin; j < pJIntEnd; j++) {
            double D2pDx2 = (discretization_->p(i-1, j) - 2 * discretization_->p(i, j) + discretization_->p(i+1, j)) / dx2;
            double D2pDy2 = (discretization_->p(i, j-1) - 2 * discretization_->p(i, j) + discretization_->p(i, j+1)) / dy2;
            // Calculate initial residuum (r₀)(i,j) = rhs(i,j) - (Ap)(i,j)
            residual_(i - pIBegin, j - pJBegin) = discretization_->rhs(i,j) - (D2pDx2 + D2pDy2);
            
            // set search direction to preconditioned defect q = z
            (*q_)(i - pIBegin, j - pJBegin) = residual_(i - pIBegin, j - pJBegin);
        }
    }

    double local_alpha = 0.0;  
    // Calculate initial alpha value
    for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
        for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
            local_alpha += residual_(i, j) * residual_(i, j); // α₀ = r₀ᵀ r₀
        }
    }
    double alpha = partitioning_->globalSum(local_alpha);

    do {

        qGhostLayer();

        // Calculate auxillary variable Aq
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                double D2qDx2 = ((*q_)(i-1, j) - 2 * (*q_)(i, j) + (*q_)(i+1, j)) / dx2;
                double D2qDy2 = ((*q_)(i, j-1) - 2 * (*q_)(i, j) + (*q_)(i, j+1)) / dy2;
                Aq_(i, j) = D2qDx2 + D2qDy2;
            }
        }

        double local_lambda = 0.0;        // λ = αₖ / qₖᵀAqₖ 
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                local_lambda += (*q_)(i, j) * Aq_(i, j);             // qₖᵀAqₖ
            }
        }
        double lambda = alpha / partitioning_->globalSum(local_lambda);
        iteration++;

        // Update variables in the search direction
        for (int i = pIIntBegin; i < pIIntEnd; i++) {
            for (int j = pJIntBegin; j < pJIntEnd; j++) {
                
                // pₖ₊₁ = pₖ + λ qₖ
                discretization_->p(i, j)+= + lambda * (*q_)(i - pIBegin, j - pJBegin); 
                
                // rₖ₊₁ = rₖ - λ Aqₖ
                residual_(i - pIBegin, j - pJBegin) = residual_(i - pIBegin, j - pJBegin) - lambda * Aq_(i - pIBegin ,j - pJBegin);
            }
        }

        double alphaold = alpha;
        local_alpha = 0.0;
        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                local_alpha += residual_(i, j) * residual_(i, j); // αₖ₊₁ = rₖ₊₁ᵀ rₖ₊₁
            }
        }
        alpha = partitioning_->globalSum(local_alpha);

        // βₖ₊₁ = αₖ₊₁ / αₖ
        double beta = alpha / alphaold;

        for (int i = pIIntBegin - pIBegin; i < pIIntEnd - pIBegin; i++) {
            for (int j = pJIntBegin - pJBegin; j < pJIntEnd - pJBegin; j++) {
                (*q_)(i, j) = residual_(i, j) + beta * (*q_)(i, j);             // qₖ₊₁ = rₖ₊₁ + β qₖ
            }
        }
        residual_norm2_ = alpha / N;

    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);

    pGhostLayer();

    iterations_ = iteration;
}

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
void ConjugateGradient::qGhostLayer() {
    const int pInteriorIBegin = discretization_->pInteriorIBegin();
    const int pInteriorIEnd = discretization_->pInteriorIEnd();
    const int pInteriorJBegin = discretization_->pInteriorJBegin();
    const int pInteriorJEnd = discretization_->pInteriorJEnd();
    const int pIBegin = discretization_->pIBegin();
    const int pJBegin = discretization_->pJBegin();
    const int pIEnd = discretization_->pIEnd();
    const int pJEnd = discretization_->pJEnd();
    const int pIIntBegin = discretization_->pInteriorIBegin();
    const int pJIntBegin = discretization_->pInteriorJBegin();
    const int pIIntEnd = discretization_->pInteriorIEnd();
    const int pJIntEnd = discretization_->pInteriorJEnd();

    int p_columnCount = pInteriorJEnd - pInteriorJBegin;
    int p_columnOffset = pInteriorJBegin;
    int p_rowCount = pInteriorIEnd - pInteriorIBegin;
    int p_rowOffset = pInteriorIBegin;

    MPI_Request request_p_rightColumn;
    MPI_Request request_p_leftColumn;
    MPI_Request request_p_topRow;
    MPI_Request request_p_bottomRow;
    std::vector<double> p_rightColumn(p_columnCount, 0);
    std::vector<double> p_leftColumn(p_columnCount, 0);
    std::vector<double> p_topRow(p_rowCount, 0);
    std::vector<double> p_bottomRow(p_rowCount, 0);
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        setQBoundaryValuesTop();
    }
    else {
        // send row on the top to top neighbour
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            p_topRow.at(i - p_rowOffset) = (*q_)(i - pIBegin, pInteriorJEnd - 1 - pJBegin);
        }
        partitioning_->isendToTop(p_topRow, request_p_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecvFromTop(p_topRow, p_rowCount, request_p_topRow);
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        setQBoundaryValuesBottom();
    }
    else {
        // send row on the bottom to bottom neighbour
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            p_bottomRow.at(i - p_rowOffset) = (*q_)(i - pIBegin, pInteriorJBegin - pJBegin);
        }
        partitioning_->isendToBottom(p_bottomRow, request_p_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecvFromBottom(p_bottomRow, p_rowCount, request_p_bottomRow);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        setQBoundaryValuesRight();
    }
    else {
        // send to column on the right to right neighbour
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            p_rightColumn.at(j - p_columnOffset) = (*q_)(pInteriorIEnd - 1 - pIBegin, j - pJBegin);
            //std::cout << "RANK " << partitioning_->ownRankNo() << " : right send no. " << j << " = " << p_leftColumn.at(j - p_columnOffset) << std::endl;
        }

        partitioning_->isendToRight(p_rightColumn, request_p_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecvFromRight(p_rightColumn, p_columnCount, request_p_rightColumn);
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        setQBoundaryValuesLeft();
    }
    else {
        // send to column on the left to left neighbour
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            p_leftColumn.at(j - p_columnOffset) = (*q_)(pInteriorIBegin - pIBegin, j - pJBegin);
        }
        partitioning_->isendToLeft(p_leftColumn, request_p_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecvFromLeft(p_leftColumn, p_columnCount, request_p_leftColumn);
    }

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        partitioning_->wait(request_p_topRow);
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            (*q_)(i - pIBegin, discretization_->pJEnd() - 1 - pJBegin) = p_topRow.at(i - p_rowOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        partitioning_->wait(request_p_bottomRow);
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            (*q_)(i - pIBegin, discretization_->pJBegin() - pJBegin) = p_bottomRow.at(i - p_rowOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        partitioning_->wait(request_p_rightColumn);
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            //std::cout << "RANK " << partitioning_->ownRankNo() << " : right recv no. " << j << " = " << p_rightColumn.at(j - p_columnOffset) << std::endl;
            (*q_)(discretization_->pIEnd() - 1 - pIBegin, j - pJBegin) = p_rightColumn.at(j - p_columnOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        partitioning_->wait(request_p_leftColumn);
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            (*q_)(discretization_->pIBegin() - pIBegin, j - pJBegin) = p_leftColumn.at(j - p_columnOffset);
        }
    }
}

void ConjugateGradient::setQBoundaryValuesBottom() {
    for (int i = 0; i < discretization_->pIEnd() - discretization_->pIBegin(); i++) {
        // copy values to bottom boundary
        (*q_)(i, 0) = (*q_)(i, 1);
    }
}

void ConjugateGradient::setQBoundaryValuesTop() {
    for (int i = 0; i < discretization_->pIEnd() - discretization_->pIBegin(); i++) {
        // copy values to top boundary
        (*q_)(i, discretization_->pJEnd() - 1 - discretization_->pJBegin()) = (*q_)(i, discretization_->pInteriorJEnd() - 1 - discretization_->pJBegin());
    }
}

void ConjugateGradient::setQBoundaryValuesLeft() {
    for (int j = 0; j < discretization_->pJEnd() - discretization_->pJBegin(); j++) {
        // copy values to left boundary
        (*q_)(0, j) = (*q_)(1, j);
    }
}

void ConjugateGradient::setQBoundaryValuesRight() {
    for (int j = 0; j < discretization_->pJEnd() - discretization_->pJBegin(); j++) {
        // copy values to right boundary
        (*q_)(discretization_->pIEnd() - 1 - discretization_->pIBegin(), j) = (*q_)(discretization_->pInteriorIEnd() - 1 - discretization_->pIBegin(), j);
    }
}