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
         std::shared_ptr<Partitioning> partitioning,
         std::shared_ptr<Array2D> residual) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations),
partitioning_(partitioning),
residual_(residual)
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
    std::unique_ptr<Array2D> q = std::make_unique<Array2D>(discretization_->pSize()); // Search direction qₖ
    std::unique_ptr<Array2D> Aq = std::make_unique<Array2D>(discretization_->pSize()); //precomputed matrix-vector product between the system matrix and the search direction Aqₖ

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

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
void RedBlack::pGhostLayer() {
    int p_columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int p_columnOffset = discretization_->pInteriorJBegin();
    int p_rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int p_rowOffset = discretization_->pInteriorIBegin();

    MPI_Request request_p_rightColumn;
    MPI_Request request_p_leftColumn;
    MPI_Request request_p_topRow;
    MPI_Request request_p_bottomRow;
    std::vector<double> p_rightColumn(p_columnCount, 0);
    std::vector<double> p_leftColumn(p_columnCount, 0);
    std::vector<double> p_topRow(p_rowCount, 0);
    std::vector<double> p_bottomRow(p_rowCount, 0);
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        //std::cout << "setBoundaryValuesTop" << std::endl;
        setBoundaryValuesTop();
    }
    else {
        //std::cout << "comm Top" << std::endl;
        // send row on the top to top neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            p_topRow.at(i - p_rowOffset) = discretization_->p(i, discretization_->pInteriorJEnd() - 1);
        }
        partitioning_->isendToTop(p_topRow, request_p_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecvFromTop(p_topRow, p_rowCount, request_p_topRow);
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //std::cout << "setBoundaryValuesBottom" << std::endl;
        setBoundaryValuesBottom();
    }
    else {
        //std::cout << "comm Bottom" << std::endl;
        // send row on the bottom to bottom neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            p_bottomRow.at(i - p_rowOffset) = discretization_->p(i, discretization_->pInteriorJBegin());
        }
        partitioning_->isendToBottom(p_bottomRow, request_p_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecvFromBottom(p_bottomRow, p_rowCount, request_p_bottomRow);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        //std::cout << "setBoundaryValuesRight" << std::endl;
        setBoundaryValuesRight();
    }
    else {
        //std::cout << "comm Right" << std::endl;
        // send to column on the right to right neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            p_rightColumn.at(j - p_columnOffset) = discretization_->p(discretization_->pInteriorIEnd() - 1, j);
            //std::cout << "RANK " << partitioning_->ownRankNo() << " : right send no. " << j << " = " << p_leftColumn.at(j - p_columnOffset) << std::endl;
        }

        partitioning_->isendToRight(p_rightColumn, request_p_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecvFromRight(p_rightColumn, p_columnCount, request_p_rightColumn);
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        //std::cout << "setBoundaryValuesLeft" << std::endl;
        setBoundaryValuesLeft();
    }
    else {
        //std::cout << "comm Left" << std::endl;
        // send to column on the left to left neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            p_leftColumn.at(j - p_columnOffset) = discretization_->p(discretization_->pInteriorIBegin(), j);
            //std::cout << "RANK " << partitioning_->ownRankNo() << " : left send no. " << j << " = " << p_leftColumn.at(j - p_columnOffset) << std::endl;
        }
        partitioning_->isendToLeft(p_leftColumn, request_p_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecvFromLeft(p_leftColumn, p_columnCount, request_p_leftColumn);
    }

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        partitioning_->wait(request_p_topRow);
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd() - 1) = p_topRow.at(i - p_rowOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        partitioning_->wait(request_p_bottomRow);
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin()) = p_bottomRow.at(i - p_rowOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        partitioning_->wait(request_p_rightColumn);
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            //std::cout << "RANK " << partitioning_->ownRankNo() << " : right recv no. " << j << " = " << p_rightColumn.at(j - p_columnOffset) << std::endl;
            discretization_->p(discretization_->pIEnd() - 1, j) = p_rightColumn.at(j - p_columnOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        partitioning_->wait(request_p_leftColumn);
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            //std::cout << "RANK " << partitioning_->ownRankNo() << " : left recv no. " << j << " = " << p_leftColumn.at(j - p_columnOffset) << std::endl;
            discretization_->p(discretization_->pIBegin(), j) = p_leftColumn.at(j - p_columnOffset);
        }
    }
}
