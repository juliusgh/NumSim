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
partitioning_(partitioning)
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
    int iteration = 0;
    do {
        iteration++;

        // red half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j += 2) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                std::cout << "RED: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayer();

        // black half step
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j += 2) {
            int iStart = discretization_->pInteriorIBegin() + (j - discretization_->pInteriorJBegin()) % 2;
            for (int i = iStart; i < discretization_->pInteriorIEnd(); i += 2) {
                std::cout << "BLACK: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayer();

        computeResidualNorm();
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
}

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
void ConjugateGradient::pGhostLayer() {
    int columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int columnOffset = discretization_->pInteriorJBegin();
    int rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int rowOffset = discretization_->pInteriorIBegin();

    MPI_Request request;
    std::vector<double> rightColumn;
    std::vector<double> leftColumn;
    std::vector<double> topRow;
    std::vector<double> bottomRow;

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        setBoundaryValuesRight();
    }
    else {
        // send to column on the right to right neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            rightColumn.push_back(discretization_->p(discretization_->pInteriorIEnd(), j));
        }

        partitioning_->isend(partitioning_->rightNeighbourRankNo(), rightColumn, request);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecv(partitioning_->rightNeighbourRankNo(), rightColumn, columnCount, request);
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        setBoundaryValuesLeft();
    }
    else {
        // send to column on the left to left neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
        }
        partitioning_->isend(partitioning_->leftNeighbourRankNo(), leftColumn, request);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecv(partitioning_->leftNeighbourRankNo(), leftColumn, columnCount, request);
    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        setBoundaryValuesTop();
    }
    else {
        // send row on the top to top neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
        }
        partitioning_->isend(partitioning_->topNeighbourRankNo(), topRow, request);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecv(partitioning_->topNeighbourRankNo(), topRow, rowCount, request);
    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        setBoundaryValuesBottom();
    }
    else {
        // send row on the bottom to bottom neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            bottomRow.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
        }
        partitioning_->isend(partitioning_->bottomNeighbourRankNo(), bottomRow, request);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecv(partitioning_->bottomNeighbourRankNo(), bottomRow, rowCount, request);
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin()) = bottomRow.at(i - rowOffset);
        }
    }

    partitioning_->wait(request);

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            discretization_->p(discretization_->pIEnd(), j) = rightColumn.at(j - columnOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            discretization_->p(discretization_->pIBegin(), j) = leftColumn.at(j - columnOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = topRow.at(i - rowOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin()) = bottomRow.at(i - rowOffset);
        }
    }
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
            double pxx = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2;
            double pyy = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2;
            residual_norm2 += pow(pxx + pyy - discretization_->rhs(i, j), 2);
        }
    }

    double residual_norm2_global = partitioning_->globalSum(residual_norm2);
    residual_norm2_ = residual_norm2_global / N;
}