#include "pressure_solver/1_red_black.h"

/**
 * Implementation of the red-black solver, a parallelisized version of the Gauss-Seidel solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations
 */

RedBlack::RedBlack(std::shared_ptr<Discretization> discretization,
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
void RedBlack::solve() {
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
                //std::cout << "RED: i = " << i << ", j = " << j << std::endl;
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
                //std::cout << "BLACK: i = " << i << ", j = " << j << std::endl;
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
 * Calculation of the global residual norm using the squared Euclidean norm
 */
void RedBlack::computeResidualNorm() {
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

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
void RedBlack::pGhostLayer() {
    int columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int columnOffset = discretization_->pInteriorJBegin();
    int rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int rowOffset = discretization_->pInteriorIBegin();

    MPI_Request request_p_rightColumn;
    MPI_Request request_p_leftColumn;
    MPI_Request request_p_topRow;
    MPI_Request request_p_bottomRow;
    std::vector<double> rightColumn(columnCount, 0);
    std::vector<double> leftColumn(columnCount, 0);
    std::vector<double> topRow(rowCount, 0);
    std::vector<double> bottomRow(rowCount, 0);

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        //std::cout << "setBoundaryValuesRight" << std::endl;
        setBoundaryValuesRight();
    }
    else {
        //std::cout << "comm Right" << std::endl;
        // send to column on the right to right neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            rightColumn.push_back(discretization_->p(discretization_->pInteriorIEnd(), j));
        }

        partitioning_->isendToRight(rightColumn, request_p_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecvFromRight(rightColumn, columnCount, request_p_rightColumn);
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        //std::cout << "setBoundaryValuesLeft" << std::endl;
        setBoundaryValuesLeft();
    }
    else {
        //std::cout << "comm Left" << std::endl;
        // send to column on the left to left neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
        }
        partitioning_->isendToLeft(leftColumn, request_p_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecvFromLeft(leftColumn, columnCount, request_p_leftColumn);
    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        //std::cout << "setBoundaryValuesTop" << std::endl;
        setBoundaryValuesTop();
    }
    else {
        //std::cout << "comm Top" << std::endl;
        // send row on the top to top neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
        }
        partitioning_->isendToTop(topRow, request_p_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecvFromTop(topRow, rowCount, request_p_topRow);
    }
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //std::cout << "setBoundaryValuesBottom" << std::endl;
        setBoundaryValuesBottom();
    }
    else {
        //std::cout << "comm Bottom" << std::endl;
        // send row on the bottom to bottom neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            bottomRow.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
        }
        partitioning_->isendToBottom(bottomRow, request_p_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecvFromBottom(bottomRow, rowCount, request_p_bottomRow);
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        partitioning_->wait(request_p_rightColumn);
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            discretization_->p(discretization_->pIEnd(), j) = rightColumn.at(j - columnOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        partitioning_->wait(request_p_leftColumn);
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            discretization_->p(discretization_->pIBegin(), j) = leftColumn.at(j - columnOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        partitioning_->wait(request_p_topRow);
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = topRow.at(i - rowOffset);
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        partitioning_->wait(request_p_bottomRow);
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin()) = bottomRow.at(i - rowOffset);
        }
    }
}
