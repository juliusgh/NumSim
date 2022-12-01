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
                std::cout << "RED: i = " << i << ", j = " << j << std::endl;
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / dx2;
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }

        pGhostLayerHorizontal();
        pGhostLayerVertical();

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

        pGhostLayerHorizontal();
        pGhostLayerVertical();

        computeResidualNorm();
    } while (residualNorm() > eps2 && iteration < maximumNumberOfIterations_);
    iterations_ = iteration;
}

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
void RedBlack::pGhostLayerHorizontal() {
    int columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int columnOffset = discretization_->pInteriorJBegin();

    // Even processes first send to left and right, then receive from left and right
    if (partitioning_->column() % 2) {
        if (partitioning_->ownPartitionContainsRightBoundary()) {
            setBoundaryValuesRight();
        }
        else {
            // send to column on the right to right neighbour
            std::vector<double> rightColumn;
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                rightColumn.push_back(discretization_->p(discretization_->pInteriorIEnd(), j));
            }
            partitioning_->sendToRight(rightColumn);

            // receive ghost layer column on the right from right neighbour
            partitioning_->recvFromRight(rightColumn, columnCount);
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pIEnd(), j) = rightColumn.at(j - columnOffset);
            }
        }
        if (partitioning_->ownPartitionContainsLeftBoundary()) {
            setBoundaryValuesLeft();
        }
        else {
            // send to column on the left to left neighbour
            std::vector<double> leftColumn;
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
            }
            partitioning_->sendToLeft(leftColumn);

            // receive ghost layer column on the left from left neighbour
            partitioning_->recvFromRight(leftColumn, columnCount);
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pIBegin(), j) = leftColumn.at(j - columnOffset);
            }
        }
    }

    // Odd processes first receive from left and right, then send to left and right
    else {
        if (partitioning_->ownPartitionContainsRightBoundary()) {
            setBoundaryValuesRight();
        }
        else {
            // receive ghost layer column on the right from right neighbour
            std::vector<double> rightColumn;
            partitioning_->recvFromRight(rightColumn, columnCount);
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pIEnd(), j) = rightColumn.at(j - columnOffset);
            }

            // send to column on the right to right neighbour
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                rightColumn.push_back(discretization_->p(discretization_->pInteriorIEnd(), j));
            }
            partitioning_->sendToRight(rightColumn);
        }
        if (partitioning_->ownPartitionContainsLeftBoundary()) {
            setBoundaryValuesLeft();
        }
        else {
            // receive ghost layer column on the right from right neighbour
            std::vector<double> leftColumn;
            partitioning_->recvFromLeft(leftColumn, columnCount);
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pIBegin(), j) = leftColumn.at(j - columnOffset);
            }

            // send to column on the left to left neighbour
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
            }
            partitioning_->sendToLeft(leftColumn);
        }
    }
}

/**
 *  Implementation of vertical communication of pressure values between neighbouring subdomains
 */
void RedBlack::pGhostLayerVertical() {
    int rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int rowOffset = discretization_->pInteriorIBegin();

    // Even processes first send to bottom and top, then receive from bottom and top
    if (partitioning_->row() % 2) {
        if (partitioning_->ownPartitionContainsTopBoundary()) {
            setBoundaryValuesTop();
        }
        else {
            // send row on the top to top neighbour
            std::vector<double> topRow;
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
            }
            partitioning_->sendToTop(topRow);

            // receive ghost layer row on the top from top neighbour
            partitioning_->recvFromTop(topRow, rowCount);
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                discretization_->p(i, discretization_->pJEnd()) = topRow.at(i - rowOffset);
            }
        }
        if (!partitioning_->ownPartitionContainsBottomBoundary()) {
            // send row on the bottom to bottom neighbour
            std::vector<double> bottomRow;
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                bottomRow.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
            }
            partitioning_->sendToBottom(bottomRow);

            // receive ghost layer row on the bottom from bottom neighbour
            partitioning_->recvFromTop(bottomRow, rowCount);
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                discretization_->p(i, discretization_->pJBegin()) = bottomRow.at(i - rowOffset);
            }
        }
    } else {
        if (partitioning_->ownPartitionContainsTopBoundary()) {
            setBoundaryValuesTop();
        }
        else {
            // receive row on the top from top neighbour
            std::vector<double> topRow;
            partitioning_->recvFromTop(topRow, rowCount);
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                discretization_->p(i, discretization_->pInteriorJEnd()) = topRow.at(i - rowOffset);
            }

            // send row on the top to top neighbour
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
            }
            partitioning_->sendToTop(topRow);
        }
        if (partitioning_->ownPartitionContainsBottomBoundary()) {
            setBoundaryValuesBottom();
        }
        else {
            // receive column on the right from right neighbour
            std::vector<double> bottomRow;
            partitioning_->recvFromTop(bottomRow, rowCount);
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                discretization_->p(i, discretization_->pInteriorJBegin()) = bottomRow.at(i - rowOffset);
            }

            // send row on the bottom to bottom neighbour
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                bottomRow.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
            }
            partitioning_->sendToBottom(bottomRow);
        }
    }
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