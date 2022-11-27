#include "pressure_solver/1_red_black.h"

/**
 * Successive over-relaxation solver for solving a linear system of equations.
 * @param discretization
 * @param epsilon
 * @param maximumNumberOfIterations
 * @param omega
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

void RedBlack::pGhostLayerHorizontal() {
    int columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int columnOffset = discretization_->pInteriorJBegin();
    // 1. step: Even processes send to left and right, odd processes receive from left and right
    if (partitioning_->ownRank() % 2) {
        if (!partitioning_->containsRightBoundary()) {
            // send to column on the right to right neighbour
            int destinationRank = partitioning_->rightRank();
            std::vector<double> rightColumn;
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                rightColumn.push_back(discretization_->p(discretization_->pInteriorIEnd(), j));
            }
            MPI_Send(
                    rightColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    destinationRank,
                    0,
                    MPI_COMM_WORLD);
        }
        if (!partitioning_->containsLeftBoundary()) {
            // send to column on the left to left neighbour
            int destinationRank = partitioning_->leftRank();
            std::vector<double> leftColumn;
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
            }
            MPI_Send(
                    leftColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    destinationRank,
                    0,
                    MPI_COMM_WORLD);
        }
    }

    // 2. step: Even processes receive from left and right, odd processes send to left and right
    else {
        if (!partitioning_->containsRightBoundary()) {
            // receive column on the right from right neighbour
            int sourceRank = partitioning_->rightRank();
            std::vector<double> rightColumn;
            MPI_Recv(
                    rightColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    sourceRank,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
            );
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pInteriorIEnd(), j) = rightColumn.at(j - columnOffset);
            }
        }
        if (!partitioning_->containsLeftBoundary()) {
            // receive column on the right from right neighbour
            int sourceRank = partitioning_->leftRank();
            std::vector<double> leftColumn;
            MPI_Recv(
                    leftColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    sourceRank,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
            );
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pInteriorIBegin(), j) = leftColumn.at(j - columnOffset);
            }
        }
    }
}

void RedBlack::pGhostLayerVertical() {
    int rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int rowOffset = discretization_->pInteriorIBegin();
    // 1. step: Even processes send to left and right, odd processes receive from left and right
    if (partitioning_->ownRank() % 2) {
        if (!partitioning_->containsTopBoundary()) {
            // send to column on the top to top neighbour
            int destinationRank = partitioning_->topRank();
            std::vector<double> topRow;
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
            }
            MPI_Send(
                    topRow.data(),
                    rowCount,
                    MPI_DOUBLE,
                    destinationRank,
                    0,
                    MPI_COMM_WORLD);
        }
        if (!partitioning_->containsBottomBoundary()) {
            // send to column on the bottom to bottom neighbour
            int destinationRank = partitioning_->bottomRank();
            std::vector<double> bottomColumn;
            for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
                bottomColumn.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
            }
            MPI_Send(
                    leftColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    destinationRank,
                    0,
                    MPI_COMM_WORLD);
        }
    } else {
        if (!partitioning_->containsRightBoundary()) {
            // receive column on the right from right neighbour
            int sourceRank = partitioning_->rightRank();
            std::vector<double> rightColumn;
            MPI_Recv(
                    rightColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    sourceRank,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
            );
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pInteriorIEnd(), j) = rightColumn.at(j - columnOffset);
            }
        }
        if (!partitioning_->containsLeftBoundary()) {
            // receive column on the right from right neighbour
            int sourceRank = partitioning_->leftRank();
            std::vector<double> leftColumn;
            MPI_Recv(
                    leftColumn.data(),
                    columnCount,
                    MPI_DOUBLE,
                    sourceRank,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
            );
            for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
                discretization_->p(discretization_->pInteriorIBegin(), j) = leftColumn.at(j - columnOffset);
            }
        }
    }
}

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

    double residual_norm2_global;
    MPI_Allreduce(&residual_norm2, &residual_norm2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    residual_norm2_ = residual_norm2_global / N;
}