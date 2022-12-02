#include "pressure_solver/0_pressure_solver.h"

/**
 *  Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 * @param discretization: pointer to distrectization implementation
 * @param epsilon: error tolerance
 * @param maximumNumberOfIterations: maximal number of iterations before ending the iteration
 */

PressureSolver::PressureSolver(std::shared_ptr <Discretization> discretization,
                               double epsilon,
                               int maximumNumberOfIterations) :
        discretization_(discretization),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations) {

}

/**
 * Set pressure solver boundary values
 */
void PressureSolver::setBoundaryValues() {
    // copy values to bottom and top boundary (lower priority)
    setBoundaryValuesBottom();
    setBoundaryValuesTop();

    // copy values to left and right boundary (higher priority)
    setBoundaryValuesLeft();
    setBoundaryValuesRight();
}

void PressureSolver::setBoundaryValuesBottom() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        // copy values to bottom boundary
        discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pInteriorJBegin());
    }
}

void PressureSolver::setBoundaryValuesTop() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        // copy values to top boundary
        discretization_->p(i, discretization_->pJEnd() - 1) = discretization_->p(i, discretization_->pInteriorJEnd() - 1);
    }
}

void PressureSolver::setBoundaryValuesLeft() {
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
        // copy values to left boundary
        discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pInteriorIBegin(), j);
    }
}

void PressureSolver::setBoundaryValuesRight() {
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
        // copy values to right boundary
        discretization_->p(discretization_->pIEnd() - 1, j) = discretization_->p(discretization_->pInteriorIEnd() - 1, j);
    }
}

/**
 * Compute squared Euclidean norm to measure pressure solver convergence performance
 * ||X||_2^2 = x^2 + y^2 for X = (x,y)^T
 */
void PressureSolver::computeResidualNorm() {
    double residual_norm2 = 0.0;
    const double dx2 = pow(discretization_->dx(),2);
    const double dy2 = pow(discretization_->dy(),2);
    const int N = discretization_->nCells()[0] * discretization_->nCells()[1];
    for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            double pxx = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2;
            double pyy = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2;
            residual_norm2 += pow(pxx + pyy - discretization_->rhs(i, j), 2);
        }
    }
    residual_norm2_ = residual_norm2 / N;
};

double PressureSolver::residualNorm() const {
    return residual_norm2_;
}

int PressureSolver::iterations() const {
    return iterations_;
}

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
virtual void RedBlack::pGhostLayer() {
    int columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int columnOffset = discretization_->pInteriorJBegin();
    int rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int rowOffset = discretization_->pInteriorIBegin();

    MPI_Request request_p_rightColumn;
    MPI_Request request_p_leftColumn;
    MPI_Request request_p_topRow;
    MPI_Request request_p_bottomRow;
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

        partitioning_->isend(partitioning_->rightNeighbourRankNo(), rightColumn, request_p_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecv(partitioning_->rightNeighbourRankNo(), rightColumn, columnCount, request_p_rightColumn);
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        setBoundaryValuesLeft();
    }
    else {
        // send to column on the left to left neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
        }
        partitioning_->isend(partitioning_->leftNeighbourRankNo(), leftColumn, request_p_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecv(partitioning_->leftNeighbourRankNo(), leftColumn, columnCount, request_p_leftColumn);
    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        setBoundaryValuesTop();
    }
    else {
        // send row on the top to top neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
        }
        partitioning_->isend(partitioning_->topNeighbourRankNo(), topRow, request_p_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecv(partitioning_->topNeighbourRankNo(), topRow, rowCount, request_p_topRow);
    }
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        setBoundaryValuesBottom();
    }
    else {
        // send row on the bottom to bottom neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            bottomRow.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
        }
        partitioning_->isend(partitioning_->bottomNeighbourRankNo(), bottomRow, request_p_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecv(partitioning_->bottomNeighbourRankNo(), bottomRow, rowCount, request_p_bottomRow);
    }

    partitioning_->wait(request_p_rightColumn);
    partitioning_->wait(request_p_leftColumn);
    partitioning_->wait(request_p_topRow);
    partitioning_->wait(request_p_bottomRow);

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