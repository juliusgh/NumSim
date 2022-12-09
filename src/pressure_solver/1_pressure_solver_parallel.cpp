#include "pressure_solver/1_pressure_solver_parallel.h"

/**
 *  Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 * @param discretization: pointer to distrectization implementation
 * @param epsilon: error tolerance
 * @param maximumNumberOfIterations: maximal number of iterations before ending the iteration
 */

PressureSolverParallel::PressureSolverParallel(std::shared_ptr<Discretization> discretization,
                               double epsilon,
                               int maximumNumberOfIterations,
                               std::shared_ptr<Partitioning> partitioning) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations),
partitioning_(partitioning)
{

}

/**
 *  Implementation of horizontal communication of pressure values between neighbouring subdomains
 */
void PressureSolverParallel::pGhostLayer() {
    const int pInteriorIBegin = discretization_->pInteriorIBegin();
    const int pInteriorIEnd = discretization_->pInteriorIEnd();
    const int pInteriorJBegin = discretization_->pInteriorJBegin();
    const int pInteriorJEnd = discretization_->pInteriorJEnd();

    int p_columnCount = pInteriorJEnd - pInteriorJBegin;
    int p_columnOffset = pInteriorJBegin;
    int p_rowCount = pInteriorIEnd - pInteriorIBegin;
    int p_rowOffset = pInteriorIBegin;

    MPI_Request request_p_rightColumn; // send to right - receive from left
    MPI_Request request_p_leftColumn; // send to left - receive from right
    MPI_Request request_p_topRow; // send to top - receive from bottom
    MPI_Request request_p_bottomRow; // send to bottom - receive from top
    std::vector<double> p_rightColumn(p_columnCount, 0);
    std::vector<double> p_leftColumn(p_columnCount, 0);
    std::vector<double> p_topRow(p_rowCount, 0);
    std::vector<double> p_bottomRow(p_rowCount, 0);
    /*
    * pressure p communication: send to and receive from subdomain directly above 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches upper domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        setBoundaryValuesTop();
    }
    else {
        // send row on the top to top neighbour
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            p_topRow.at(i - p_rowOffset) = discretization_->p(i, pInteriorJEnd - 1);
        }
        partitioning_->isendToTop(p_topRow, request_p_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecvFromTop(p_topRow, p_rowCount, request_p_topRow);
    }
    /*
    * pressure p communication: send to and receive from subdomain directly below 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches lower domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //std::cout << "setBoundaryValuesBottom" << std::endl;
        setBoundaryValuesBottom();
    }
    else {
        // send row on the bottom to bottom neighbour
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            p_bottomRow.at(i - p_rowOffset) = discretization_->p(i, pInteriorJBegin);
        }
        partitioning_->isendToBottom(p_bottomRow, request_p_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecvFromBottom(p_bottomRow, p_rowCount, request_p_bottomRow);
    }
    /*
    * pressure p communication: send to and receive from subdomain directly right of 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches right domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsRightBoundary()) {
        setBoundaryValuesRight();
    }
    else {
        // send to column on the right to right neighbour
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            p_rightColumn.at(j - p_columnOffset) = discretization_->p(pInteriorIEnd - 1, j);
        }

        partitioning_->isendToRight(p_rightColumn, request_p_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecvFromRight(p_rightColumn, p_columnCount, request_p_rightColumn);
    }
    /*
    * pressure p communication: send to and receive from subdomain directly left 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches left domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        setBoundaryValuesLeft();
    }
    else {
        // send to column on the left to left neighbour
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            p_leftColumn.at(j - p_columnOffset) = discretization_->p(pInteriorIBegin, j);
        }
        partitioning_->isendToLeft(p_leftColumn, request_p_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecvFromLeft(p_leftColumn, p_columnCount, request_p_leftColumn);
    }
   /* 
   * Set subdomain ghost layer at top for pressure
   */
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        partitioning_->wait(request_p_topRow);
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            discretization_->p(i, discretization_->pJEnd() - 1) = p_topRow.at(i - p_rowOffset);
        }
    }
   
   /* 
   * Set subdomain ghost layer at bottom for pressure
   */
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        partitioning_->wait(request_p_bottomRow);
        for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
            discretization_->p(i, discretization_->pJBegin()) = p_bottomRow.at(i - p_rowOffset);
        }
    }
   /* 
   * Set subdomain ghost layer at right for pressure
   */
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        partitioning_->wait(request_p_rightColumn);
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            discretization_->p(discretization_->pIEnd() - 1, j) = p_rightColumn.at(j - p_columnOffset);
        }
    }
   /* 
   * Set subdomain ghost layer at left for pressure
   */
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        partitioning_->wait(request_p_leftColumn);
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            discretization_->p(discretization_->pIBegin(), j) = p_leftColumn.at(j - p_columnOffset);
        }
    }
}

/**
 * Calculation of the global residual norm using the squared Euclidean norm
 */
void PressureSolverParallel::computeResidualNorm() {
    const int pInteriorIBegin = discretization_->pInteriorIBegin();
    const int pInteriorIEnd = discretization_->pInteriorIEnd();
    const int pInteriorJBegin = discretization_->pInteriorJBegin();
    const int pInteriorJEnd = discretization_->pInteriorJEnd();

    double residual_norm2 = 0.0;
    const double dx2 = pow(discretization_->dx(),2);
    const double dy2 = pow(discretization_->dy(),2);
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    for (int i = pInteriorIBegin; i < pInteriorIEnd; i++) {
        for (int j = pInteriorJBegin; j < pInteriorJEnd; j++) {
            double pxx = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2;
            double pyy = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2;
            residual_norm2 += pow(pxx + pyy - discretization_->rhs(i, j), 2);
        }
    }
    double residual_norm2_global = partitioning_->globalSum(residual_norm2);
    residual_norm2_ = residual_norm2_global / N;
}
