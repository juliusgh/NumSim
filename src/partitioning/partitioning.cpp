#include <vector>
#include <iostream>
#include <cmath>
#include "partitioning/partitioning.h"

/**
 * Partitioning encapsulate functionality corresponding to subdomain handling.
 * This includes the ability to tell the own rank and the ranks of the neighbouring processes.
 * It can check whether a process lies on a  global grid boundary.
 * Additionally, the number of cells in a local staggered grid can be queried.
 *
 * @param nCellsGlobal total number of cells in the global domain
 */
Partitioning::Partitioning(std::array<int, 2> nCellsGlobal)
        : nCellsGlobal_(nCellsGlobal), nDomains_(std::array<int, 2>()) {
#ifndef NPARALLEL
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);

    // Partition the domain
    /*
     * The following routine is only optimal for a quadratic domain.
     * However, for some reason we achieved better results than by using
     * our optimization procedure partitionDomainEqual.
     */
    MPI_Dims_create(nRanks_, 2, nDomains_.data());
    // partitionDomainEqual(nRanks_);



#ifndef NDEBUG
    if (ownRankNo() == 0) {
        std::cout << "RANK " << ownRankNo() << " | nDomains[0]: " << nDomains_[0] << ", nDomains[1]: " << nDomains_[1]
                  << std::endl;
        std::cout << "RANK " << ownRankNo() << " | number of processes: " << nRanks_ << std::endl;
    }
#endif
#else
    nRanks_ = 1;
    ownRankNo_ = 0;
    nDomains_ = {1, 1};
#endif

    domainColumn_ = computeColumn(ownRankNo_);
    domainRow_ = computeRow(ownRankNo_);

    // Compute nCells_, nCellsGlobal_
    nCellsLocal_ = std::array<int, 2>{nCellsGlobal_[0] / nDomains_[0], nCellsGlobal_[1] / nDomains_[1]};
    int cellsColumnRemainder = nCellsGlobal_[0] % nDomains_[0];
    int cellsRowRemainder = nCellsGlobal_[1] % nDomains_[1];

    int columnOffset = 0;
    for (int i = columnsBegin(); i < domainColumn_; i++) {
        columnOffset += nCellsLocal_[0];
        if (i <= cellsColumnRemainder)
            columnOffset++;
    }

    int rowOffset = 0;
    for (int j = rowsBegin(); j < domainRow_; j++) {
        rowOffset += nCellsLocal_[1];
        if (j <= cellsRowRemainder)
            rowOffset++;
    }
    if (domainColumn_ <= cellsColumnRemainder)
        nCellsLocal_[0]++;
    if (domainRow_ <= cellsRowRemainder)
        nCellsLocal_[1]++;

    nodeOffset_ = {columnOffset, rowOffset};
#ifndef NDEBUG
    std::cout << "RANK " << ownRankNo() << " | nCellsLocal_[0]: " << nCellsLocal_[0] << ", nCellsLocal_[1]: "
              << nCellsLocal_[1] << std::endl;
    std::cout << "RANK " << ownRankNo() << " | nodeOffset_[0]: " << nodeOffset_[0] << ", nodeOffset_[1]: "
              << nodeOffset_[1] << std::endl;
#endif
}

/**
 * this method partitions the global domain in nRanks subdomains  
 * @param nRanks number of processes
*/
void Partitioning::partitionDomain(int nRanks) {
    int optimalX = nRanks;
    int optimalY = 1;

    // Ratio between the length of the subdomains in x and y direction should be near 1
    double best_cost =
            (2 * (nCellsGlobal_[0] / optimalX) + 2 * (nCellsGlobal_[1] / optimalY)) / (nCellsGlobal_[0] / optimalX) *
            (nCellsGlobal_[1] / optimalY);

    // Iterate over all possible combinations of partitionings
    for (int testY = 1; testY < nRanks + 1; testY++) {

        if ((nRanks % testY) == 0) {
            int testX = nRanks / testY;

            double cost =
                    (2 * (nCellsGlobal_[0] / testX) + 2 * (nCellsGlobal_[1] / testY)) / (nCellsGlobal_[0] / testX) *
                    (nCellsGlobal_[1] / testY);

            if (cost < best_cost) {
                best_cost = cost;
                optimalY = testY;
                optimalX = testX;
            }
        }
    }
    nDomains_ = {optimalX, optimalY};
}

/**
 * this method partitions the global domain in nRanks subdomains  
 * @param nRanks number of processes
*/
void Partitioning::partitionDomainEqual(int nRanks) {
    int optimalX = nRanks;
    int optimalY = 1;

    // Ratio between the length of the subdomains in x and y direction should be near 1
    double best_cost = fabs((nCellsGlobal_[0] / optimalX) / (nCellsGlobal_[1] / optimalY) - 1.0);

    // Iterrate over all possible combinations of partitionings
    for (int testY = 1; testY < nRanks + 1; testY++) {

        if ((nRanks % testY) == 0) {
            int testX = nRanks / testY;

            double cost = (nCellsGlobal_[0] / testX) / (nCellsGlobal_[1] / testY);

            if (fabs(cost - 1.0) < best_cost) {
                best_cost = cost;
                optimalY = testY;
                optimalX = testX;
            }
        }
    }

    nDomains_ = {optimalX, optimalY};
}


/**
 * get column index of subdomain
 * @return column length of process
 */
int Partitioning::column() const {
    return computeColumn(ownRankNo());
}

/**
 * get row index of subdomain
 * @return row length of process
 */
int Partitioning::row() const {
    return computeRow(ownRankNo());
}

/**
 * get discrete subdomain starting index in x-direction
 * @return discrete subdomain starting index in x-direction
 */
int Partitioning::columnsBegin() const {
    return 1;
}

/**
 * get last discrete subdomain index in x-direction
 * @return last discrete subdomain index in x-direction
 */
int Partitioning::columnsEnd() const {
    return nDomains_[0];
}

/**
 * get discrete subdomain starting index in y-direction
 * @return discrete subdomain starting index in y-direction
 */
int Partitioning::rowsBegin() const {
    return 1;
}

/**
 * get last discrete subdomain index in y-direction
 * @return last discrete subdomain index in y-direction
 */
int Partitioning::rowsEnd() const {
    return nDomains_[1];
}

/**
 * check whether current subdomain touches global left boundary
 * @return  whether current subdomain touches global left boundary
 */
bool Partitioning::ownPartitionContainsLeftBoundary() const {
    return domainColumn_ == columnsBegin();
}

/**
 * check whether current subdomain touches global right boundary
 * @return  whether current subdomain touches global right boundary
 */
bool Partitioning::ownPartitionContainsRightBoundary() const {
    return domainColumn_ == columnsEnd();
}

/**
 * check whether current subdomain touches global bottom boundary
 * @return  whether current subdomain touches global bottom boundary
 */
bool Partitioning::ownPartitionContainsBottomBoundary() const {
    return domainRow_ == rowsBegin();
}

/**
 * check whether current subdomain touches global top boundary
 * @return  whether current subdomain touches global top boundary
 */
bool Partitioning::ownPartitionContainsTopBoundary() const {
    return domainRow_ == rowsEnd();
}

/**
 * get own rank
 * @return own rank
 */
int Partitioning::ownRankNo() const {
    return ownRankNo_;
}

/**
 * get number of ranks (world size)
 * @return world size
 */
int Partitioning::nRanks() const {
    return nRanks_;
}

/**
 * get rank of left neighbour process.
 * In the case that the current boundary touches the left boundary, its own rank is returned.
 * @return rank of left neighbour process
 */
int Partitioning::leftNeighbourRankNo() const {
    if (ownPartitionContainsLeftBoundary()) {
        // subdomain boundary is domain boundary
        return ownRankNo_;
    }
    return computeRank(domainColumn_ - 1, domainRow_);
}

/**
 * get rank of right neighbour process.
 * In the case that the current boundary touches the right boundary, its own rank is returned.
 * @return rank of right neighbour process
 */
int Partitioning::rightNeighbourRankNo() const {
    if (ownPartitionContainsRightBoundary()) {
        // subdomain boundary is domain boundary
        return ownRankNo_;
    }
    return computeRank(domainColumn_ + 1, domainRow_);
}

/**
 * get rank of bottom neighbour process.
 * In the case that the current boundary touches the bottom boundary, its own rank is returned.
 * @return rank of bottom neighbour process
 */
int Partitioning::bottomNeighbourRankNo() const {
    if (ownPartitionContainsBottomBoundary()) {
        // subdomain boundary is domain boundary
        return ownRankNo_;
    }
    return computeRank(domainColumn_, domainRow_ - 1);
}

/**
 * get rank of top neighbour process.
 * In the case that the current boundary touches the top boundary, its own rank is returned.
 * @return rank of top neighbour process
 */
int Partitioning::topNeighbourRankNo() const {
    if (ownPartitionContainsTopBoundary()) {
        // subdomain boundary is domain boundary
        return ownRankNo_;
    }
    return computeRank(domainColumn_, domainRow_ + 1);
}

/**
 * Implementation of call of the MPI-send command
 * @param destinationRank rank of the process, data is send to
 * @param data data to be send
 */
void Partitioning::send(int destinationRank, std::vector<double> &data) {
    std::cout << "own Rank: " << ownRankNo() << ", destination Rank: " << destinationRank << std::endl;
    MPI_Send(
            data.data(),
            data.size(),
            MPI_DOUBLE,
            destinationRank,
            0,
            MPI_COMM_WORLD
    );
}

/**
 * Implementation of call of the MPI-receive command
 * @param sourceRank rank of the process, data is received from
 * @param data data to be received
 * @param count size of data to be received
 */
void Partitioning::recv(int sourceRank, std::vector<double> &data, int count) {
    MPI_Recv(
            data.data(),
            count,
            MPI_DOUBLE,
            sourceRank,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
    );
}

/**
 * Implementation of call of the MPI-send command
 * @param destinationRank rank of the process, data is send to
 * @param data data to be send
 */
void Partitioning::isend(int destinationRank, std::vector<double> &data, MPI_Request &request) {
    MPI_Isend(
            data.data(),
            data.size(),
            MPI_DOUBLE,
            destinationRank,
            0,
            MPI_COMM_WORLD,
            &request
    );
}

/**
 * Implementation of call of the MPI-receive command
 * @param sourceRank rank of the process, data is received from
 * @param data data to be received
 * @param count size of data to be received
 */
void Partitioning::irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request) {
    MPI_Irecv(
            data.data(),
            count,
            MPI_DOUBLE,
            sourceRank,
            0,
            MPI_COMM_WORLD,
            &request
    );
}

/**
 * Implementation of call of the MPI-wait command
 * @param request request the program should wait for before continuing
*/
void Partitioning::wait(MPI_Request &request) {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
}

/**
 * method used to send information to the subdomain left of the current subdomain
 * @param data information to be send to the left
 */
void Partitioning::sendToLeft(std::vector<double> &data) {
    send(leftNeighbourRankNo(), data);
}

/**
 * method used to send information to the subdomain right of the current subdomain
 * @param data information to be send to the right
 */
void Partitioning::sendToRight(std::vector<double> &data) {
    send(rightNeighbourRankNo(), data);
}

/**
 * method used to send information to the subdomain below the current subdomain
 * @param data information to be send down
 */
void Partitioning::sendToBottom(std::vector<double> &data) {
    send(bottomNeighbourRankNo(), data);
}

/**
 * method used to send information to the subdomain above the current subdomain
 * @param data information to be send up
 */
void Partitioning::sendToTop(std::vector<double> &data) {
    send(topNeighbourRankNo(), data);
}

/**
 * method used to receive information from the subdomain left of the current subdomain
 * @param data information to be received from the left
 */
void Partitioning::recvFromLeft(std::vector<double> &data, int count) {
    recv(leftNeighbourRankNo(), data, count);
}

/**
 * method used to receive information from the subdomain right of the current subdomain
 * @param data information to be received from the right
 */
void Partitioning::recvFromRight(std::vector<double> &data, int count) {
    recv(rightNeighbourRankNo(), data, count);
}

/**
 * method used to receive information from the subdomain below of the current subdomain
 * @param data information to be received from below
 */
void Partitioning::recvFromBottom(std::vector<double> &data, int count) {
    recv(bottomNeighbourRankNo(), data, count);
}

/**
 * method used to receive information from the subdomain above of the current subdomain
 * @param data information to be received from the top
 */
void Partitioning::recvFromTop(std::vector<double> &data, int count) {
    recv(topNeighbourRankNo(), data, count);
}


/**
 * method used to send information to the subdomain left of the current subdomain
 * @param data information to be send to the left
 */
void Partitioning::isendToLeft(std::vector<double> &data, MPI_Request &request) {
    isend(leftNeighbourRankNo(), data, request);
}

/**
 * method used to send information to the subdomain right of the current subdomain
 * @param data information to be send to the right
 */
void Partitioning::isendToRight(std::vector<double> &data, MPI_Request &request) {
    isend(rightNeighbourRankNo(), data, request);
}

/**
 * method used to send information to the subdomain below the current subdomain
 * @param data information to be send down
 */
void Partitioning::isendToBottom(std::vector<double> &data, MPI_Request &request) {
    isend(bottomNeighbourRankNo(), data, request);
}

/**
 * method used to send information to the subdomain above the current subdomain
 * @param data information to be send up
 */
void Partitioning::isendToTop(std::vector<double> &data, MPI_Request &request) {
    isend(topNeighbourRankNo(), data, request);
}

/**
 * method used to receive information from the subdomain left of the current subdomain
 * @param data information to be received from the left
 */
void Partitioning::irecvFromLeft(std::vector<double> &data, int count, MPI_Request &request) {
    irecv(leftNeighbourRankNo(), data, count, request);
}

/**
 * method used to receive information from the subdomain right of the current subdomain
 * @param data information to be received from the right
 */
void Partitioning::irecvFromRight(std::vector<double> &data, int count, MPI_Request &request) {
    irecv(rightNeighbourRankNo(), data, count, request);
}

/**
 * method used to receive information from the subdomain below of the current subdomain
 * @param data information to be received from below
 */
void Partitioning::irecvFromBottom(std::vector<double> &data, int count, MPI_Request &request) {
    irecv(bottomNeighbourRankNo(), data, count, request);
}

/**
 * method used to receive information from the subdomain above of the current subdomain
 * @param data information to be received from the top
 */
void Partitioning::irecvFromTop(std::vector<double> &data, int count, MPI_Request &request) {
    irecv(topNeighbourRankNo(), data, count, request);
}

/**
 * get the column position of the subdomain
 * @param rank
 * @return column position of the subdomain
 */
int Partitioning::computeColumn(int rank) const {
    return rank % nDomains_[0] + 1;
}

/**
 * get the row position of the subdomain
 * @param rank unique number of process
 * @return row position of the subdomain
 */
int Partitioning::computeRow(int rank) const {
    return (int) (rank / nDomains_[0]) + 1;
}

/**
 * get the rank of a subdomain based on its position in the global grid
 * @param column positions of subdomain in horizontal direction 
 *               on the grid compared to other subdomains
 * @param row positions of subdomain in vertical direction 
*             on the grid compared to other subdomains
 * @return subdomain rank
 */
int Partitioning::computeRank(int column, int row) const {
    return (column - 1) + nDomains_[0] * (row - 1);
}

/**
 * get local number of cells in current subdomain
 * @return number of cells in current subdomain
 */
const std::array<int, 2> Partitioning::nCellsLocal() const {
    return nCellsLocal_;
}

/**
 * get global number of cells in domain
 * @return number of cells in global domain
 */
const std::array<int, 2> Partitioning::nCellsGlobal() const {
    return nCellsGlobal_;
}

/**
 * sum local values over multiple subdomains
 * @param localValue local values on subdomains
 * @return sum over values of multiple subdomains
 */
double Partitioning::globalSum(double localValue) {
    return allReduce(localValue, MPI_SUM);
}

/**
 *  get maximum of a value over multiple subdomains
 * @param localValue  local values
 * @return global maximum
 */
double Partitioning::globalMax(double localValue) {
    return allReduce(localValue, MPI_MAX);
}

/**
 *  get minimum of a value over multiple subdomains
 * @param localValue  local values
 * @return global minimum
 */
double Partitioning::globalMin(double localValue) {
    return allReduce(localValue, MPI_MIN);
}

/**
 * Implementation of call of MPI-allReduce command,
 * which combines values from all processes and distributes the result back to all processes
 * @param localValue values on subdomains
 * @param op MPI operation to be performed in allReduce
 * @return globalValue combined value
 */
double Partitioning::allReduce(double localValue, MPI_Op op) {
#ifdef NPARALLEL
    return localValue;
#else
    double globalValue = 0.0;
    MPI_Allreduce(&localValue,
                  &globalValue,
                  1,
                  MPI_DOUBLE,
                  op,
                  MPI_COMM_WORLD
    );
    return globalValue;
#endif
}

/**
 * get the offset values for counting local nodes in x and y direction.
 * (i_local,j_local) + nodeOffset = (i_global,j_global)
 * used in OutputWriterParaviewParallel
 * @return offset of current subdomain indices compared to the global position
 */
std::array<int, 2> Partitioning::nodeOffset() const {
    return nodeOffset_;
}

/**
 * Method for debugging to print out what a process is currently doing
*/
void Partitioning::log(const char *message) {
    std::cout << "RANK " << ownRankNo() << " : " << message << std::endl;
}
