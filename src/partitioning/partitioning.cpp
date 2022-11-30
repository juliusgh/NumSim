#include <vector>
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
: nCellsGlobal_(nCellsGlobal), nDomains_(std::array<int, 2>())
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize_);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);

    // Partition the domain
    MPI_Dims_create(worldSize_, 2, nDomains_.data());
    domainColumn_ = computeColumn(ownRank_);
    domainRow_ = computeRow(ownRank_);
    // compute nCells_, nCellsGlobal_
    nCells_ = std::array<int, 2>{ nCellsGlobal_[0] / nDomains_[0], nCellsGlobal_[1] / nDomains_[1] };
    int cellsColumnRemainder = nCellsGlobal_[0] % nDomains_[0];
    int cellsRowRemainder = nCellsGlobal_[1] % nDomains_[1];
    if (domainColumn_ <= cellsColumnRemainder)
        nCells_[0]++;
    if (domainRow_ <= cellsRowRemainder)
        nCells_[1]++;
}
/**
 * get column length of process
 * @return column length of process
 */
int Partitioning::column() const {
    return computeColumn(ownRank());
}

/**
 * get column row of process
 * @return row length of process
 */
int Partitioning::row() const {
    return computeRow(ownRank());
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
bool Partitioning::containsLeftBoundary() const {
    return domainColumn_ == columnsBegin();
}
/**
 * check whether current subdomain touches global right boundary
 * @return  whether current subdomain touches global right boundary
 */
bool Partitioning::containsRightBoundary() const {
    return domainColumn_ == columnsEnd();
}
/**
 * check whether current subdomain touches global bottom boundary
 * @return  whether current subdomain touches global bottom boundary
 */
bool Partitioning::containsBottomBoundary() const {
    return domainRow_ == rowsBegin();
}
/**
 * check whether current subdomain touches global top boundary
 * @return  whether current subdomain touches global top boundary
 */
bool Partitioning::containsTopBoundary() const {
    return domainColumn_ == rowsEnd();
}

/**
 * get own rank
 * @return own rank
 */
int Partitioning::ownRank() const {
    return ownRank_;
}

/**
 * get rank of left neighbour process.
 * In the case that the current boundary touches the left boundary, its own rank is returned.
 * @return rank of left neighbour process
 */
int Partitioning::leftRank() const {
    if (containsLeftBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_ - 1, domainRow_);
}
/**
 * get rank of right neighbour process.
 * In the case that the current boundary touches the right boundary, its own rank is returned.
 * @return rank of right neighbour process
 */
int Partitioning::rightRank() const {
    if (containsRightBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_ + 1, domainRow_);
}
/**
 * get rank of bottom neighbour process.
 * In the case that the current boundary touches the bottom boundary, its own rank is returned.
 * @return rank of bottom neighbour process
 */
int Partitioning::bottomRank() const {
    if (containsBottomBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_, domainRow_ - 1);
}

/**
 * get rank of top neighbour process.
 * In the case that the current boundary touches the top boundary, its own rank is returned.
 * @return rank of top neighbour process
 */
int Partitioning::topRank() const {
    if (containsTopBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_, domainRow_ + 1);
}
/**
 * Implementation of call of the MPI-send command
 * @param destinationRank rank of the process, data is send to
 * @param data data to be send
 */
void Partitioning::send(int destinationRank, std::vector<double> &data) {
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
 * method used to send information to the subdomain left of the current subdomain
 * @param data information to be send to the left
 */
void Partitioning::sendToLeft(std::vector<double> &data) const {
    send(leftRank(), data);
}
/**
 * method used to send information to the subdomain right of the current subdomain
 * @param data information to be send to the right
 */
void Partitioning::sendToRight(std::vector<double> &data) const {
    send(rightRank(), data);
}

/**
 * method used to send information to the subdomain below the current subdomain
 * @param data information to be send down
 */
void Partitioning::sendToBottom(std::vector<double> &data) const {
    send(bottomRank(), data);
}

/**
 * method used to send information to the subdomain above the current subdomain
 * @param data information to be send up
 */
void Partitioning::sendToTop(std::vector<double> &data) const {
    send(topRank(), data);
}

/**
 * method used to receive information from the subdomain left of the current subdomain
 * @param data information to be received from the left
 */
void Partitioning::recvFromLeft(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

/**
 * method used to receive information from the subdomain right of the current subdomain
 * @param data information to be received from the right
 */
void Partitioning::recvFromRight(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

/**
 * method used to receive information from the subdomain below of the current subdomain
 * @param data information to be received from below
 */
void Partitioning::recvFromBottom(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

/**
 * method used to receive information from the subdomain above of the current subdomain
 * @param data information to be received from the top
 */
void Partitioning::recvFromTop(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

/**
 * get the column position of the subdomain
 * @param rank
 * @return column position of the subdomain
 */
int Partitioning::computeColumn(int rank) const {
    return rank % nDomains_[0];
}

/**
 * get the row position of the subdomain
 * @param rank
 * @return row position of the subdomain
 */
int Partitioning::computeRow(int rank) const {
    return (int)((rank - 1) / nDomains_[0]) + 1;
}

/**
 * get the rank of a subdomain based on its position in the global grid
 * @param column
 * @param row
 * @return
 */
int Partitioning::computeRank(int column, int row) const {
    return column + nDomains_[0] * row;
}

/**
 * get local number of cells in current subdomain
 * @return number of cells in current subdomain
 */
const std::array<int, 2> Partitioning::nCells() const {
    return nCells_;
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
    double globalValue;
    MPI_Allreduce(&localValue,
                  &globalValue,
                  1,
                  MPI_DOUBLE,
                  op,
                  MPI_COMM_WORLD
    );
    return globalValue;
}
