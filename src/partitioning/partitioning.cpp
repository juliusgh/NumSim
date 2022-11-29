#include <vector>
#include "partitioning/partitioning.h"


Partitioning::Partitioning(int ownRank, int worldSize, std::array<int, 2> nCellsGlobal)
: ownRank_(ownRank), worldSize_(worldSize), nCellsGlobal_(nCellsGlobal), nDomains_(std::array<int, 2>())
{
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

int Partitioning::column() const {
    return computeColumn(ownRank());
}

int Partitioning::row() const {
    return computeRow(ownRank());
}

int Partitioning::columnsBegin() const {
    return 1;
}

int Partitioning::columnsEnd() const {
    return nDomains_[0];
}

int Partitioning::rowsBegin() const {
    return 1;
}

int Partitioning::rowsEnd() const {
    return nDomains_[1];
}

bool Partitioning::containsLeftBoundary() const {
    return domainColumn_ == columnsBegin();
}

bool Partitioning::containsRightBoundary() const {
    return domainColumn_ == columnsEnd();
}

bool Partitioning::containsBottomBoundary() const {
    return domainRow_ == rowsBegin();
}

bool Partitioning::containsTopBoundary() const {
    return domainColumn_ == rowsEnd();
}

int Partitioning::ownRank() const {
    return ownRank_;
}

int Partitioning::leftRank() const {
    if (containsLeftBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_ - 1, domainRow_);
}

int Partitioning::rightRank() const {
    if (containsRightBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_ + 1, domainRow_);
}

int Partitioning::bottomRank() const {
    if (containsBottomBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_, domainRow_ - 1);
}

int Partitioning::topRank() const {
    if (containsTopBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_, domainRow_ + 1);
}

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

void Partitioning::sendToLeft(std::vector<double> &data) const {
    send(leftRank(), data);
}

void Partitioning::sendToRight(std::vector<double> &data) const {
    send(rightRank(), data);
}

void Partitioning::sendToBottom(std::vector<double> &data) const {
    send(bottomRank(), data);
}

void Partitioning::sendToTop(std::vector<double> &data) const {
    send(topRank(), data);
}

void Partitioning::recvFromLeft(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

void Partitioning::recvFromRight(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

void Partitioning::recvFromBottom(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

void Partitioning::recvFromTop(std::vector<double> &data, int count) const {
    recv(leftRank(), data, count);
}

int Partitioning::computeColumn(int rank) const {
    return rank % nDomains_[0];
}

int Partitioning::computeRow(int rank) const {
    return (int)((rank - 1) / nDomains_[0]) + 1;
}

int Partitioning::computeRank(int column, int row) const {
    return column + nDomains_[0] * row;
}

const std::array<int, 2> Partitioning::nCells() const {
    return nCells_;
}

const std::array<int, 2> Partitioning::nCellsGlobal() const {
    return nCellsGlobal_;
}

double Partitioning::globalSum(double localValue) {
    return allReduce(localValue, MPI_SUM);
}

double Partitioning::globalMax(double localValue) {
    return allReduce(localValue, MPI_MAX);
}

double Partitioning::globalMin(double localValue) {
    return allReduce(localValue, MPI_MIN);
}

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
