#include "partitioning/partitioning.h"


Partitioning::Partitioning(int ownRank, int worldSize, std::array<int, 2> nCells)
: ownRank_(ownRank), worldSize_(worldSize), nCells_(nCells), nDomains_(std::array<int, 2>())
{
    MPI_Dims_create(worldSize_, 2, nDomains_.data());
    domainColumn_ = computeColumn(ownRank_);
    domainRow_ = computeRow(ownRank_);
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

int Partitioning::getLeftRank() const {
    if (containsLeftBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_ - 1, domainRow_);
}

int Partitioning::getRightRank() const {
    if (containsRightBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_ + 1, domainRow_);
}

int Partitioning::getBottomRank() const {
    if (containsBottomBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_, domainRow_ - 1);
}

int Partitioning::getTopRank() const {
    if (containsTopBoundary()) {
        // subdomain boundary is domain boundary
        return ownRank_;
    }
    return computeRank(domainColumn_, domainRow_ + 1);
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

const std::array<int, 2> Partitioning::getCellNumbers() const {
    return std::array<int, 2>();
}