#include "partitioning/partitioning.h"


Partitioning::Partitioning(std::array<int, 2> nCells, int worldSize)
: nCells_(nCells), worldSize_(worldSize), nDomains_(std::array<int, 2>())
{
    // TODO: determine nDomains_;

}

int Partitioning::getLeftRank(int rank) const {
    int i = rank2i(rank);
    int j = rank2j(rank);
    if (i == 1) {
        // subdomain boundary is domain boundary
        return rank;
    }
    return ij2rank(i - 1, j);
}

int Partitioning::getRightRank(int rank) const {
    int i = rank2i(rank);
    int j = rank2j(rank);
    if (i == nDomains_[0]) {
        // subdomain boundary is domain boundary
        return rank;
    }
    return ij2rank(i + 1, j);
}

int Partitioning::getTopRank(int rank) const {
    int i = rank2i(rank);
    int j = rank2j(rank);
    if (j == nDomains_[1]) {
        // subdomain boundary is domain boundary
        return rank;
    }
    return ij2rank(i, j + 1);
}

int Partitioning::getBottomRank(int rank) const {
    int i = rank2i(rank);
    int j = rank2j(rank);
    if (j == 0) {
        // subdomain boundary is domain boundary
        return rank;
    }
    return ij2rank(i, j - 1);
}

int Partitioning::rank2i(int rank) const {
    return rank % nDomains_[0];
}

int Partitioning::rank2j(int rank) const {
    return (int)((rank - 1) / nDomains_[0]) + 1;
}

int Partitioning::ij2rank(int i, int j) const {
    return i + nDomains_[0] * j;
}

const std::array<int, 2> Partitioning::getCellNumbers(int rank) const {
    return std::array<int, 2>();
}