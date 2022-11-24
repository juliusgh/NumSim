#pragma once

#include <array>
#include <mpi.h>
/**
 * encapsulate functionality corresponding to subdomain handling
 *
 * should be able to tell the own rank number and the rank numbers of neighbouring processes,
 * know whether the own subdomain touches one of the boundaries left, right, top or bottom of the global domain
 * and it should know the number of cells in the local staggered grid of the own subdomain.
 */

class Partitioning {
public:
    Partitioning(int ownRank, int worldSize, std::array<int, 2> nCells);
    bool containsLeftBoundary() const;
    bool containsRightBoundary() const;
    bool containsBottomBoundary() const;
    bool containsTopBoundary() const;
    int getOwnRank() const;
    int getLeftRank() const;
    int getRightRank() const;
    int getBottomRank() const;
    int getTopRank() const;
    const std::array<int, 2> getCellNumbers() const;
private:
    int columnsBegin() const;
    int columnsEnd() const;
    int rowsBegin() const;
    int rowsEnd() const;
    int computeColumn(int rank) const;
    int computeRow(int rank) const;
    int computeRank(int i, int j) const;
    std::array<int, 2> nDomains_;
    int domainColumn_;
    int domainRow_;
    const std::array<int, 2> nCells_;
    int ownRank_;
    int worldSize_;
};