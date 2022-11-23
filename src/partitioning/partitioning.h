#pragma once

#include <array>
/**
 * encapsulate functionality corresponding to subdomain handling
 *
 * should be able to tell the own rank number and the rank numbers of neighbouring processes,
 * know whether the own subdomain touches one of the boundaries left, right, top or bottom of the global domain
 * and it should know the number of cells in the local staggered grid of the own subdomain.
 */

class Partitioning {
public:
    int getLeftRank(int rank) const;
    int getRightRank(int rank) const;
    int getTopRank(int rank) const;
    int getBottomRank(int rank) const;
protected:
    int rank2i(int rank) const;
    int rank2j(int rank) const;
    int ij2rank(int i, int j) const;
    const std::array<int, 2> nDomains_;

};