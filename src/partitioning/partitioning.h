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
    Partitioning(int ownRank, int worldSize, std::array<int, 2> nCellsGlobal);
    int column() const;
    int row() const;
    bool containsLeftBoundary() const;
    bool containsRightBoundary() const;
    bool containsBottomBoundary() const;
    bool containsTopBoundary() const;
    int ownRank() const;
    int leftRank() const;
    int rightRank() const;
    int bottomRank() const;
    int topRank() const;
    void sendToLeft(std::vector<double> &data) const;
    void sendToRight(std::vector<double> &data) const;
    void sendToBottom(std::vector<double> &data) const;
    void sendToTop(std::vector<double> &data) const;
    void recvFromLeft(std::vector<double> &data, int count) const;
    void recvFromRight(std::vector<double> &data, int count) const;
    void recvFromBottom(std::vector<double> &data, int count) const;
    void recvFromTop(std::vector<double> &data, int count) const;
    double globalSum(double localValue);
    const std::array<int, 2> nCells() const;
    const std::array<int, 2> nCellsGlobal() const;
private:
    int columnsBegin() const;
    int columnsEnd() const;
    int rowsBegin() const;
    int rowsEnd() const;
    int computeColumn(int rank) const;
    int computeRow(int rank) const;
    int computeRank(int i, int j) const;
    static void send(int destinationRank, std::vector<double> &data);
    static void recv(int sourceRank, std::vector<double> &data, int count);
    static double allReduce(double localValue, MPI_Op op);
    std::array<int, 2> nDomains_;
    int domainColumn_;
    int domainRow_;
    std::array<int, 2> nCells_;
    std::array<int, 2> nCellsGlobal_;
    int ownRank_;
    int worldSize_;
};