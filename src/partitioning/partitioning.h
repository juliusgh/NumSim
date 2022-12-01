#pragma once

#include <array>
#include <mpi.h>

/**
 * Partitioning encapsulate functionality corresponding to subdomain handling.
 * This includes the ability to tell the own rank and the ranks of the neighbouring processes.
 * It can check whether a process lies on a  global grid boundary.
 * Additionally, the number of cells in a local staggered grid can be queried.
 *
 */

class Partitioning {
public:
    Partitioning(std::array<int, 2> nCellsGlobal);
    int column() const;
    int row() const;

    //! if the own partition has part of the bottom boundary of the whole domain
    bool ownPartitionContainsBottomBoundary() const;

    //! if the own partition has part of the top boundary of the whole domain
    //! used in OutputWriterParaviewParallel
    bool ownPartitionContainsTopBoundary() const;

    //! if the own partition has part of the left boundary of the whole domain
    bool ownPartitionContainsLeftBoundary() const;

    //! if the own partition has part of the right boundary of the whole domain
    //! used in OutputWriterParaviewParallel
    bool ownPartitionContainsRightBoundary() const;

    //! number of MPI ranks
    int nRanks() const;

    //! get the own MPI rank no
    //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
    int ownRankNo() const;

    //! get the rank no of the left neighbouring rank
    int leftNeighbourRankNo() const;

    //! get the rank no of the right neighbouring rank
    int rightNeighbourRankNo() const;

    //! get the rank no of the bottom neighbouring rank
    int bottomNeighbourRankNo() const;

    //! get the rank no of the top neighbouring rank
    int topNeighbourRankNo() const;

    void sendToLeft(std::vector<double> &data) const;
    void sendToRight(std::vector<double> &data) const;
    void sendToBottom(std::vector<double> &data) const;
    void sendToTop(std::vector<double> &data) const;
    void recvFromLeft(std::vector<double> &data, int count) const;
    void recvFromRight(std::vector<double> &data, int count) const;
    void recvFromBottom(std::vector<double> &data, int count) const;
    void recvFromTop(std::vector<double> &data, int count) const;
    double globalSum(double localValue);
    double globalMax(double localValue);
    double globalMin(double localValue);
    const std::array<int, 2> nCellsLocal() const;
    const std::array<int, 2> nCellsGlobal() const;
    std::array<int, 2> nodeOffset() const;
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
    std::array<int, 2> nCellsLocal_;
    std::array<int, 2> nCellsGlobal_;
    std::array<int, 2> nodeOffset_;
    int ownRankNo_;
    int nRanks_;
};