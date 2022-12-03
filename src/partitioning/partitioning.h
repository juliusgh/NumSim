#pragma once

#include <memory>
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

    void sendToLeft(std::vector<double> &data);
    void sendToRight(std::vector<double> &data);
    void sendToBottom(std::vector<double> &data);
    void sendToTop(std::vector<double> &data);
    void recvFromLeft(std::vector<double> &data, int count);
    void recvFromRight(std::vector<double> &data, int count);
    void recvFromBottom(std::vector<double> &data, int count);
    void recvFromTop(std::vector<double> &data, int count);
    void isendToLeft(std::vector<double> &data, MPI_Request &request);
    void isendToRight(std::vector<double> &data, MPI_Request &request);
    void isendToBottom(std::vector<double> &data, MPI_Request &request);
    void isendToTop(std::vector<double> &data, MPI_Request &request);
    void irecvFromLeft(std::vector<double> &data, int count, MPI_Request &request);
    void irecvFromRight(std::vector<double> &data, int count, MPI_Request &request);
    void irecvFromBottom(std::vector<double> &data, int count, MPI_Request &request);
    void irecvFromTop(std::vector<double> &data, int count, MPI_Request &request);
    void send(int destinationRank, std::vector<double> &data);
    void recv(int sourceRank, std::vector<double> &data, int count);
    void isend(int destinationRank, std::vector<double> &data, MPI_Request &request);
    void irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request);
    void wait(MPI_Request &request);
    double allReduce(double localValue, MPI_Op op);
    double globalSum(double localValue);
    double globalMax(double localValue);
    double globalMin(double localValue);
    const std::array<int, 2> nCellsLocal() const;
    const std::array<int, 2> nCellsGlobal() const;
    std::array<int, 2> nodeOffset() const;
    void log(const char* message);
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
    std::array<int, 2> nCellsLocal_;
    std::array<int, 2> nCellsGlobal_;
    std::array<int, 2> nodeOffset_;
    int ownRankNo_;
    int nRanks_;
};