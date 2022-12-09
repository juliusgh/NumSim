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
    
    /**
    * get column index of subdomain
    * @return column length of process
    */
    int column() const;
    
    /**
     * get row index of subdomain
     * @return row length of process
     */    
    int row() const;

    /**
     * check whether current subdomain touches global bottom boundary
     * @return  whether current subdomain touches global bottom boundary
     */    
    bool ownPartitionContainsBottomBoundary() const;

    /**
     * check whether current subdomain touches global top boundary
     * @return  whether current subdomain touches global top boundary
    */
    bool ownPartitionContainsTopBoundary() const;

   /**
    * check whether current subdomain touches global left boundary
    * @return  whether current subdomain touches global left boundary
   */
   bool ownPartitionContainsLeftBoundary() const;

    /**
     * check whether current subdomain touches global right boundary
     * @return  whether current subdomain touches global right boundary
    */
    bool ownPartitionContainsRightBoundary() const;

    /**
     * get number of ranks (world size)
     * @return world size
    */
    int nRanks() const;

    /**
     * get own rank
     * @return own rank
    */
    int ownRankNo() const;

    /**
     * get rank of left neighbour process.
     * In the case that the current boundary touches the left boundary, its own rank is returned.
     * @return rank of left neighbour process
    */
    int leftNeighbourRankNo() const;

    /**
     * get rank of right neighbour process.
     * In the case that the current boundary touches the right boundary, its own rank is returned.
     * @return rank of right neighbour process
    */
    int rightNeighbourRankNo() const;

    /**
     * get rank of bottom neighbour process.
     * In the case that the current boundary touches the bottom boundary, its own rank is returned.
     * @return rank of bottom neighbour process
    */
    int bottomNeighbourRankNo() const;

    /**
     * get rank of top neighbour process.
     * In the case that the current boundary touches the top boundary, its own rank is returned.
     * @return rank of top neighbour process
    */
    int topNeighbourRankNo() const;

    /**
     * method used to send information to the subdomain left of the current subdomain
     * @param data information to be send to the left
    */
    void sendToLeft(std::vector<double> &data);

    /**
     * method used to send information to the subdomain right of the current subdomain
     * @param data information to be send to the right
    */
    void sendToRight(std::vector<double> &data);

    /**
     * method used to send information to the subdomain below the current subdomain
     * @param data information to be send down
    */    
    void sendToBottom(std::vector<double> &data);

    /**
     * method used to send information to the subdomain above the current subdomain
     * @param data information to be send up
    */    
    void sendToTop(std::vector<double> &data);

    /**
     * method used to receive information from the subdomain left of the current subdomain
     * @param data information to be received from the left
    */    
    void recvFromLeft(std::vector<double> &data, int count);

    /**
     * method used to receive information from the subdomain right of the current subdomain
     * @param data information to be received from the right
    */    
    void recvFromRight(std::vector<double> &data, int count);

    /**
     * method used to receive information from the subdomain below of the current subdomain
     * @param data information to be received from below
    */    
    void recvFromBottom(std::vector<double> &data, int count);

    /**
     * method used to receive information from the subdomain above of the current subdomain
     * @param data information to be received from the top
    */    
    void recvFromTop(std::vector<double> &data, int count);

    /**
     * method used to send information to the subdomain left of the current subdomain
     * @param data information to be send to the left
    */    
    void isendToLeft(std::vector<double> &data, MPI_Request &request);

    /**
     * method used to send information to the subdomain right of the current subdomain
     * @param data information to be send to the right
    */    
    void isendToRight(std::vector<double> &data, MPI_Request &request);

    /**
    * method used to send information to the subdomain below the current subdomain
    * @param data information to be send down
    */    
    void isendToBottom(std::vector<double> &data, MPI_Request &request);
 
    /**
    * method used to send information to the subdomain above the current subdomain
    * @param data information to be send up
    */   
    void isendToTop(std::vector<double> &data, MPI_Request &request);
 
    /**
    * method used to receive information from the subdomain left of the current subdomain
    * @param data information to be received from the left
    */   
    void irecvFromLeft(std::vector<double> &data, int count, MPI_Request &request);

    /**
    * method used to receive information from the subdomain right of the current subdomain
    * @param data information to be received from the right
    */    
    void irecvFromRight(std::vector<double> &data, int count, MPI_Request &request);

    /**
    * method used to receive information from the subdomain below of the current subdomain
    * @param data information to be received from below
    */    
    void irecvFromBottom(std::vector<double> &data, int count, MPI_Request &request);

    /**
    * method used to receive information from the subdomain above of the current subdomain
    * @param data information to be received from the top
    */    
    void irecvFromTop(std::vector<double> &data, int count, MPI_Request &request);

    /**
    * Implementation of call of the MPI-send command
    * @param destinationRank rank of the process, data is send to
    * @param data data to be send
    */    
    void send(int destinationRank, std::vector<double> &data);

    /**
    * Implementation of call of the MPI-receive command
    * @param sourceRank rank of the process, data is received from
    * @param data data to be received
    * @param count size of data to be received
    */    
    void recv(int sourceRank, std::vector<double> &data, int count);

    /**
    * Implementation of call of the MPI-send command
    * @param destinationRank rank of the process, data is send to
    * @param data data to be send
    */    
    void isend(int destinationRank, std::vector<double> &data, MPI_Request &request);

    /**
    * Implementation of call of the MPI-receive command
    * @param sourceRank rank of the process, data is received from
    * @param data data to be received
    * @param count size of data to be received
    */    
    void irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request);

    /**
    * Implementation of call of the MPI-wait command
    * @param request request the program should wait for before continuing
    */    
    void wait(MPI_Request &request);

    /**
    * Implementation of call of MPI-allReduce command,
    * which combines values from all processes and distributes the result back to all processes
    * @param localValue values on subdomains
    * @param op MPI operation to be performed in allReduce
    * @return globalValue combined value
    */    
    double allReduce(double localValue, MPI_Op op);

    /**
    * sum local values over multiple subdomains
    * @param localValue local values on subdomains
    * @return sum over values of multiple subdomains
    */    
    double globalSum(double localValue);

    /**
    *  get maximum of a value over multiple subdomains
    * @param localValue  local values
    * @return global maximum
    */    
    double globalMax(double localValue);

    /**
    *  get minimum of a value over multiple subdomains
    * @param localValue  local values
    * @return global minimum
    */    
    double globalMin(double localValue);

    /**
    * get local number of cells in current subdomain
    * @return number of cells in current subdomain
    */    
    const std::array<int, 2> nCellsLocal() const;

    /**
    * get global number of cells in domain
    * @return number of cells in global domain
    */    
    const std::array<int, 2> nCellsGlobal() const;

    /**
    * get the offset values for counting local nodes in x and y direction.
    * (i_local,j_local) + nodeOffset = (i_global,j_global)
    * used in OutputWriterParaviewParallel
    * @return offset of current subdomain indices compared to the global position
    */    
    std::array<int, 2> nodeOffset() const;

    /**
    * Method for debugging to print out what a process is currently doing
    */    
    void log(const char* message);
private:
    /**
    * this method partitions the global domain in nRanks subdomains  
    * @param nRanks number of processes
    */
    void partitionDomain(int nRanks);
    
    /**
     * this method partitions the global domain in nRanks subdomains  
     * @param nRanks number of processes
    */  
    void partitionDomainEqual(int nRanks);

    /**
    * get discrete subdomain starting index in x-direction
    * @return discrete subdomain starting index in x-direction
    */
    int columnsBegin() const;

    /**
    * get last discrete subdomain index in x-direction
    * @return last discrete subdomain index in x-direction
    */    
    int columnsEnd() const;

    /**
    * get discrete subdomain starting index in y-direction
    * @return discrete subdomain starting index in y-direction
    */   
    int rowsBegin() const;

    /**
    * get last discrete subdomain index in y-direction
    * @return last discrete subdomain index in y-direction
    */    
    int rowsEnd() const;

    /**
    * get the column position of the subdomain
    * @param rank
    * @return column position of the subdomain
    */    
    int computeColumn(int rank) const;

    /**
    * get the row position of the subdomain
    * @param rank unique number of process
    * @return row position of the subdomain
    */    
    int computeRow(int rank) const;

    /**
    * get the rank of a subdomain based on its position in the global grid
    * @param column positions of subdomain in horizontal direction 
    *               on the grid compared to other subdomains
    * @param row positions of subdomain in vertical direction 
    *             on the grid compared to other subdomains
    * @return subdomain rank
    */    
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