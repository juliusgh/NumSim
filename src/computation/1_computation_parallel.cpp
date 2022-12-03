#include "computation/1_computation_parallel.h"
#include "pressure_solver/1_gauss_seidel.h"
#include "pressure_solver/1_sor.h"
#include "pressure_solver/1_red_black.h"
#include "output_writer/output_writer_text_parallel.h"
#include "output_writer/output_writer_paraview_parallel.h"

/**
 * Initialize the ComputationParallel object for a parallel simulation
 * 
 * Parse the settings from the parameter file that is given as the command line argument
 * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
 */
void ComputationParallel::initialize(string filename)
{
    settings_ = Settings();
    // Load settings from file
    settings_.loadFromFile(filename);
    // Print settings
#ifndef NDEBUG
    //settings_.printSettings();
#endif

    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);

    // Initialize discretization
    for (int i = 0; i < 2; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha);
    }
    else {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_);
    }

    // Initialize solver
    pressureSolver_ = std::make_unique<RedBlack>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);

    // Initialize output writers
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
};

/**
 * Run the whole simulation until tend
 */
void ComputationParallel::runSimulation() {
    int t_iter = 0;
    double time = 0.0;
    while (time < settings_.endTime){
        t_iter++;

        /*
        * 1) Apply boundary values (for u, v, F, G)
        */
        partitioning_->log("applyBoundaryValues");
        applyBoundaryValues();
        partitioning_->log("applyPreliminaryBoundaryValues");
        applyPreliminaryBoundaryValues();

        /*
        * 2) Compute the next time step width
        */
        partitioning_->log("computeTimeStepWidth");
        computeTimeStepWidth();
        // endTime should be reached exactly:
        if (time + dt_ > settings_.endTime) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;

        /*
        * 3) Compute preliminary velocities (F, G)
        */
        computePreliminaryVelocities();

        /*
        * 4) Compute the right hand side (rhs)
        */
        computeRightHandSide();

        /*
        * 5) Compute the pressure (p) by solving the Poisson equation
        */
        computePressure();

        /*
        * 6) Update the velocities (u, v)
        */
        computeVelocities(); 

        /*
        * 7) Output debug information and simulation results
        */
#ifndef NDEBUG
        cout << "time step " << t_iter << ", t: " << time << "/" << settings_.endTime << ", dt: " << dt_ <<
             ", res. " << pressureSolver_->residualNorm() << ", solver iterations: " << pressureSolver_->iterations() << endl;
        //outputWriterText_->writePressureFile();
        outputWriterText_->writeFile(time);
#endif
        outputWriterParaview_->writeFile(time);
    }
};

/**
 * Set the boundary values of the velocities (u, v)
 * 
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void ComputationParallel::applyBoundaryValues() {
    int u_columnCount = discretization_->uInteriorJEnd() - discretization_->uInteriorJBegin();
    int u_columnOffset = discretization_->uInteriorJBegin();
    int u_rowCount = discretization_->uInteriorIEnd() - discretization_->uInteriorIBegin();
    int u_rowOffset = discretization_->uInteriorIBegin();
    int v_columnCount = discretization_->vInteriorJEnd() - discretization_->vInteriorJBegin();
    int v_columnOffset = discretization_->vInteriorJBegin();
    int v_rowCount = discretization_->vInteriorIEnd() - discretization_->vInteriorIBegin();
    int v_rowOffset = discretization_->vInteriorIBegin();

    MPI_Request request_v_rightColumn; // send to right - receive from left
    MPI_Request request_v_leftColumn; // send to left - receive from right
    MPI_Request request_v_topRow; // send to top - receive from bottom
    MPI_Request request_v_bottomRow; //senf to bottom - reiceive from top
    std::vector<double> v_rightColumn(v_columnCount, 0);
    std::vector<double> v_leftColumn(v_columnCount, 0);
    std::vector<double> v_topRow(v_rowCount, 0);
    std::vector<double> v_bottomRow(v_rowCount, 0);

    MPI_Request request_u_rightColumn;
    MPI_Request request_u_leftColumn;
    MPI_Request request_u_topRow;
    MPI_Request request_u_bottomRow;
    std::vector<double> u_rightColumn(u_columnCount, 0);
    std::vector<double> u_leftColumn(u_columnCount, 0);
    std::vector<double> u_topRow(u_rowCount, 0);
    std::vector<double> u_bottomRow(u_rowCount, 0);

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        applyBoundaryValuesTop();
    }
    else {
        // v: send second to last row on the top to top neighbour
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            v_topRow.at(i - v_rowOffset) = discretization_->v(i, discretization_->vInteriorJEnd() - 2);
        }
        partitioning_->isendToTop(v_topRow, request_v_topRow);

        // u: send last row on the top to top neighbour
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            u_topRow.at(i - u_rowOffset) = discretization_->u(i, discretization_->uInteriorJEnd() - 1);
        }
        partitioning_->isendToTop(u_topRow, request_u_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecvFromTop(v_topRow, v_rowCount, request_v_topRow);
        partitioning_->irecvFromTop(u_topRow, u_rowCount, request_u_topRow);
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        applyBoundaryValuesBottom();
    }
    else {
        // v: send second row on the bottom to bottom neighbour
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            v_bottomRow.at(i - v_rowOffset) = discretization_->v(i, discretization_->vInteriorJBegin() + 1);
        }
        partitioning_->isendToBottom(v_bottomRow, request_v_bottomRow);

        // u: send first row on the bottom to bottom neighbour
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            u_bottomRow.at(i - u_rowOffset) = discretization_->u(i, discretization_->uInteriorJBegin());
        }
        partitioning_->isendToBottom(u_bottomRow, request_u_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecvFromBottom(v_bottomRow, v_rowCount, request_v_bottomRow);
        partitioning_->irecvFromBottom(u_bottomRow, u_rowCount, request_u_bottomRow);
    }

    // after that set left and right boundary values (high priority)

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        applyBoundaryValuesRight();
    }
    else {
        // v: send last column on the right to right neighbour
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            partitioning_->log("send Right: j - v_columnOffset =");
            std::cout << j - v_columnOffset << std::endl;
            v_rightColumn.at(j - v_columnOffset) = discretization_->v(discretization_->vInteriorIEnd() - 1, j);
        }
        partitioning_->isendToRight(v_rightColumn, request_v_rightColumn);

        // u: send second to last column on the right to right neighbour
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            partitioning_->log("send Right: j - u_columnOffset =");
            std::cout << j - u_columnOffset << std::endl;
            u_rightColumn.at(j - u_columnOffset) = discretization_->u(discretization_->uInteriorIEnd() - 2, j);
        }
        partitioning_->isendToRight(u_rightColumn, request_u_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecvFromRight( v_rightColumn, v_columnCount, request_v_rightColumn);
        partitioning_->irecvFromRight( u_rightColumn, u_columnCount, request_u_rightColumn);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        applyBoundaryValuesLeft();
    }
    else {
        // v: send first column on the left to left neighbour
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            partitioning_->log("send Left: j - v_columnOffset =");
            std::cout << j - v_columnOffset << std::endl;
            v_leftColumn.at(j - v_columnOffset) = discretization_->v(discretization_->vInteriorIBegin(), j);
        }
        partitioning_->isendToLeft(v_leftColumn, request_v_leftColumn);

        // u: send second column on the left to left neighbour
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            partitioning_->log("send Left: j - u_columnOffset =");
            std::cout << j - u_columnOffset << std::endl;
            u_leftColumn.at(j - u_columnOffset) = discretization_->u(discretization_->uInteriorIBegin() + 1, j);
        }
        partitioning_->isendToLeft(u_leftColumn, request_u_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecvFromLeft(v_leftColumn, v_columnCount, request_v_leftColumn);
        partitioning_->irecvFromLeft(u_leftColumn, u_columnCount, request_u_leftColumn);
    }

    // set top and bottom ghost layers first (low priority)

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        // wait until all MPI requests from the top neighbour are finished
        partitioning_->wait(request_v_topRow);
        partitioning_->wait(request_u_topRow);

        // write values from top neighbour to top ghost layer
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            discretization_->v(i, discretization_->vJEnd() - 1) = v_topRow.at(i - v_rowOffset);
        }
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            discretization_->u(i, discretization_->uJEnd() - 1) = u_topRow.at(i - u_rowOffset);
        }
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        // wait until all MPI requests from the bottom neighbour are finished
        partitioning_->wait(request_v_bottomRow);
        partitioning_->wait(request_u_bottomRow);

        // write values from bottom neighbour to bottom ghost layer
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            discretization_->v(i, discretization_->vJBegin()) = v_bottomRow.at(i - v_rowOffset);
        }
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            discretization_->u(i, discretization_->uJBegin()) = u_bottomRow.at(i - u_rowOffset);
        }
    }

    // after that set left and right boundary values (high priority)

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        // wait until all MPI requests from the right neighbour are finished
        partitioning_->wait(request_v_rightColumn);
        partitioning_->wait(request_u_rightColumn);

        // write values from right neighbour to right ghost layer
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            partitioning_->log("recv Right: j - v_columnOffset =");
            std::cout << j - v_columnOffset << std::endl;
            discretization_->v(discretization_->vIEnd() - 1, j) = v_rightColumn.at(j - v_columnOffset);
        }
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            partitioning_->log("recv Right: j - u_columnOffset =");
            std::cout << j - u_columnOffset << std::endl;
            discretization_->u(discretization_->uIEnd() - 1, j) = u_rightColumn.at(j - u_columnOffset);
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        // wait until all MPI requests from the left neighbour are finished
        partitioning_->wait(request_v_leftColumn);
        partitioning_->wait(request_u_leftColumn);

        // write values from left neighbour to left ghost layer
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            partitioning_->log("recv Left: j - v_columnOffset =");
            std::cout << j - v_columnOffset << std::endl;
            discretization_->v(discretization_->vIBegin(), j) = v_leftColumn.at(j - v_columnOffset);
        }
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            partitioning_->log("recv Left: j - u_columnOffset =");
            std::cout << j - u_columnOffset << std::endl;
            discretization_->u(discretization_->uIBegin(), j) = u_leftColumn.at(j - u_columnOffset);
        }
    }
}

/**
 * Compute the time step width dt based on the maximum velocities
 */
void ComputationParallel::computeTimeStepWidth() {
    // Compute maximal time step width regarding the diffusion
    double dt_diff = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );

    // Compute maximal time step width regarding the convection u
    double u_absMax_local = discretization_->u().absMax();
    partitioning_->log("u_absMax_local:");
    std::cout << u_absMax_local << std::endl;
    discretization_->u().print();
    double u_absMax = partitioning_->globalMax(u_absMax_local);
    partitioning_->log("u_absMax:");
    std::cout << u_absMax << std::endl;
    double dt_conv_u = std::numeric_limits<double>::max();
    if (u_absMax > 0.0)
        dt_conv_u = discretization_->dx() / u_absMax;


    // Compute maximal time step width regarding the convection v
    double v_absMax_local = discretization_->v().absMax();
    partitioning_->log("v_absMax_local:");
    std::cout << v_absMax_local << std::endl;
    double v_absMax = partitioning_->globalMax(v_absMax_local);
    partitioning_->log("v_absMax:");
    std::cout << v_absMax << std::endl;
    double dt_conv_v = std::numeric_limits<double>::max();
    if (v_absMax > 0.0)
        dt_conv_v = discretization_->dy() / v_absMax;

    // Set the appropriate time step width by using a security factor tau
    dt_ = settings_.tau * std::min({dt_diff, dt_conv_u, dt_conv_v});
}
