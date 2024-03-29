#include "computation/1_computation_parallel.h"
#include "pressure_solver/1_gauss_seidel.h"
#include "pressure_solver/1_sor.h"
#include "pressure_solver/2_red_black.h"
#include "pressure_solver/2_conjugate_gradient_parallel.h"
#include "output_writer/output_writer_text_parallel.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "pressure_solver/2_red_black_sor.h"
#include <ctime>

/**
 * Initialize the ComputationParallel object for a parallel simulation
 * 
 * Parse the settings from the parameter file that is given as the command line argument
 * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
 */
void ComputationParallel::initialize(string filename) {
    std::cout << "I am parallel" << std::endl;
    settings_ = Settings();
    // Load settings from file
    settings_.loadFromFile(filename);

#ifndef NDEBUG
    // Print settings
    settings_.printSettings();
#endif

    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);

    // Initialize discretization
    for (int i = 0; i < 2; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha, settings_.gamma, settings_);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_, settings_);
    }


    // Initialize solver
#ifndef NDEBUG
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<RedBlackSOR>(discretization_, settings_.epsilon,
                                                        settings_.maximumNumberOfIterations, settings_.omega,
                                                        partitioning_);
    } else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<RedBlack>(discretization_, settings_.epsilon,
                                                     settings_.maximumNumberOfIterations, partitioning_);
    } else if (settings_.pressureSolver == "CG") {
        pressureSolver_ = std::make_unique<ConjugateGradientParallel>(discretization_, settings_.epsilon,
                                                                      settings_.maximumNumberOfIterations, partitioning_);
    } else {
        std::cout << "Solver not found!" << std::endl;
    }
#else
    pressureSolver_ = std::make_unique<ConjugateGradientParallel>(discretization_, settings_.epsilon,
                                                     settings_.maximumNumberOfIterations, partitioning_);
#endif



    // Initialize output writers
#ifndef NDEBUG
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);
#endif
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
}

/**
 * Run the whole simulation until tend
 */
void ComputationParallel::runSimulation() {
    int t_iter = 0;
    double time = 0.0;
#ifndef NDEBUG
    const clock_t sim_start_time = clock();
    float sim_duration;
#else
    double time_steps_print = 1;
    double time_last_printed = -time_steps_print;
#endif
    while (time < settings_.endTime) {
        t_iter++;

        /*
        * 1) Apply boundary values (for u, v, F, G)
        */
        //partitioning_->log("applyBoundaryValues");
        applyBoundaryValues();
        //partitioning_->log("applyPreliminaryBoundaryValues");
        applyPreliminaryBoundaryValues();

        /*
        * 2) Compute the next time step width
        */
        //partitioning_->log("computeTimeStepWidth");
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
        if (partitioning_->ownRankNo() == 0) {
            cout << "time step " << t_iter << ", t: " << time << "/" << settings_.endTime << ", dt: " << dt_ <<
                 ", res. " << pressureSolver_->residualNorm() << ", solver iterations: "
                 << pressureSolver_->iterations() << endl;
        }
        outputWriterText_->writeFile(time);
        outputWriterParaview_->writeFile(time);
        sim_duration = float(clock() - sim_start_time) / CLOCKS_PER_SEC;
        if (sim_duration >= 600) {
            std::cout << "TERMINATED (cause: runTimeError)" << std::endl;
            break;
        }
#else

        if (time - time_last_printed >= time_steps_print){
            outputWriterParaview_->writeFile(time);
            time_last_printed += time_steps_print;
        }

#endif

    }
#ifndef NDEBUG
    std::cout << "Total computation time: " << sim_duration << "s" << std::endl;
#endif
}


/**
 * Set the boundary values of the velocities (u, v)
 * 
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void ComputationParallel::applyBoundaryValues() {
    const int vIBegin = discretization_->vIBegin();
    const int vJBegin = discretization_->vJBegin();
    const int uIBegin = discretization_->uIBegin();
    const int uJBegin = discretization_->uJBegin();
    const int vIEnd = discretization_->vIEnd();
    const int vJEnd = discretization_->vJEnd();
    const int uIEnd = discretization_->uIEnd();
    const int uJEnd = discretization_->uJEnd();
    const int vInteriorIBegin = discretization_->vInteriorIBegin();
    const int vInteriorJBegin = discretization_->vInteriorJBegin();
    const int uInteriorIBegin = discretization_->uInteriorIBegin();
    const int uInteriorJBegin = discretization_->uInteriorJBegin();
    const int vInteriorIEnd = discretization_->vInteriorIEnd();
    const int vInteriorJEnd = discretization_->vInteriorJEnd();
    const int uInteriorIEnd = discretization_->uInteriorIEnd();
    const int uInteriorJEnd = discretization_->uInteriorJEnd();
    const int u_columnCount = uInteriorJEnd - uInteriorJBegin;
    int u_columnOffset = uInteriorJBegin;
    int u_rowCount = uInteriorIEnd - uInteriorIBegin;
    int u_rowOffset = uInteriorIBegin;
    int v_columnCount = vInteriorJEnd - vInteriorJBegin;
    int v_columnOffset = vInteriorJBegin;
    int v_rowCount = vInteriorIEnd - vInteriorIBegin;
    int v_rowOffset = vInteriorIBegin;

    MPI_Request request_v_rightColumn; // send to right - receive from left
    MPI_Request request_v_leftColumn; // send to left - receive from right
    MPI_Request request_v_topRow; // send to top - receive from bottom
    MPI_Request request_v_bottomRow; //send to bottom - reiceive from top
    std::vector<double> v_rightColumn(v_columnCount, 0);
    std::vector<double> v_leftColumn(v_columnCount, 0);
    std::vector<double> v_topRow(v_rowCount, 0);
    std::vector<double> v_bottomRow(v_rowCount, 0);

    MPI_Request request_u_rightColumn; // send to right - receive from left
    MPI_Request request_u_leftColumn; // send to left - receive from right
    MPI_Request request_u_topRow; // send to top - receive from bottom
    MPI_Request request_u_bottomRow; //send to bottom - reiceive from top
    std::vector<double> u_rightColumn(u_columnCount, 0);
    std::vector<double> u_leftColumn(u_columnCount, 0);
    std::vector<double> u_topRow(u_rowCount, 0);
    std::vector<double> u_bottomRow(u_rowCount, 0);

    /*
    * velocities u and v communication: send to and receive from subdomain directly above 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches upper domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        applyBoundaryValuesTop();
    } else {
        // v: send second to last row on the top to top neighbour
        for (int i = vInteriorIBegin; i < vInteriorIEnd; i++) {
            v_topRow.at(i - v_rowOffset) = discretization_->v(i, vInteriorJEnd - 2);
        }
        partitioning_->isendToTop(v_topRow, request_v_topRow);

        // u: send last row on the top to top neighbour
        for (int i = uInteriorIBegin; i < uInteriorIEnd; i++) {
            u_topRow.at(i - u_rowOffset) = discretization_->u(i, uInteriorJEnd - 1);
        }
        partitioning_->isendToTop(u_topRow, request_u_topRow);

        // receive ghost layer row on the top from top neighbour
        partitioning_->irecvFromTop(v_topRow, v_rowCount, request_v_topRow);
        partitioning_->irecvFromTop(u_topRow, u_rowCount, request_u_topRow);
    }

    /*
    * velocities u and v communication: send to and receive from subdomain directly below 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches lower domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        applyBoundaryValuesBottom();
    } else {
        // v: send second row on the bottom to bottom neighbour
        for (int i = vInteriorIBegin; i < vInteriorIEnd; i++) {
            v_bottomRow.at(i - v_rowOffset) = discretization_->v(i, vInteriorJBegin + 1);
        }
        partitioning_->isendToBottom(v_bottomRow, request_v_bottomRow);

        // u: send first row on the bottom to bottom neighbour
        for (int i = uInteriorIBegin; i < uInteriorIEnd; i++) {
            u_bottomRow.at(i - u_rowOffset) = discretization_->u(i, uInteriorJBegin);
        }
        partitioning_->isendToBottom(u_bottomRow, request_u_bottomRow);

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecvFromBottom(v_bottomRow, v_rowCount, request_v_bottomRow);
        partitioning_->irecvFromBottom(u_bottomRow, u_rowCount, request_u_bottomRow);
    }

    /*
    * after that set left and right boundary values (high priority)
    */
    /*
    * velocities u and v communication: send to and receive from subdomain directly right of 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches right domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsRightBoundary()) {
        applyBoundaryValuesRight();
    } else {
        // v: send last column on the right to right neighbour
        for (int j = vInteriorJBegin; j < vInteriorJEnd; j++) {
            v_rightColumn.at(j - v_columnOffset) = discretization_->v(vInteriorIEnd - 1, j);
        }
        partitioning_->isendToRight(v_rightColumn, request_v_rightColumn);

        // u: send second to last column on the right to right neighbour
        for (int j = uInteriorJBegin; j < uInteriorJEnd; j++) {

            u_rightColumn.at(j - u_columnOffset) = discretization_->u(uInteriorIEnd - 2, j);
        }
        partitioning_->isendToRight(u_rightColumn, request_u_rightColumn);

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecvFromRight(v_rightColumn, v_columnCount, request_v_rightColumn);
        partitioning_->irecvFromRight(u_rightColumn, u_columnCount, request_u_rightColumn);
    }
    /*
    * velocities u and v communication: send to and receive from subdomain directly left of 
    * current subdomain if current subdomain is not on upper boundary. 
    * If subdomain touches left domain boundary, boundary conditions are applied.
    */
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        applyBoundaryValuesLeft();
    } else {
        // v: send first column on the left to left neighbour
        for (int j = vInteriorJBegin; j < vInteriorJEnd; j++) {

            v_leftColumn.at(j - v_columnOffset) = discretization_->v(vInteriorIBegin, j);
        }
        partitioning_->isendToLeft(v_leftColumn, request_v_leftColumn);

        // u: send second column on the left to left neighbour
        for (int j = uInteriorJBegin; j < uInteriorJEnd; j++) {
            //partitioning_->log("send Left: j - u_columnOffset =");
            //std::cout << j - u_columnOffset << std::endl;
            u_leftColumn.at(j - u_columnOffset) = discretization_->u(uInteriorIBegin + 1, j);
        }
        partitioning_->isendToLeft(u_leftColumn, request_u_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecvFromLeft(v_leftColumn, v_columnCount, request_v_leftColumn);
        partitioning_->irecvFromLeft(u_leftColumn, u_columnCount, request_u_leftColumn);
    }

    /*
    set top and bottom ghost layers first (low priority)
    */
    /*
    * Set subdomain ghost layer at top for velocities
    */
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        // wait until all MPI requests from the top neighbour are finished
        partitioning_->wait(request_v_topRow);
        partitioning_->wait(request_u_topRow);

        // write values from top neighbour to top ghost layer
        for (int i = vInteriorIBegin; i < vInteriorIEnd; i++) {
            discretization_->v(i, vJEnd - 1) = v_topRow.at(i - v_rowOffset);
        }
        for (int i = uInteriorIBegin; i < uInteriorIEnd; i++) {
            discretization_->u(i, uJEnd - 1) = u_topRow.at(i - u_rowOffset);
        }
    }

    /* 
    * Set subdomain ghost layer at bottom for velocities
    */
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        // wait until all MPI requests from the bottom neighbour are finished
        partitioning_->wait(request_v_bottomRow);
        partitioning_->wait(request_u_bottomRow);

        // write values from bottom neighbour to bottom ghost layer
        for (int i = vInteriorIBegin; i < vInteriorIEnd; i++) {
            discretization_->v(i, vJBegin) = v_bottomRow.at(i - v_rowOffset);
        }
        for (int i = uInteriorIBegin; i < uInteriorIEnd; i++) {
            discretization_->u(i, uJBegin) = u_bottomRow.at(i - u_rowOffset);
        }
    }

    /*
    * after that set left and right boundary values (high priority)
    */
    /* 
    * Set subdomain ghost layer at right side for velocities
    */
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        // wait until all MPI requests from the right neighbour are finished
        partitioning_->wait(request_v_rightColumn);
        partitioning_->wait(request_u_rightColumn);

        // write values from right neighbour to right ghost layer
        for (int j = vInteriorJBegin; j < vInteriorJEnd; j++) {

            discretization_->v(vIEnd - 1, j) = v_rightColumn.at(j - v_columnOffset);
        }
        for (int j = uInteriorJBegin; j < uInteriorJEnd; j++) {

            discretization_->u(uIEnd - 1, j) = u_rightColumn.at(j - u_columnOffset);
        }
    }
    /* 
    * Set subdomain ghost layer at left side for velocities
    */
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        // wait until all MPI requests from the left neighbour are finished
        partitioning_->wait(request_v_leftColumn);
        partitioning_->wait(request_u_leftColumn);

        // write values from left neighbour to left ghost layer
        for (int j = vInteriorJBegin; j < vInteriorJEnd; j++) {

            discretization_->v(vIBegin, j) = v_leftColumn.at(j - v_columnOffset);
        }
        for (int j = uInteriorJBegin; j < uInteriorJEnd; j++) {

            discretization_->u(uIBegin, j) = u_leftColumn.at(j - u_columnOffset);
        }
    }
}

/**
 * Compute the time step width dt based on the maximum velocities 
 * over all subdomains 
 */
void ComputationParallel::computeTimeStepWidth() {
    const double dx = discretization_->dx();
    const double dy = discretization_->dy();

    // Compute maximal time step width regarding the diffusion
    double dt_diff = settings_.re / 2 / (1 / (dx * dx) + 1 / (dy * dy));

    // Compute maximal time step width regarding the convection u
    double u_absMax_local = discretization_->u().absMax();

    double u_absMax = partitioning_->globalMax(u_absMax_local);

    double dt_conv_u = std::numeric_limits<double>::max();
    if (u_absMax > 0.0)
        dt_conv_u = dx / u_absMax;


    // Compute maximal time step width regarding the convection v
    double v_absMax_local = discretization_->v().absMax();

    double v_absMax = partitioning_->globalMax(v_absMax_local);

    double dt_conv_v = std::numeric_limits<double>::max();
    if (v_absMax > 0.0)
        dt_conv_v = dy / v_absMax;

    // Set the appropriate time step width by using a security factor tau
    dt_ = settings_.tau * std::min({dt_diff, dt_conv_u, dt_conv_v});
}
