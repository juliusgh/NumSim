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
    settings_.printSettings();
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
        applyBoundaryValues();
        applyPreliminaryBoundaryValues();

        /*
        * 2) Compute the next time step width
        */
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
    int columnCount = discretization_->pInteriorJEnd() - discretization_->pInteriorJBegin();
    int columnOffset = discretization_->pInteriorJBegin();
    int rowCount = discretization_->pInteriorIEnd() - discretization_->pInteriorIBegin();
    int rowOffset = discretization_->pInteriorIBegin();

    MPI_Request request_p_rightColumn;
    MPI_Request request_p_leftColumn;
    MPI_Request request_p_topRow;
    MPI_Request request_p_bottomRow;
    std::vector<double> p_rightColumn;
    std::vector<double> p_leftColumn;
    std::vector<double> p_topRow;
    std::vector<double> p_bottomRow;

    MPI_Request request_v_rightColumn;
    MPI_Request request_v_leftColumn;
    MPI_Request request_v_topRow;
    MPI_Request request_v_bottomRow;
    std::vector<double> v_rightColumn;
    std::vector<double> v_leftColumn;
    std::vector<double> v_topRow;
    std::vector<double> v_bottomRow;

    MPI_Request request_u_rightColumn;
    MPI_Request request_u_leftColumn;
    MPI_Request request_u_topRow;
    MPI_Request request_u_bottomRow;
    std::vector<double> u_rightColumn;
    std::vector<double> u_leftColumn;
    std::vector<double> u_topRow;
    std::vector<double> u_bottomRow;


    if (partitioning_->ownPartitionContainsRightBoundary()) {
        applyBoundaryValuesRight();
    }
    else {
        // p: send last column on the right to right neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            p_rightColumn.push_back(discretization_->p(discretization_->pInteriorIEnd(), j));
        }
        partitioning_->isend(partitioning_->rightNeighbourRankNo(), p_rightColumn, request_p_rightColumn);

        // v: send last column on the right to right neighbour
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            v_rightColumn.push_back(discretization_->v(discretization_->vInteriorIEnd(), j));
        }
        partitioning_->isend(partitioning_->rightNeighbourRankNo(), v_rightColumn, request_v_rightColumn);
        // u: send second to last column on the right to right neighbour
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            u_rightColumn.push_back(discretization_->u(discretization_->uInteriorIEnd()-1, j));
        }
        partitioning_->isend(partitioning_->rightNeighbourRankNo(), u_rightColumn, request_u_rightColumn);
        

        // receive ghost layer column on the right from right neighbour
        partitioning_->irecv(partitioning_->rightNeighbourRankNo(), p_rightColumn, columnCount, request_p_rightColumn);
        partitioning_->irecv(partitioning_->rightNeighbourRankNo(), v_rightColumn, columnCount, request_v_rightColumn);
        partitioning_->irecv(partitioning_->rightNeighbourRankNo(), u_rightColumn, columnCount, request_u_rightColumn);
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        applyBoundaryValuesLeft();
    }
    else {
        // p: send first column on the left to left neighbour
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            p_leftColumn.push_back(discretization_->p(discretization_->pInteriorIBegin(), j));
        }
        partitioning_->isend(partitioning_->leftNeighbourRankNo(), p_leftColumn, request_p_leftColumn);
        // v: send first column on the left to left neighbour
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vpInteriorJEnd(); j++) {
            v_leftColumn.push_back(discretization_->v(discretization_->vInteriorIBegin(), j));
        }
        partitioning_->isend(partitioning_->leftNeighbourRankNo(), v_leftColumn, request_v_leftColumn);
        // u: send second column on the left to left neighbour
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            p_leftColumn.push_back(discretization_->u(discretization_->uInteriorIBegin()+1, j));
        }
        partitioning_->isend(partitioning_->leftNeighbourRankNo(), u_leftColumn, request_u_leftColumn);

        // receive ghost layer column on the left from left neighbour
        partitioning_->irecv(partitioning_->leftNeighbourRankNo(), p_leftColumn, columnCount, request_p_leftColumn);
        partitioning_->irecv(partitioning_->leftNeighbourRankNo(), v_leftColumn, columnCount, request_v_leftColumn);
        partitioning_->irecv(partitioning_->leftNeighbourRankNo(), u_leftColumn, columnCount, request_u_leftColumn);

    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        applyBoundaryValuesTop();
    }
    else {
        // p: send last row on the top to top neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            p_topRow.push_back(discretization_->p(i, discretization_->pInteriorJEnd()));
        }
        partitioning_->isend(partitioning_->topNeighbourRankNo(), p_topRow, request_p_topRow);
        // v: send second to last row on the top to top neighbour
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            v_topRow.push_back(discretization_->v(i, discretization_->vInteriorJEnd()-1));
        }
        partitioning_->isend(partitioning_->topNeighbourRankNo(), v_topRow, request_v_topRow);
        // u: send last row on the top to top neighbour
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            u_topRow.push_back(discretization_->u(i, discretization_->uInteriorJEnd()));
        }
        partitioning_->isend(partitioning_->topNeighbourRankNo(), u_topRow, request_u_topRow);


        // receive ghost layer row on the top from top neighbour
        partitioning_->irecv(partitioning_->topNeighbourRankNo(), p_topRow, rowCount, request_p_topRow);
        partitioning_->irecv(partitioning_->topNeighbourRankNo(), v_topRow, rowCount, request_v_topRow);
        partitioning_->irecv(partitioning_->topNeighbourRankNo(), u_topRow, rowCount, request_u_topRow);
    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        applyBoundaryValuesBottom();
    }
    else {
        // p: send first row on the bottom to bottom neighbour
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            p_bottomRow.push_back(discretization_->p(i, discretization_->pInteriorJBegin()));
        }
        partitioning_->isend(partitioning_->bottomNeighbourRankNo(), p_bottomRow, request_p_bottomRow);
        // v: send second row on the bottom to bottom neighbour
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            v_bottomRow.push_back(discretization_->v(i, discretization_->vInteriorJBegin()+1));
        }
        partitioning_->isend(partitioning_->bottomNeighbourRankNo(), v_bottomRow, request_v_bottomRow);
        // u: send first row on the bottom to bottom neighbour
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            u_bottomRow.push_back(discretization_->u(i, discretization_->uInteriorJBegin()));
        }
        partitioning_->isend(partitioning_->bottomNeighbourRankNo(), u_bottomRow, request_u_bottomRow);
        
        

        // receive ghost layer row on the bottom from bottom neighbour
        partitioning_->irecv(partitioning_->bottomNeighbourRankNo(), p_bottomRow, rowCount, request_p_bottomRow);
        partitioning_->irecv(partitioning_->bottomNeighbourRankNo(), v_bottomRow, rowCount, request_v_bottomRow);
        partitioning_->irecv(partitioning_->bottomNeighbourRankNo(), u_bottomRow, rowCount, request_u_bottomRow);
    }


    partitioning_->wait(request_p_rightColumn);
    partitioning_->wait(request_v_rightColumn);
    partitioning_->wait(request_u_rightColumn);
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            discretization_->p(discretization_->pIEnd(), j) = p_rightColumn.at(j - columnOffset);
        }
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            discretization_->v(discretization_->vIEnd(), j) = v_rightColumn.at(j - columnOffset);
        }
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            discretization_->u(discretization_->uIEnd(), j) = u_rightColumn.at(j - columnOffset);
        }
    }
    partitioning_->wait(request_p_leftColumn);
    partitioning_->wait(request_v_leftColumn);
    partitioning_->wait(request_u_leftColumn);
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = discretization_->pInteriorJBegin(); j < discretization_->pInteriorJEnd(); j++) {
            discretization_->p(discretization_->pIBegin(), j) = p_leftColumn.at(j - columnOffset);
        }
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            discretization_->v(discretization_->vIBegin(), j) = v_leftColumn.at(j - columnOffset);
        }
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            discretization_->u(discretization_->uIBegin(), j) = u_leftColumn.at(j - columnOffset);
        }
    }
    partitioning_->wait(request_p_topRow);
    partitioning_->wait(request_v_topRow);
    partitioning_->wait(request_u_topRow);
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = p_topRow.at(i - rowOffset);
        }
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            discretization_->v(i, discretization_->vJEnd()) = v_topRow.at(i - rowOffset);
        }
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            discretization_->u(i, discretization_->uJEnd()) = u_topRow.at(i - rowOffset);
        }
    }
    partitioning_->wait(request_p_bottomRow);
    partitioning_->wait(request_v_bottomRow);
    partitioning_->wait(request_u_bottomRow);
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = discretization_->pInteriorIBegin(); i < discretization_->pInteriorIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin()) = p_bottomRow.at(i - rowOffset);
        }
        for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
            discretization_->v(i, discretization_->vJBegin()) = v_bottomRow.at(i - rowOffset);
        }
        for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
            discretization_->u(i, discretization_->uJBegin()) = u_bottomRow.at(i - rowOffset);
        }
    }
};


/**
 * Compute the time step width dt based on the maximum velocities
 */
void ComputationParallel::computeTimeStepWidth() {
    // Compute maximal time step width regarding the diffusion
    double dt_diff = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );

    // Compute maximal time step width regarding the convection u
    double u_absMax_local = discretization_->u().absMax();
    double u_absMax = partitioning_->globalMax(u_absMax_local);
    double dt_conv_u = discretization_->dx() / u_absMax;


    // Compute maximal time step width regarding the convection v
    double v_absMax_local = discretization_->v().absMax();
    double v_absMax = partitioning_->globalMax(v_absMax_local);
    double dt_conv_v = discretization_->dy() / v_absMax;

    // Set the appropriate time step width by using a security factor tau
    dt_ = settings_.tau * std::min({dt_diff, dt_conv_u, dt_conv_v});
};



void ComputationParallel::applyBoundaryValuesTop(){
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for u at top side
        discretization_->u(i, discretization_->uJEnd() - 1) =
                2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uInteriorJEnd() - 1);
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        ;
        // set boundary values for v at top side
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
    }
};

void ComputationParallel::applyBoundaryValuesBottom(){
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for u at bottom side
        discretization_->u(i, discretization_->uJBegin()) =
                2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uInteriorJBegin());
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // set boundary values for v at bottom side
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
    }  
};

void ComputationParallel::applyBoundaryValuesLeft(){
    // set boundary values for u at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for u at left side
        discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
    }

    // set boundary values for v at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for v at left side
        discretization_->v(discretization_->vIBegin(), j) =
                2.0 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
    }
    
};

void ComputationParallel::applyBoundaryValuesRight(){
    // set boundary values for u at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for u at right side
        discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
    }

    // set boundary values for v at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for v at right side
        discretization_->v(discretization_->vIEnd() - 1, j) =
                2.0 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vInteriorIEnd() - 1, j);
    }
    
};
