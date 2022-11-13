#include "computation/computation.h"
#include "pressure_solver/1_gauss_seidel.h"
#include "pressure_solver/1_sor.h"

//! initialize the computation object
//! parse the settings from file that is given as the only command line argument
void Computation::initialize(string filename)
{
    settings_ = Settings();
    // load settings from file
    settings_.loadFromFile(filename);
    // print settings
    settings_.printSettings();

    // Initialize discretization
    for (int i = 0; i < 2; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    }
    else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    // Initialize solver
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations, settings_.omega);
    }
    else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations);

    } else {
        std::cout << "Solver not found!" << std::endl;
    }

    // Initialize output writers
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
};

//! run the whole simulation until tend
void Computation::runSimulation() {
    int t_iter = 0;
    double time = 0.0;
    // TODO: reach last endTime exactly!
    while (time < settings_.endTime){
        t_iter++;
        applyBoundaryValues();
        computeTimeStepWidth();
        applyPreliminaryBoundaryValues();
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        if (time + dt_ > settings_.endTime) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;
        //output
        cout << "time step " << t_iter << ", t: " << time << "/" << settings_.endTime << ", dt: " << dt_ <<
            ", res. " << pressureSolver_->residualNorm() << ", solver iterations: " << pressureSolver_->iterations() << endl;
        outputWriterText_->writeFile(time);
        outputWriterText_->writePressureFile();
        outputWriterParaview_->writeFile(time);
    }
};

//! set boundary values of u and v to correct values
void Computation::applyBoundaryValues() {
    // set boundary values for u at bottom and top side
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for u at bottom side
        discretization_->u(i, discretization_->uJBegin()) =
                2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        // set boundary values for u at top side
        discretization_->u(i, discretization_->uJEnd() - 1) =
                2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 2);
    }

    // set boundary values for v at bottom and top side
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // set boundary values for v at bottom side
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        // set boundary values for v at top side
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
    }

    // set boundary values for u at left and right side
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for u at left side
        discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        // set boundary values for u at right side
        discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
    }

    // set boundary values for v at left and right side
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for v at left side
        discretization_->v(discretization_->vIBegin(), j) =
                2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        // set boundary values for v at right side
        discretization_->v(discretization_->vIEnd() - 1, j) =
                2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 2, j);
    }
};

//! set boundary values of F and G to correct values
// TODO: verify is this correct
void Computation::applyPreliminaryBoundaryValues(){
    // set boundary values for F at bottom and top side
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for F at bottom side
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
        // set boundary values for F at top side
        discretization_->f(i, discretization_->uJEnd() - 1) = discretization_->u(i, discretization_->uJEnd() - 1);
    }

    // set boundary values for G at bottom and top side
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // set boundary values for v at bottom side
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
        // set boundary values for v at top side
        discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
    }

    // set boundary values for F at left and right side
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for F at left side
        discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
        // set boundary values for F at right side
        discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
    }

    // set boundary values for G at left and right side
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for G at left side
        discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
        // set boundary values for G at right side
        discretization_->g(discretization_->vIEnd() - 1, j) = discretization_->v(discretization_->vIEnd() - 1, j);
    }
};

//! compute the preliminary velocities, F and G
void Computation::computePreliminaryVelocities(){
    // compute F and G in the interior of the domain
    for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            double lap_u = discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j);
            double conv_u = discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j);
            discretization_->f(i,j) = discretization_->u(i,j) + dt_ * (lap_u / settings_.re - conv_u + settings_.g[0]);
        }
    }
    for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            double lap_v = discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j);
            double conv_v = discretization_->computeDv2Dy(i,j) + discretization_->computeDuvDx(i,j);
            discretization_->g(i,j) = discretization_->v(i,j) + dt_ * (lap_v / settings_.re - conv_v + settings_.g[1]);
        }
    }
};

//! solve the Poisson equation for the pressure
void Computation::computePressure() {
    pressureSolver_->solve();
};

//! compute the right hand side of the Poisson equation for the pressure
void Computation::computeRightHandSide() {
    for (int i = discretization_->rhsInteriorIBegin(); i < discretization_->rhsInteriorIEnd(); i++) {
        for (int j = discretization_->rhsInteriorJBegin(); j < discretization_->rhsInteriorJEnd(); j++) {
            double fx = (discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx();
            double gy = (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy();
            discretization_->rhs(i, j) = (fx + gy) / dt_;
        }
    }
};

//! compute the time step width dt from maximum velocities
void Computation::computeTimeStepWidth() {
    double dt_diff = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );

    double dt_conv_u = discretization_->dx() / discretization_->u().absMax();
    double dt_conv_v = discretization_->dy() / discretization_->v().absMax();

    dt_ = settings_.tau * std::min({dt_diff, dt_conv_u, dt_conv_v});
};

//! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
void Computation::computeVelocities() {
    for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }
    for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
};
