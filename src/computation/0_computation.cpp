#include "computation/0_computation.h"
#include "pressure_solver/1_gauss_seidel.h"
#include "pressure_solver/1_sor.h"
#include "pressure_solver/2_conjugate_gradient.h"

/**
 * Initialize the computation object for a sequential simulation
 * 
 * Parse the settings from the parameter file that is given as the command line argument
 * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
 */
void Computation::initialize(string filename) {
    std::cout << "I am not parallel" << std::endl;
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
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha, settings_.gamma, std::shared_ptr<Settings>(&settings_));
    } else {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_, std::shared_ptr<Settings>(&settings_));
    }

    // Initialize solver
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations, settings_.omega);
    } else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
                                                        settings_.maximumNumberOfIterations);

    } else if (settings_.pressureSolver == "CG") {
        pressureSolver_ = std::make_unique<ConjugateGradient>(discretization_, settings_.epsilon,
                                                      settings_.maximumNumberOfIterations, partitioning_);
    } else {
        std::cout << "Solver not found!" << std::endl;
    }

    // Initialize output writers
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
};

/**
 * Run the whole simulation until tend
 */
void Computation::runSimulation() {
#ifndef NDEBUG
    std::cout << "Running simulation ..." << std::endl;
#endif
    int t_iter = 0;
    double time = 0.0;
    setInitialValues();
#ifndef NDEBUG
    std::cout << "initialized" << std::endl;
#endif
    while (time < settings_.endTime) {
        t_iter++;

        /*
        * 1) Apply boundary values (for u, v, F, G)
        */
        applyBoundaryValues();
        //applyPreliminaryBoundaryValues();

#ifndef NDEBUG
        //std::cout << "Preliminary Boundary values applied" << std::endl;
#endif

        /*
        * 2) Compute the next time step width
        */
        computeTimeStepWidth();
        // endTime should be reached exactly:
        if (time + dt_ > settings_.endTime - dt_ / 100) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;

        /*
         * 2.5) Compute Temperature t
         */
        computeTemperature();

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
        updateLastVelocities();
        computeVelocities();
        discretization_->setObstacleValues();

        /*
        * 7) Output debug information and simulation results
        */
#ifndef NDEBUG
        cout << "time step " << t_iter << ", t: " << time << "/" << settings_.endTime << ", dt: " << dt_ <<
             ", res. " << pressureSolver_->residualNorm() << ", solver iterations: " << pressureSolver_->iterations()
             << endl;
        outputWriterText_->writePressureFile();
        outputWriterText_->writeFile(time);
#endif
        outputWriterParaview_->writeFile(time);
    }
};

/**
 * Set the initial values of the temperature t
 */
void Computation::applyBoundaryValues() {
    discretization_->applyBoundaryVelocities();
    discretization_->applyBoundaryTemperature();
    setExternalHeat();
};

/**
 * Set the boundary values of the velocities (u, v)
 * 
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void Computation::setInitialValues() {
    for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++) {
        for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++) {
            discretization_->t(i, j) = settings_.initialTemp;
        }
    }
};

/**
 * Set the boundary values of the velocities (u, v)
 *
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void Computation::setExternalHeat() {
    /*int mI = (discretization_->tInteriorIBegin() + discretization_->tInteriorIEnd()) / 2;
    int mJ = (discretization_->tInteriorJBegin() + discretization_->tInteriorJEnd()) / 2;

    for (int i = mI - 2; i < mI + 2; i++) {
        for (int j = mJ - 2; j < mJ + 2; j++) {
            discretization_->q(i, j) = 0.0;
        }
    }
    for (int i = discretization_->tInteriorIBegin(); i < discretization_->tInteriorIEnd(); i++) {
        discretization_->q(i, discretization_->tInteriorJBegin()) = 100;
    }*/
};

void Computation::applyBoundaryValuesTop() {
    // set boundary values for u at top side
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        if (settings_.outflowTop) {
            discretization_->u(i, discretization_->uJEnd() - 1) =
                    discretization_->u(i, discretization_->uInteriorJEnd() - 1);
        } else {
            discretization_->u(i, discretization_->uJEnd() - 1) =
                    2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uInteriorJEnd() - 1);
        }
    }

    // set boundary values for v at top side
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        if (settings_.outflowTop) {
            discretization_->v(i, discretization_->vJEnd() - 1) =
                    discretization_->v(i, discretization_->vInteriorJEnd() - 1);
        } else {
            discretization_->v(i, discretization_->vJEnd() - 1) =
                    settings_.dirichletBcTop[1];
        }
    }

    // set boundary values for t at top side
    for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++) {
        if (settings_.setFixedTempTop) {
            discretization_->t(i, discretization_->tJEnd() - 1) =
                    2.0 * settings_.tempBcTop
                    - discretization_->t(i, discretization_->tInteriorJEnd() - 1);
        } else {
            discretization_->t(i, discretization_->tJEnd() - 1) =
                    discretization_->t(i, discretization_->tInteriorJEnd() - 1)
                    - discretization_->dy() * settings_.tempBcTop;
        }
    }
};

void Computation::applyBoundaryValuesBottom() {
    // set boundary values for u at bottom side
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        if (settings_.outflowBottom) {
            discretization_->u(i, discretization_->uJBegin()) =
                    discretization_->u(i, discretization_->uInteriorJBegin());
        } else {
            discretization_->u(i, discretization_->uJBegin()) =
                    2.0 * settings_.dirichletBcBottom[0]
                    - discretization_->u(i, discretization_->uInteriorJBegin());
        }
    }

    // set boundary values for v at bottom side
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        if (settings_.outflowBottom) {
            discretization_->v(i, discretization_->vJBegin()) =
                    discretization_->v(i, discretization_->vInteriorJBegin());
        } else {
            discretization_->v(i, discretization_->vJBegin()) =
                    settings_.dirichletBcBottom[1];
        }
    }

    // set boundary values for t at bottom side
    for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++) {
        if (settings_.setFixedTempBottom) {
            //std::cout << "tempBcBottom = " << settings_.tempBcBottom << std::endl;
            discretization_->t(i, discretization_->tJBegin()) =
                    2.0 * settings_.tempBcBottom - discretization_->t(i, discretization_->tInteriorJBegin());
            //std::cout << "t(i, discretization_->tJBegin()) = " << discretization_->t(i, discretization_->tJBegin()) << std::endl;
        } else {
            discretization_->t(i, discretization_->tJBegin()) =
                    discretization_->t(i, discretization_->tInteriorJBegin())
                    - discretization_->dy() * settings_.tempBcBottom;
        }
    }
};

void Computation::applyBoundaryValuesLeft() {
    // set boundary values for u at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for u at left side
        if (settings_.outflowLeft) {
            discretization_->u(discretization_->uIBegin(), j) =
                    discretization_->u(discretization_->uInteriorIBegin(), j);
        } else {
            discretization_->u(discretization_->uIBegin(), j) =
                    settings_.dirichletBcLeft[0];
        }
    }

    // set boundary values for v at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for v at left side
        if (settings_.outflowLeft) {
            discretization_->v(discretization_->vIBegin(), j) =
                    discretization_->v(discretization_->vInteriorIBegin(), j);
        } else {
            discretization_->v(discretization_->vIBegin(), j) =
                    2.0 * settings_.dirichletBcLeft[1]
                    - discretization_->v(discretization_->vIBegin() + 1, j);
        }
    }

    // set boundary values for t at left side
    for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++) {
        if (settings_.setFixedTempLeft) {
            discretization_->t(discretization_->tIBegin(), j) =
                    2.0 * settings_.tempBcLeft
                    - discretization_->t(discretization_->tInteriorIBegin(), j);
        } else {
            discretization_->t(discretization_->tIBegin(), j) =
                    discretization_->t(discretization_->tInteriorIBegin(), j)
                    - discretization_->dx() * settings_.tempBcLeft;
        }
    }
};

void Computation::applyBoundaryValuesRight() {
    // set boundary values for u at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for u at right side
        if (settings_.outflowRight) {
            discretization_->u(discretization_->uIEnd() - 1, j) =
                    discretization_->u(discretization_->uInteriorIEnd() - 1, j);
        } else {
            discretization_->u(discretization_->uIEnd() - 1, j) =
                    settings_.dirichletBcRight[0];
        }
    }

    // set boundary values for v at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for v at right side
        if (settings_.outflowRight) {
            discretization_->v(discretization_->vIBegin(), j) =
                    discretization_->v(discretization_->vInteriorIBegin(), j);
        } else {
            discretization_->v(discretization_->vIEnd() - 1, j) =
                    2.0 * settings_.dirichletBcRight[1]
                    - discretization_->v(discretization_->vInteriorIEnd() - 1, j);
        }
    }

    // set boundary values for t at right side
    for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++) {
        if (settings_.setFixedTempRight) {
            discretization_->t(discretization_->tIEnd() - 1, j) =
                    2.0 * settings_.tempBcRight
                    - discretization_->t(discretization_->tInteriorIEnd() - 1, j);
        } else {
            discretization_->t(discretization_->tIEnd() - 1, j) =
                    discretization_->t(discretization_->tInteriorIEnd() - 1, j)
                    - discretization_->dx() * settings_.tempBcRight;
        }
    }
};

/**
 * Set the boundary values of the preliminary velocities (u, v)
 * 
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void Computation::applyPreliminaryBoundaryValues() {
    // set boundary values for F at bottom and top side (lower priority)
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for F at bottom side
        if (settings_.outflowBottom) {
            discretization_->f(i, discretization_->uJBegin()) =
                    2 * discretization_->u(i, discretization_->uJBegin())
                    - discretization_->uLast(i, discretization_->uJBegin());
        } else {
            discretization_->f(i, discretization_->uJBegin()) =
                    discretization_->u(i, discretization_->uJBegin());
        }
        // set boundary values for F at top side
        if (settings_.outflowTop) {
            discretization_->f(i, discretization_->uJEnd() - 1) =
                    2 * discretization_->u(i, discretization_->uJEnd() - 1)
                    - discretization_->uLast(i, discretization_->uJEnd() - 1);
        } else {
            discretization_->f(i, discretization_->uJEnd() - 1) =
                    discretization_->u(i, discretization_->uJEnd() - 1);
        }
    }

    // set boundary values for G at bottom and top side (lower priority)
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // set boundary values for v at bottom side
        if (settings_.outflowBottom) {
            discretization_->g(i, discretization_->vJBegin()) =
                    2 * discretization_->v(i, discretization_->vJBegin())
                    - discretization_->vLast(i, discretization_->vJBegin());
        } else {
            discretization_->g(i, discretization_->vJBegin()) =
                    discretization_->v(i, discretization_->vJBegin());
        }
        // set boundary values for v at top side
        if (settings_.outflowBottom) {
            discretization_->g(i, discretization_->vJEnd() - 1) =
                    2 * discretization_->v(i, discretization_->vJEnd() - 1)
                    - discretization_->vLast(i, discretization_->vJEnd() - 1);
        } else {
            discretization_->g(i, discretization_->vJEnd() - 1) =
                    discretization_->v(i, discretization_->vJEnd() - 1);
        }
    }

    // set boundary values for F at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for F at left side
        discretization_->f(discretization_->uIBegin(), j) =
                discretization_->u(discretization_->uIBegin(), j);
        // set boundary values for F at right side
        discretization_->f(discretization_->uIEnd() - 1, j) =
                discretization_->u(discretization_->uIEnd() - 1, j);
    }

    // set boundary values for G at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for G at left side
        if (settings_.outflowLeft) {
            discretization_->g(discretization_->vIBegin(), j) =
                    2 * discretization_->v(discretization_->vIBegin(), j)
                    - discretization_->vLast(discretization_->vIBegin(), j);
        } else {
            discretization_->g(discretization_->vIBegin(), j) =
                    discretization_->v(discretization_->vIBegin(), j);
        }

        // set boundary values for G at right side
        if (settings_.outflowRight) {
            discretization_->g(discretization_->vIEnd() - 1, j) =
                    2 * discretization_->v(discretization_->vIEnd() - 1, j)
                    - discretization_->vLast(discretization_->vIEnd() - 1, j);
        } else {
            discretization_->g(discretization_->vIEnd() - 1, j) =
                    discretization_->v(discretization_->vIEnd() - 1, j);
        }
    }
};

/**
 * Compute the preliminary velocities (F, G) using finite differences
 */
void Computation::computePreliminaryVelocities() {
    // Compute F in the interior of the domain
    for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            double lap_u = discretization_->computeD2uDx2(i, j) + discretization_->computeD2uDy2(i, j);
            double conv_u = discretization_->computeDu2Dx(i, j) + discretization_->computeDuvDy(i, j);
            double f_tilde = discretization_->u(i, j) + dt_ * (lap_u / settings_.re - conv_u + settings_.g[0]);
            double t_interp_right = (discretization_->t(i, j) - discretization_->t(i + 1, j)) / 2.0;
            discretization_->f(i, j) = f_tilde - dt_ * settings_.beta * settings_.g[0] * t_interp_right;
        }
    }

    // Compute G in the interior of the domain
    for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            double lap_v = discretization_->computeD2vDx2(i, j) + discretization_->computeD2vDy2(i, j);
            double conv_v = discretization_->computeDv2Dy(i, j) + discretization_->computeDuvDx(i, j);
            double g_tilde = discretization_->v(i, j) + dt_ * (lap_v / settings_.re - conv_v + settings_.g[1]);
            double t_interp_up = (discretization_->t(i, j) - discretization_->t(i, j + 1)) / 2.0;
            discretization_->g(i, j) = g_tilde - dt_ * settings_.beta * settings_.g[1] * t_interp_up;
        }
    }
};

/**
 * Compute the pressure p by solving the Poisson equation
 */
void Computation::computePressure() {
    pressureSolver_->solve();
};

/**
 * Compute the right hand side rhs of the pressure Poisson equation 
 */
void Computation::computeRightHandSide() {
    // Compute rhs in the interior of the domain using finite differences
    for (int i = discretization_->rhsInteriorIBegin(); i < discretization_->rhsInteriorIEnd(); i++) {
        for (int j = discretization_->rhsInteriorJBegin(); j < discretization_->rhsInteriorJEnd(); j++) {
            double fx = (discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx();
            double gy = (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy();
            discretization_->rhs(i, j) = (fx + gy) / dt_;
        }
    }
};

/**
 * Compute the time step width dt based on the maximum velocities
 */
void Computation::computeTimeStepWidth() {
    // Compute maximal time step width regarding the diffusion
    double dt_diff = settings_.re * settings_.pr / 2 / (1 / (discretization_->dx() * discretization_->dx()) +
                                                        1 / (discretization_->dy() * discretization_->dy()));

    // Compute maximal time step width regarding the convection u
    double dt_conv_u = discretization_->dx() / discretization_->u().absMax();

    // Compute maximal time step width regarding the convection v
    double dt_conv_v = discretization_->dy() / discretization_->v().absMax();

    // Set the appropriate time step width by using a security factor tau
    double computed_dt = settings_.tau * std::min({dt_diff, dt_conv_u, dt_conv_v});
    //std::cout << "dt_diff = " << dt_diff << ", dt_conv_u = " << dt_conv_u << ", dt_conv_v = " << dt_conv_v << std::endl;
    dt_ = std::min({settings_.maximumDt, computed_dt});
};

/**
 * Update the last velocities (uLast, vLast) with the current velocities (u, v)
 */
void Computation::updateLastVelocities() {
    // Update uLast with u
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->uLast(i, j) = discretization_->u(i, j);
        }
    }

    // Update vLast with v
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->vLast(i, j) = discretization_->v(i, j);
        }
    }
};

/**
 * Compute the new velocities (u, v) based on the preliminary velocities (F, G) and the pressure (p)
 */
void Computation::computeVelocities() {
    // Compute u in the interior of the domain
    for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }

    // Compute v in the interior of the domain
    for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
};

void Computation::computeTemperature() {
    for (int i = discretization_->tInteriorIBegin(); i < discretization_->tInteriorIEnd(); i++) {
        for (int j = discretization_->tInteriorJBegin(); j < discretization_->tInteriorJEnd(); j++) {
            discretization_->t(i, j) = discretization_->t(i, j) + dt_ * (1 / (settings_.pr * settings_.re)
                                                                         * (discretization_->computeD2tD2x(i, j) +
                                                                            discretization_->computeD2tD2y(i, j))
                                                                         - discretization_->computeDutDx(i, j) -
                                                                         discretization_->computeDvtDy(i, j)
                                                                         + discretization_->q(i, j));
        }
    }
};
