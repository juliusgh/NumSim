#pragma once
#include "settings.h"
#include <memory>
#include <algorithm>
#include "discretization/1_discretization.h"
#include "discretization/2_donor_cell.h"
#include "discretization/2_central_differences.h"
#include "pressure_solver/0_pressure_solver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"


/** This class handles the main simulation.
* It implements the time stepping scheme, computes all the terms and calls the pressure solver.
*/

class Computation {
public:
    //! initialize the computation object
    //! parse the settings from file that is given as the only command line argument
    void initialize(string filename);

    //! run the whole simulation until tend
    void runSimulation();

private:
    //! set boundary values of u and v to correct values
    void applyBoundaryValues();

    //! set boundary values of F and G to correct values
    void applyPreliminaryBoundaryValues();

    //! compute the preliminary velocities, F and G
    void computePreliminaryVelocities();

    //! solve the Poisson equation for the pressure
    void computePressure();

    //! compute the right hand side of the Poisson equation for the pressure
    void computeRightHandSide();
    //! compute the time step width dt from maximum velocities
    void computeTimeStepWidth();

    //! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
    void computeVelocities();

    Settings settings_;
    std::shared_ptr<Discretization> discretization_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 2> meshWidth_;
    double dt_;
};