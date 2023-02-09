#pragma once

#include <memory>
#include <algorithm>
#include "settings.h"
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
    /**
     * Initialize the computation object
     * 
     * Parse the settings from the parameter file that is given as the command line argument
     * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
     */
    virtual void initialize(string filename);

    /**
     * Run the whole simulation until tend
     */
    virtual void runSimulation();

protected:

    /**
     * Set the boundary values of the velocities (u, v)
     * 
     * Left and right boundaries should overwrite bottom and top boundaries
     */

    void applyBoundaryValues();

    void applyBoundaryValuesTop();

    void applyBoundaryValuesBottom();

    void applyBoundaryValuesLeft();

    void applyBoundaryValuesRight();

    /**
     * Set the boundary values of the preliminary velocities (u, v)
     * 
     * Left and right boundaries should overwrite bottom and top boundaries
     */
    virtual void applyPreliminaryBoundaryValues();

    /**
     * Compute the preliminary velocities (F, G) using finite differences
     */
    virtual void computePreliminaryVelocities();

    /**
     * Compute the pressure p by solving the Poisson equation
     */
    virtual void computePressure();

    /**
     * Compute the right hand side rhs of the pressure Poisson equation 
     */
    virtual void computeRightHandSide();

    /**
     * Compute the time step width dt based on the maximum velocities
     */
    virtual void computeTimeStepWidth();

    virtual void updateLastVelocities();

    /**
     * Compute the new velocities (u, v) based on the preliminary velocities (F, G) and the pressure (p)
     */
    virtual void computeVelocities();

    virtual void trackParticles();

    /**
     * Compute new temperature t
     */
    virtual void computeTemperature();

    Settings settings_;
    std::shared_ptr<Discretization> discretization_;
    std::shared_ptr<Partitioning> partitioning_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 2> meshWidth_;
    double dt_;
};