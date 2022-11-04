#include "computation/computation.h"


//! initialize the computation object
//! parse the settings from file that is given as the only command line argument
void Computation::initialize(int argc, char *argv[])
{

};

//! run the whole simulation until tend
void Computation::runSimulation() {
    int t_iter = 0;
    double t = 0;
    while (t <= settings_.endTime){
        t_iter++;
        applyBoundaryValues();
        computeTimeStepWidth();
        // set boundary fg
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        //output
        t += dt_;
    }
};
//! set boundary values of u and v to correct values

void Computation::applyBoundaryValues(){

};

//! compute the preliminary velocities, F and G
void Computation::computePreliminaryVelocities(){
    for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++) {
            double lap_u = discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j);
            double conv_u = discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j);
            discretization_->f(i,j) = discretization_->u(i,j) + dt_ * (lap_u / settings_.re - conv_u + settings_.g[0]);
        }
    }
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd(); j++) {
            double lap_v = discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j);
            double conv_v = discretization_->computeDv2Dy(i,j) + discretization_->computeDuvDx(i,j);
            discretization_->g(i,j) = discretization_->v(i,j) + dt_ * (lap_v / settings_.re - conv_v + settings_.g[1]);
        }
    }
};

//! solve the Poisson equation for the pressure
void Computation::computePressure() {
    auto sor = SOR(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    sor.solve();
};

//! compute the right hand side of the Poisson equation for the pressure
void Computation::computeRightHandSide() {
    for (int i = discretization_->rhsIBegin() + 1; i < discretization_->rhsIEnd(); i++) {
        for (int j = discretization_->rhsJBegin() + 1; j < discretization_->rhsJEnd(); j++) {
            double fx = (discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx();
            double gy = (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy();
            discretization_->rhs(i, j) = (fx + gy) / dt_;
        }
    }
};
//! compute the time step width dt from maximum velocities
void Computation::computeTimeStepWidth() {
    double dt_diff = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );
    double dt_conv_u = discretization_->dx() / discretization_->uMax();
    double dt_conv_v = discretization_->dy() / discretization_->vMax();
    dt_ = min(dt_diff, min(dt_conv_u, dt_conv_v));
};

//! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
void Computation::computeVelocities() {
    for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd(); j++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
};
