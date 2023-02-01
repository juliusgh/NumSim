#pragma once

#include <iostream>
#include <array>

using namespace std;

/** All settings that parametrize a simulation run.
 */
struct Settings {
    array<int, 2> nCells;          //< number of cells in x and y direction
    array<double, 2> physicalSize; //< physical size of the domain
    double re = 1000;                  //< Reynolds number
    double pr = 1;                   //< Prandtl number
    double beta = 0;                  //< volume expansion coefficient
    double endTime = 10.0;             //< end time of the simulation
    double tau = 0.5;                  //< safety factor for time step width
    double maximumDt = 0.1;            //< maximum time step width

    array<double, 2> g{0., 0.};    //< external forces

    bool useDonorCell = false;         //< if the donor cell scheme should be used
    double alpha = 0.5;                //< factor for donor-cell scheme
    double gamma = alpha;                //< factor for donor-cell for temperature

    // velocity boundary conditions
    array<double, 2> dirichletBcBottom;  //< prescribed values of u,v at bottom of domain
    array<double, 2> dirichletBcTop;     //< prescribed values of u,v at top of domain
    array<double, 2> dirichletBcLeft;    //< prescribed values of u,v at left of domain
    array<double, 2> dirichletBcRight;   //< prescribed values of u,v at right of domain
    bool outflowBottom = false;
    bool outflowTop = false;
    bool outflowLeft = false;
    bool outflowRight = false;

    // temperature initial condition
    double initialTemp = 0; // in Kelvin, 293K = 20°C

    // heat source magnitude
    double heatMagnitude = 10.0;

    // temperature boundary conditions
    double obstacleHot = 1.0;
    double obstacleCold = 0.0;
    bool setFixedTempBottom = false;
    double tempBcBottom = 0.0;  //< prescribed values of tb or tn at bottom of domain
    bool setFixedTempTop = false;
    double tempBcTop = 0.0;  //< prescribed values of tb or tn at top of domain
    bool setFixedTempLeft = false;
    double tempBcLeft = 0.0;  //< prescribed values of tb or tn at left of domain
    bool setFixedTempRight = false;
    double tempBcRight = 0.0;  //< prescribed values of tb or tn at right of domain

    string pressureSolver = "SOR";      //< which pressure solver to use, "GaussSeidel" or "SOR"
    double omega = 1.0;                //< overrelaxation factor
    double epsilon = 1e-5;             //< tolerance for the residual in the pressure solver
    int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver

    string domainfile_path;
    //! parse a text file with settings, each line contains "<parameterName> = <value>"
    void loadFromFile(string filename);

    //! output all settings to console
    void printSettings();
};
