//
// Created by nhornischer on 28.01.23.
//
#pragma once

#include "pressure_solver/0_pressure_solver.h"

class Multigrid : public PressureSolver {

public:
    Multigrid(std::shared_ptr<Discretization> discretization,
              double epsilon,
              int maximuNumberOfIterations,
              int gamma);

    void solve();

private:
    int gamma_;

    void smoothing();

    void defect_calculation();

    void restriction();

    void prolongation();

    void correction();

    void MG_Cycle(int layer);

};