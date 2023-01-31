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
              int maximumNumberOfLayers,
              int gamma,
              int mu,
              int nu,
              double theta);

    void solve() override;

private:
    int gamma_;

    int mu_;

    int nu_;

    int maximumNumberOfLayers_;

    double theta_;

    void smoothing(Array2D& v, int smoothing_steps, bool convergenceCheck = false);


    Array2D defect_calculation(Array2D v_L , Array2D rhs);

    Array2D restriction(Array2D d_L);

    Array2D prolongation(Array2D);

    void correction(Array2D c_L, Array2D& v_L);

    Array2D MG_Cycle(int layer, Array2D& v, Array2D rhs);

};