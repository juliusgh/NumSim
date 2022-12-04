#pragma once

#include "pressure_solver/1_red_black.h"

/**
 * Parallel red black solver for solving a linear system of equations.
 */

class RedBlackSOR : public RedBlack {
public:
    /**
    * constructor
    * @param discretization
    * @param epsilon
    * @param maximumNumberOfIterations
    * @param omega
    */
    RedBlackSOR(std::shared_ptr <Discretization> discretization,
        double epsilon,
        int maximumNumberOfIterations,
        double omega,
        std::shared_ptr<Partitioning> partitioning);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    void solve() override;
private:
    double omega_;
};
