//
// Created by nhornischer on 28.01.23.
//

#include "1_multigrid.h"

Multigrid::Multigrid(std::shared_ptr<Discretization> discretization,
                     double epsilon,
                     int maximuNumberOfIterations,
                     int gamma):
        PressureSolver(discretization,epsilon,maximuNumberOfIterations),
        gamma_(gamma){

}

void Multigrid::solve(){

};

void Multigrid::smoothing() {

};

void Multigrid::defect_calculation() {

};

void Multigrid::restriction() {

};

void Multigrid::prolongation() {

};

void Multigrid::correction() {

};

void Multigrid::MG_Cycle(int layer) {
    // Presmoothing
    // defect calculation
    // restriction
    //residual equation(MGCycle)
    for (int i = 1; i < gamma_; i++) {
        MG_Cycle(layer - 1);
    }
    //Prolongation
    //Correction
    //Postsmoothing
};