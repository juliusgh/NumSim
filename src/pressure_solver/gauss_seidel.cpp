#include "pressure_solver/gauss_seidel.h"

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization,
                double epsilon,
                int maximumNumberOfIterations) :
discretization_(discretization),
epsilon_(epsilon),
maximumNumberOfIterations_(maximumNumberOfIterations);

GaussSeidel::solve() {

};