#include "pressure_solver/sor.h"

SOR::SOR(std::shared_ptr<Discretization> discretization,
                double epsilon,
                int maximumNumberOfIterations,
                double omega) :
discretization_(discretization),
epsilon_(epsilon),
maximumNumberOfIterations_(maximumNumberOfIterations),
omega_(omega);

SOR::solve() {

};