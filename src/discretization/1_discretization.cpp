#include "1_discretization.h

double Discretization::computeDpDx(int i, int j) const {
    return (p_(i+1, j) - p_(i, j)) / dx();
};