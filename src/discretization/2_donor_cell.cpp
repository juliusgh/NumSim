#include "discretization/2_donor_cell.h"

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha) :
Discretization(nCells, meshWidth)
{

}

//! compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j) const {
    return 0;
}

//! compute the 1st derivative ∂ v^2 / ∂y
double DonorCell::computeDv2Dy(int i, int j) const {
    return 0;
}

//! compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const {
    return 0;
}

//! compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const {
    return 0;
}