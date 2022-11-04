#include "discretization/2_central_differences.h"

CentralDifferences::CentralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth) :
    Discretization(nCells, meshWidth)
{

}

//! compute the 1st derivative ∂ u^2 / ∂x
double CentralDifferences::computeDu2Dx(int i, int j) const {
    return (u(i+1, j) * u(i+1, j) + 2 * u(i, j) * (u(i+1, j) - u(i-1, j)) - u(i-1, j) * u(i-1, j) ) / (4 * Discretization::dx());
}

//! compute the 1st derivative ∂ v^2 / ∂y
double CentralDifferences::computeDv2Dy(int i, int j) const {
    return (v(i+1, j) * v(i+1, j) + 2 * v(i, j) * (v(i+1, j) - v(i-1, j)) - v(i-1, j) * v(i-1, j) ) / (4 * Discretization::dy());
}

//! compute the 1st derivative ∂ (uv) / ∂x
double CentralDifferences::computeDuvDx(int i, int j) const {
    return ((v(i, j) + v(i+1, j)) * (u(i, j+1) + u(i, j)) - (v(i-1, j) + v(i, j)) * (u(i-1, j+1) + u(i-1, j))) / (4 * Discretization::dx());
}

//! compute the 1st derivative ∂ (uv) / ∂y
double CentralDifferences::computeDuvDy(int i, int j) const {
    return ((u(i, j+1) + u(i, j)) * (v(i+1, j) + v(i, j)) - (u(i, j) + u(i, j-1)) * (v(i+1, j-1) + v(i, j-1))) / (4 * Discretization::dy());
}
