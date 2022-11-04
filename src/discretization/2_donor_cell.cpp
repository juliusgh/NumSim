#include "discretization/2_donor_cell.h"

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha) :
    Discretization(nCells, meshWidth),
    alpha_(alpha)
{
    
}

//! compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j) const {
    double A = ((u(i,j) + u(i+1,j)) * (u(i,j) + u(i+1,j)) - (u(i-1,j) + u(i,j)) * (u(i-1,j) + u(i,j))) / (4 * Discretization::dx());
    double B = (abs(u(i,j) + u(i+1,j)) * (u(i,j) - u(i+1,j)) - abs(u(i-1,j) + u(i,j)) * (u(i-1,j) - u(i,j))) / (4 * Discretization::dx());
    return A + DonorCell::alpha_ * B;
}

//! compute the 1st derivative ∂ v^2 / ∂y
double DonorCell::computeDv2Dy(int i, int j) const {
    double A = ((v(i,j+1) + v(i,j)) * (v(i,j+1) + v(i,j)) - (v(i,j) + v(i,j-1)) * (v(i,j) + v(i,j-1))) / (4 * Discretization::dy());
    double B = ((v(i,j) - v(i,j+1)) * abs(v(i,j+1) + v(i,j)) - (v(i,j-1) - v(i,j)) * abs(v(i,j) + v(i,j-1))) / (4 * Discretization::dy());
    return A+ DonorCell::alpha_;
}

//! compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const {
    double A = ((v(i+1,j) + v(i,j)) * (u(i,j+1) + u(i,j)) - (v(i,j) + v(i-1,j)) * (u(i-1,j+1) + u(i-1,j))) / (4 * Discretization::dx());
    double B = ((v(i,j) - v(i+1,j)) * abs(u(i,j+1) + u(i,j)) - (v(i-1,j) - v(i,j)) * abs(u(i-1,j+1) + u(i-1,j))) / (4 * Discretization::dx());
    return A+ DonorCell::alpha_;
}

//! compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const {
    double A = ((u(i,j+1) + u(i,j)) * (v(i+1,j) + v(i,j)) - (u(i,j) + u(i,j-1)) * (v(i+1,j-1) + v(i,j-1))) / (4 * Discretization::dy());
    double B = ((u(i,j) - u(i,j+1)) * abs(v(i+1,j) + v(i,j)) - (u(i,j-1) - u(i,j-1)) * abs(v(i+1,j-1) + v(i,j-1))) / (4 * Discretization::dy());
    return A + DonorCell::alpha_ * B;
}