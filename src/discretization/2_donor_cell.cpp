#include <cmath>
#include "discretization/2_donor_cell.h"

/**
 *  calculate derivatives needed for pressure calculations using the donor cell approach
 * @param nCells
 * @param meshWidth
 * @param alpha
 */

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha) :
    Discretization(nCells, meshWidth),
    alpha_(alpha)
{
    
}

//! compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j) const {
    double A = ((u(i,j) + u(i+1,j)) * (u(i,j) + u(i+1,j)) - (u(i-1,j) + u(i,j)) * (u(i-1,j) + u(i,j))) / (4.0 * dx());
    double B = (fabs(u(i,j) + u(i+1,j)) * (u(i,j) - u(i+1,j)) - fabs(u(i-1,j) + u(i,j)) * (u(i-1,j) - u(i,j))) / (4.0 * dx());
    return A + alpha_ * B;
}

//! compute the 1st derivative ∂ v^2 / ∂y
double DonorCell::computeDv2Dy(int i, int j) const {
    double A = ((v(i,j+1) + v(i,j)) * (v(i,j+1) + v(i,j)) - (v(i,j) + v(i,j-1)) * (v(i,j) + v(i,j-1))) / (4.0 * dy());
    double B = ((v(i,j) - v(i,j+1)) * fabs(v(i,j+1) + v(i,j)) - (v(i,j-1) - v(i,j)) * fabs(v(i,j) + v(i,j-1))) / (4.0 * dy());
    return A + alpha_ * B;
}

//! compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const {
    double A = ((v(i+1,j) + v(i,j)) * (u(i,j+1) + u(i,j)) - (v(i,j) + v(i-1,j)) * (u(i-1,j+1) + u(i-1,j))) / (4.0 * dx());
    double B = ((v(i,j) - v(i+1,j)) * fabs(u(i,j+1) + u(i,j)) - (v(i-1,j) - v(i,j)) * fabs(u(i-1,j+1) + u(i-1,j))) / (4.0 * dx());
    return A + alpha_ * B;
}

//! compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const {
    double A = ((u(i,j+1) + u(i,j)) * (v(i+1,j) + v(i,j)) - (u(i,j) + u(i,j-1)) * (v(i+1,j-1) + v(i,j-1))) / (4.0 * dy());
    double B = ((u(i,j) - u(i,j+1)) * fabs(v(i+1,j) + v(i,j)) - (u(i,j-1) - u(i,j)) * fabs(v(i+1,j-1) + v(i,j-1))) / (4.0 * dy());
    return A + alpha_ * B;
}