#include <cmath>
#include "discretization/2_donor_cell.h"

  /**
   * calculate derivatives needed for pressure calculations using the donor cell approach
   * @param nCells: number of inner cells
   * @param meshWidth: width of a cell in both directions
   * @param alpha: donor cell weight parameter
   */

DonorCell::DonorCell(std::shared_ptr<Partitioning> partitioning, std::array<double, 2> meshWidth, double alpha, double gamma) :
    Discretization(partitioning, meshWidth),
    alpha_(alpha),
    gamma_(gamma)
{
    
}

/**
 * compute the 1st derivative ∂ u^2 / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return donor cell derivative approximation of the derivative stated above
 */
double DonorCell::computeDu2Dx(int i, int j) const {
    const double u_interp_right = (u(i+1,j) + u(i,j)) / 2.0;
    const double u_interp_left = (u(i,j) + u(i-1,j)) / 2.0;

    const double u_interp_right_donor = (u(i,j) - u(i+1,j)) / 2.0;
    const double u_interp_left_donor = (u(i-1,j) - u(i,j)) / 2.0;

    const double A = (pow(u_interp_right, 2) - pow(u_interp_left, 2)) / dx();
    const double B = (fabs(u_interp_right) * u_interp_right_donor - fabs(u_interp_left) * u_interp_left_donor) / dx();

    return A + alpha_ * B;
}

/**
 * compute the 1st derivative ∂ v^2 / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return donor cell derivative approximation of the derivative stated above
 */
double DonorCell::computeDv2Dy(int i, int j) const {
    const double v_interp_up = (v(i,j+1) + v(i,j)) / 2.0;
    const double v_interp_down = (v(i,j) + v(i,j-1)) / 2.0;

    const double v_interp_up_donor = (v(i,j) - v(i,j+1)) / 2.0;
    const double v_interp_down_donor = (v(i,j-1) - v(i,j)) / 2.0;

    const double A = (pow(v_interp_up, 2) - pow(v_interp_down, 2)) / dy();
    const double B = (v_interp_up_donor * fabs(v_interp_up) - v_interp_down_donor * fabs(v_interp_down)) / dy();

    return A + alpha_ * B;
}

/**
 * compute the 1st derivative ∂ (uv) / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return donor cell derivative approximation of the derivative stated above
 */
double DonorCell::computeDuvDx(int i, int j) const {
    const double u_interp_up = (u(i,j+1) + u(i,j)) / 2.0;
    const double uLeft_interp_up =  (u(i-1,j+1) + u(i-1,j)) / 2.0;

    const double v_interp_right = (v(i,j) + v(i+1,j)) / 2.0;
    const double v_interp_left = (v(i-1,j) + v(i,j)) / 2.0;

    const double v_interp_right_donor = (v(i,j) - v(i+1,j)) / 2.0;
    const double v_interp_left_donor = (v(i-1,j) - v(i,j)) / 2.0;

    const double A = (u_interp_up * v_interp_right - uLeft_interp_up * v_interp_left) / dx();
    const double B = (fabs(u_interp_up) * v_interp_right_donor - fabs(uLeft_interp_up) * v_interp_left_donor) / dx();

    return A + alpha_ * B;
}

/**
 * compute the 1st derivative ∂ (uv) / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return donor cell derivative approximation of the derivative stated above
 */
double DonorCell::computeDuvDy(int i, int j) const {
    const double u_interp_up = (u(i,j+1) + u(i,j)) / 2.0;
    const double u_interp_down = (u(i,j) + u(i,j-1)) / 2.0;

    const double v_interp_right = (v(i+1,j) + v(i,j)) / 2.0;
    const double vDown_interp_right = (v(i+1,j-1) + v(i,j-1)) / 2.0;

    const double u_interp_up_donor = (u(i,j) - u(i,j+1)) / 2.0;
    const double u_interp_down_donor = (u(i,j-1) - u(i,j)) / 2.0;

    const double A = (u_interp_up * v_interp_right - u_interp_down * vDown_interp_right) / dy();
    const double B = (u_interp_up_donor * fabs(v_interp_right) -  u_interp_down_donor * fabs(vDown_interp_right)) / dy();

    return A + alpha_ * B;
}

/**
 * compute the 1st derivative ∂ (ut) / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double DonorCell::computeDutDx(int i, int j) const {
    const double t_interp_right = (t(i+1,j) + t(i,j)) / 2.0;
    const double t_interp_left = (t(i,j) + t(i-1,j)) / 2.0;

    const double t_interp_right_donor = (t(i,j) - t(i+1,j)) / 2.0;
    const double t_interp_left_donor = (t(i-1,j) - t(i,j)) / 2.0;

    const double A = (u(i,j) * t_interp_right - u(i-1,j) * t_interp_left) / dx();
    const double B = (fabs(u(i,j)) * t_interp_right_donor - fabs(u(i-1,j)) * t_interp_left_donor) / dx();

    return A + gamma_ * B;
}

/**
 * compute the 1st derivative ∂ (vt) / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double DonorCell::computeDvtDy(int i, int j) const {
    const double t_interp_up = (t(i, j + 1) + t(i, j)) / 2.0;
    const double t_interp_down = (t(i, j) + t(i, j - 1)) / 2.0;

    const double t_interp_up_donor = (t(i, j) - t(i, j+1)) / 2.0;
    const double t_interp_down_donor = (t(i, j-1) - t(i, j)) / 2.0;

    const double A = (v(i, j) * t_interp_up - v(i, j - 1) * t_interp_down) / dy();
    const double B = (fabs(v(i,j)) * t_interp_up_donor - fabs(v(i,j-1)) * t_interp_down_donor) / dy();

    return A + gamma_ * B;
}