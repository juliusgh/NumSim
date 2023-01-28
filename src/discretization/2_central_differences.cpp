#include "discretization/2_central_differences.h"
#include <cmath>

/**
 *  calculate derivatives needed for pressure calculations using the central differences approach
 * @param nCells: number of inner cells
 * @param meshWidth: width of cells in all directions
 */

CentralDifferences::CentralDifferences(std::shared_ptr<Partitioning> partitioning,
                                       std::array<double, 2> meshWidth,
                                       std::shared_ptr<Settings> settings) :
        Discretization(partitioning, meshWidth, settings) {

}

/**
 * compute the 1st derivative ∂ u^2 / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDu2Dx(int i, int j) const {
    const double u_interp_right = (u(i + 1, j) + u(i, j)) / 2.0;
    const double u_interp_left = (u(i, j) + u(i - 1, j)) / 2.0;

    return (pow(u_interp_right, 2) - pow(u_interp_left, 2)) / dx();
}

/**
 * compute the 1st derivative ∂ v^2 / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDv2Dy(int i, int j) const {
    const double v_interp_up = (v(i, j + 1) + v(i, j)) / 2.0;
    const double v_interp_down = (v(i, j) + v(i, j - 1)) / 2.0;

    return (pow(v_interp_up, 2) - pow(v_interp_down, 2)) / dy();
}

/**
 * compute the 1st derivative ∂ (uv) / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDuvDx(int i, int j) const {
    const double u_interp_up = (u(i, j + 1) + u(i, j)) / 2.0;
    const double uLeft_interp_up = (u(i - 1, j + 1) + u(i - 1, j)) / 2.0;

    const double v_interp_right = (v(i, j) + v(i + 1, j)) / 2.0;
    const double v_interp_left = (v(i - 1, j) + v(i, j)) / 2.0;

    return (u_interp_up * v_interp_right - uLeft_interp_up * v_interp_left) / dx();
}

/**
 * compute the 1st derivative ∂ (uv) / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return donor cell derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDuvDy(int i, int j) const {
    const double u_interp_up = (u(i, j + 1) + u(i, j)) / 2.0;
    const double u_interp_down = (u(i, j) + u(i, j - 1)) / 2.0;

    const double v_interp_right = (v(i + 1, j) + v(i, j)) / 2.0;
    const double vDown_interp_right = (v(i + 1, j - 1) + v(i, j - 1)) / 2.0;

    return (u_interp_up * v_interp_right - u_interp_down * vDown_interp_right) / dy();
}

/**
 * compute the 1st derivative ∂ (ut) / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDutDx(int i, int j) const {
    const double t_interp_right = (t(i + 1, j) + t(i, j)) / 2.0;
    const double t_interp_left = (t(i, j) + t(i - 1, j)) / 2.0;

    return (u(i, j) * t_interp_right - u(i - 1, j) * t_interp_left) / dx();
}

/**
 * compute the 1st derivative ∂ (vt) / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDvtDy(int i, int j) const {
    const double t_interp_up = (t(i, j + 1) + t(i, j)) / 2.0;
    const double t_interp_down = (t(i, j) + t(i, j - 1)) / 2.0;

    return (v(i, j) * t_interp_up - v(i, j - 1) * t_interp_down) / dy();
}

