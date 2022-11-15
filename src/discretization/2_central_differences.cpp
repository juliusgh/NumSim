#include "discretization/2_central_differences.h"

/**
 *  calculate derivatives needed for pressure calculations using the central differences approach
 * @param nCells: number of inner cells
 * @param meshWidth: width of cells in all directions
 */

CentralDifferences::CentralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth) :
    Discretization(nCells, meshWidth)
{

}

/**
 * compute the 1st derivative ∂ u^2 / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDu2Dx(int i, int j) const {
    return ((u(i+1,j) + u(i,j)) * (u(i+1,j) + u(i,j)) - (u(i,j) + u(i-1,j)) * (u(i,j) + u(i-1,j))) / (4.0 * dx());
}

/**
 * compute the 1st derivative ∂ v^2 / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDv2Dy(int i, int j) const {
    return ((v(i,j+1) + v(i,j)) * (v(i,j+1) + v(i,j)) - (v(i,j) + v(i,j-1)) * (v(i,j) + v(i,j-1))) / (4.0 * dy());
}

/**
 * compute the 1st derivative ∂ (uv) / ∂x
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return central differences derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDuvDx(int i, int j) const {
    return ((v(i,j) + v(i+1,j)) * (u(i,j+1) + u(i,j)) - (v(i-1,j) + v(i,j)) * (u(i-1,j+1) + u(i-1,j))) / (4.0 * dx());
}

/**
 * compute the 1st derivative ∂ (uv) / ∂y
 * @param i: discretized position in x direcetion
 * @param j: discretiszed position in y direction
 * @return donor cell derivative approximation of the derivative stated above
 */
double CentralDifferences::computeDuvDy(int i, int j) const {
    return ((u(i,j+1) + u(i,j)) * (v(i+1,j) + v(i,j)) - (u(i,j) + u(i,j-1)) * (v(i+1,j-1) + v(i,j-1))) / (4.0 * dy());
}
