#include "1_discretization.h"

/**
 * Calculate derivatives needed for velocity calculations
 * @param partitioning: encapsulate functionality corresponding to subdomain handling
 * @param meshWidth: width of grid cell in one direction
 * @param settings: information about the settings received from parameter file
 */

Discretization::Discretization(std::shared_ptr<Partitioning> partitioning,
                               std::array<double, 2> meshWidth,
                               Settings settings) :
        StaggeredGrid(partitioning, meshWidth, settings) {

};

/**
* compute the 2nd derivative ∂^2 u / ∂x^2
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeD2uDx2(int i, int j) const {
    return (u(i + 1, j) - 2.0 * u(i, j) + u(i - 1, j)) / (dx() * dx());
};

/**
* compute the 2nd derivative ∂^2 u / ∂y^2
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeD2uDy2(int i, int j) const {
    return (u(i, j + 1) - 2.0 * u(i, j) + u(i, j - 1)) / (dy() * dy());
};

/**
* compute the 2nd derivative ∂^2 v / ∂x^2
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeD2vDx2(int i, int j) const {
    return (v(i + 1, j) - 2.0 * v(i, j) + v(i - 1, j)) / (dx() * dx());
};

/**
* compute the 2nd derivative ∂^2 v / ∂y^2
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeD2vDy2(int i, int j) const {
    return (v(i, j + 1) - 2.0 * v(i, j) + v(i, j - 1)) / (dy() * dy());
};

/**
* compute the 1st derivative ∂p / ∂x
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeDpDx(int i, int j) const {
    return (p(i + 1, j) - p(i, j)) / dx();
};

/**
* compute the 1st derivative ∂p / ∂y
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeDpDy(int i, int j) const {
    return (p(i, j + 1) - p(i, j)) / dy();
};

/**
* compute the 2nd derivative ∂²t / ∂x²
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeD2tD2x(int i, int j) const {
    return (t(i + 1, j) - 2.0 * t(i, j) + t(i - 1, j)) / (dx() * dx());
};

/**
* compute the 2nd derivative ∂²t / ∂y²
* @param i: discretized position in x direcetion
* @param j: discretized position in y direction
* @return derivative approximation of the derivative stated above
*/
double Discretization::computeD2tD2y(int i, int j) const {
    return (t(i, j + 1) - 2.0 * t(i, j) + t(i, j - 1)) / (dy() * dy());
};
