#include "storage/fieldvariable.h"
#include "storage/array2d.h"
#include <cmath>
#include <cassert>

/**
 * A field variable is the discretization of a scalar function f(x) with x in the computational domain.
 *  More specifically, a scalar value is stored at discrete nodes/points. The nodes are arranged in an equidistant mesh
 *  with specified mesh width.
 * @param size: number of cells
 * @param origin: origin of the coordinate system
 * @param meshWidth: width of cells in both directions
 */

FieldVariable::FieldVariable(std::array<int, 2> size,
                             std::array<double, 2> origin,
                             std::array<double, 2> meshWidth) :
    Array2D(size),
    origin_(origin),
    meshWidth_(meshWidth)
{
    
}

/**
 * Bilinear interpolation in the interior of the domain
 *
 * The left bottom corner of the domain has the coordinates x = 0, y = 0
 * @param x: x coordinate of the desired point to interpolate
 * @param y: y coordinate of the desired point to interpolate
 * @return interpolated value at the specified point (x, y)
 */
double FieldVariable::interpolateAt(double x, double y) const
{
    // Assert that the specified point is part of the domain
    #ifndef NDEBUG
    assert((0.0 <= x) && (x <= size_[0] * meshWidth_[0]));
    assert((0.0 <= y) && (y <= size_[1] * meshWidth_[1]));
    #endif

    // Determine i and j indices of the corresponding cell (shifted by origin)
    int i = (x - origin_[0]) / meshWidth_[0] + 1;
    int j = (y - origin_[1]) / meshWidth_[1] + 1;

    // Special case: If we are on the upper of right boundary, use the cell in the interior
    if (i == size_[0] - 1)
        i--;
    if (j == size_[1] - 1)
        j--;

    // Obtain the values of the four neighbouring interpolation points
    double valueLeftBottom = (*this)(i, j);
    double valueLeftTop = (*this)(i, j + 1);
    double valueRightBottom = (*this)(i + 1, j);
    double valueRightTop = (*this)(i + 1, j + 1);

    // Determine the coordinates of the four neighbouring interpolation points
    double xLeft = origin_[0] + (i - 1) * meshWidth_[0];
    double xRight = xLeft + meshWidth_[0];
    double yBottom = origin_[1] + (j - 1) * meshWidth_[1];
    double yTop = yBottom + meshWidth_[1];

    /*
     * Bilinear interpolation:
     * We implement it as a repeated linear interpolation (first along x axis, then along y axis)
     * See also https://en.wikipedia.org/wiki/Bilinear_interpolation#Repeated_linear_interpolation
     */
    // 1) Use linear interpolation in x between left and right edge (each on the bottom and top edge)
    double interpBottom = ((xRight - x) * valueLeftBottom + (x - xLeft) * valueRightBottom) / (xRight - xLeft);
    double interpTop = ((xRight - x) * valueLeftTop + (x - xLeft) * valueRightTop) / (xRight - xLeft);
    // 2) Use linear interpolation in y between bottom and top edge
    double interp = ((yTop - y) * interpBottom + (y - yBottom) * interpTop) / (yTop - yBottom);

    return interp;
}
/**
 * Compute absolute maximal value needed in computation to determine optimal time step
 * @return maximal absolute value of the field variable given
 */
double FieldVariable::absMax() const
{
    double abs_max = 0;
    for (int i = 0; i < size_[0]; i++) {
        for (int j = 0; j < size_[1]; j++) {
            if (fabs((*this)(i,j)) > abs_max)
                abs_max = fabs((*this)(i,j));
        }
    }
    return abs_max;
};
