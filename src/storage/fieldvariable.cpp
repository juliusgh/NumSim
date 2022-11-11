#include "storage/fieldvariable.h"
#include "storage/array2d.h"
#include <cassert>
#include <cmath>
#include <iostream>

FieldVariable::FieldVariable(std::array<int, 2> size,
                             std::array<double, 2> origin,
                             std::array<double, 2> meshWidth) :
    Array2D(size),
    origin_(origin),
    meshWidth_(meshWidth)
{
    
}

double FieldVariable::interpolateAt(double x, double y) const
{
    // FIXME: probably does not work
    // TODO: is this assert too strict (boundary) ?
    assert((0.0 <= x) && (x <= size_[0] * meshWidth_[0]));
    assert((0.0 <= y) && (y <= size_[1] * meshWidth_[1]));

    std::cout << x << "," << y;

    x += meshWidth_[0];
    y += meshWidth_[1];

    // determine i and j indices
    int i = (x - origin_[0]) / meshWidth_[0];
    int j = (y - origin_[1]) / meshWidth_[1];

    //std::cout << "i:" << i << ", j:" << j << std::endl;
    if (i == size_[0] - 1)
        i--;
    if (j == size_[1] - 1)
        j--;

    std::cout << " results in field "<< i << "," << j;

    double q00 = (*this)(i, j);
    double q01 = (*this)(i, j + 1);
    double q10 = (*this)(i + 1, j);
    double q11 = (*this)(i + 1, j + 1);

    // determine interpolation points
    i--;
    j--;
    double x0 = origin_[0] + i * meshWidth_[0];
    double x1 = x0 + meshWidth_[0];
    double y0 = origin_[1] + j * meshWidth_[1];
    double y1 = y0 + meshWidth_[1];

    std::cout << " interpolation points: " << x0 << "," <<  x1 <<  "," << y0 <<  "," << y1;

    std::cout << " interpolation values: " << q00 << "," <<  q01 <<  "," << q10 <<  "," << q11<< std::endl;

    x -= meshWidth_[0];
    y -= meshWidth_[1];

    // bilinear interpolation
    double f0 = (x1 - x) / (x1 - x0) * q00 + (x - x0) / (x1 - x0) * q10;
    double f1 = (x1 - x) / (x1 - x0) * q01 + (x - x0) / (x1 - x0) * q11;
    double f = (y1 - y) / (y1 - y0) * f0 + (y - y0) / (y1 - y0) * f1;
    return f;

    /*const double dx = meshWidth_[0]; // mesh width in x dir.
    const double dy = meshWidth_[1]; // mesh width in y dir.

    // indicies of (cell (i,j) in which the point (x,y) lies
    int iLeftEdge  = (int) std::floor((x + origin_[0]) / dx); // -1 for  neg idx??
    int jLowerEdge = (int) std::floor((y + origin_[1]) / dy);



    //  shift right and upper boundaries so that they don't use cells outside of the grid
    if (iLeftEdge >= (*this).size_[0] -1)
    {
        iLeftEdge = iLeftEdge -1; // shift it one column to the right
    }
    if (jLowerEdge >= (*this).size_[1] -1)
    {
        jLowerEdge = jLowerEdge -1; // shift it down
    }



    // relative position of x and y in the cell
    // one cell: |<-xr1-> x <-xr2->|
    //           |<--    dx      ->|
    const double xr1 = x  - ((meshWidth_[0]*iLeftEdge) - origin_[0]);   // relative position of x from left edge
    const double yr1 = y  - ((meshWidth_[1]*jLowerEdge) - origin_[1]);   // relative poistion of y from lower edge
    const double xr2 = dx - xr1;   // distance right_edge - x
    const double yr2 = dy - yr1;   // distance upper edge - y

    // transform to x, y coordinates when directly accessing the array2D
    int xLeftEdge = iLeftEdge;
    int yLowerEdge = jLowerEdge;

    // get values at corner points
    const double f_lowerLeft  = (*this)(xLeftEdge,     yLowerEdge);
    const double f_upperLeft  = (*this)(xLeftEdge,     yLowerEdge + 1);
    const double f_lowerRight = (*this)(xLeftEdge + 1, yLowerEdge);
    const double f_upperRight = (*this)(xLeftEdge + 1, yLowerEdge + 1);

    // bilinear interpolation
    const double f_intp = (f_lowerLeft * xr2 * yr2
                           + f_lowerRight * xr1 * yr2
                           + f_upperLeft  * xr2 * yr1
                           + f_upperRight * xr1 * yr1) / (dx * dy);

    return f_intp;*/
}

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
