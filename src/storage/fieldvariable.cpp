#include "storage/fieldvariable.h"
#include "storage/array2d.h"
#include <cassert>
#include <cmath>

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
    // TODO: is this assert too strict (boundary) ?
    assert((0.0 <= x) && (x <= size_[0] * meshWidth_[0]));
    assert((0.0 <= y) && (y <= size_[1] * meshWidth_[1]));

    // determine i and j indices
    int i = (x - origin_[0]) / meshWidth_[0];
    int j = (y - origin_[1]) / meshWidth_[1];
    if (i == size_[0] - 1) {
        i--;
    }
    if (j == size_[1] - 1) {
        j--;
    }

    // determine interpolation points
    double x0 = origin_[0] + i * meshWidth_[0];
    double x1 = x0 + meshWidth_[0];
    double y0 = origin_[1] + j * meshWidth_[1];
    double y1 = y0 + meshWidth_[1];
    double q00 = (*this)(i, j);
    double q01 = (*this)(i, j + 1);
    double q10 = (*this)(i + 1, j);
    double q11 = (*this)(i + 1, j + 1);

    // bilinear interpolation
    double f0 = (x1 - x) / (x1 - x0) * q00 + (x - x0) / (x1 - x0) * q10;
    double f1 = (x1 - x) / (x1 - x0) * q01 + (x - x0) / (x1 - x0) * q11;
    double f = (y1 - y) / (y1 - y0) * f0 + (y - y0) / (y1 - y0) * f1;
    return f;
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
