#include "storage/fieldvariable.h"
#include "storage/array2d.h"

FieldVariable::FieldVariable(std::array<int, 2> size,
                             std::array<double, 2> origin,
                             std::array<double, 2> meshWidth)
: Array2D(size)
{
    origin_ = origin;
    meshWidth_ = meshWidth;
}

double FieldVariable::interpolateAt(double x, double y) const
{
    // TODO: implement
}
