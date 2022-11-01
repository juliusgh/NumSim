#include "1_discretization.h"

Discretization::Discretization(std::array<int,2> nCells, std::array<double,2> meshWidth) :
        StaggeredGrid(nCells, meshWidth)
{

};

double Discretization::computeD2uDx2(int i, int j) const {
    return 0.0;//(p_(i+1, j) - p_(i, j)) / Discretization::dx();
};

double Discretization::computeD2uDy2(int i, int j) const {
    return 0.0;//(p_(i+1, j) - p_(i, j)) / Discretization::dx();
};

double Discretization::computeD2vDx2(int i, int j) const {
    return 0.0;//(p_(i+1, j) - p_(i, j)) / Discretization::dx();
};

double Discretization::computeD2vDy2(int i, int j) const {
    return 0.0;//(p_(i+1, j) - p_(i, j)) / Discretization::dx();
};

double Discretization::computeDpDx(int i, int j) const {
    return 0.0;//(p_(i+1, j) - p_(i, j)) / Discretization::dx();
};

double Discretization::computeDpDy(int i, int j) const {
    return 0.0;//(p_(i, j+1) - p_(i, j)) / Discretization::dy();
};

double Discretization::computeDu2Dx(int i, int j) const {
    return 0.0;//(p_(i, j+1) - p_(i, j)) / Discretization::dy();
};

double Discretization::computeDuvDx(int i, int j) const {
    return 0.0;//(p_(i, j+1) - p_(i, j)) / Discretization::dy();
};

double Discretization::computeDuvDy(int i, int j) const {
    return 0.0;//(p_(i, j+1) - p_(i, j)) / Discretization::dy();
};

double Discretization::computeDv2Dy(int i, int j) const {
    return 0.0;//(p_(i, j+1) - p_(i, j)) / Discretization::dy();
};