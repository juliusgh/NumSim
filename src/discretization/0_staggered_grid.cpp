#include "discretization/0_staggered_grid.h"
#include <cassert>

//! constructor
// TODO: initialization of f, g and rhs is probably not correct yet
// TODO: is origin initialization correct?
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells,
                             std::array<double, 2> meshWidth) :
    nCells_(nCells),
    meshWidth_(meshWidth),
    f_(nCells, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    g_(nCells, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    p_({nCells[0] + 1, nCells[1] + 1}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    rhs_(nCells, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    u_({nCells[0] + 1, nCells[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    v_({nCells[0], nCells[1] + 1}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth)
{
    
};

/*
 * mesh information
 */

//! get the mesh width
const std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
};

//! get number of cells
const std::array<int, 2> StaggeredGrid::nCells() const
{
    return nCells_;
};

//! get mesh width in x-direction
double StaggeredGrid::dx() const
{
    return meshWidth_[0];
};

//! get mesh width in y-direction
double StaggeredGrid::dy() const
{
    return meshWidth_[1];
};

/*
 * TODO: what are F and G?
 */

//! evaluate value of F in an element (i,j)
double StaggeredGrid::f(int i, int j) const
{
    // TODO: implement
    return f_(i, j);
};

//! evaluate value of G in an element (i,j)
double StaggeredGrid::g(int i, int j) const
{
    // TODO: implement
    return g_(i, j);
};

/*
 * pressure variable p
 */

//! first valid index for p in x direction
int StaggeredGrid::pIBegin() const
{
    return 1;
};

//! one after last valid index for p in x direction
int StaggeredGrid::pIEnd() const
{
    return nCells_[0];
};

//! first valid index for p in y direction
int StaggeredGrid::pJBegin() const
{
    return 1;
};

//! one after last valid index for p in y direction
int StaggeredGrid::pJEnd() const
{
    return nCells_[1];
};

//! get reference to field variable p
const FieldVariable &StaggeredGrid::p() const
{
    return p_;
};

//! evaluate field variable p in an element (i,j)
double StaggeredGrid::p(int i, int j) const
{
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
    return p_(i - 1, j - 1);
};

//! evaluate field variable p in an element (i,j)
double &StaggeredGrid::p(int i, int j)
{
    return StaggeredGrid::p(i, j);
};

/*
 * right hand side rhs
 */

//! access value of rhs in element (i,j)
double &StaggeredGrid::rhs(int i, int j)
{
    // TODO: implement
    double rhs = 0.0;
    return rhs;
};

/*
 * velocity in x-direction u
 */

//! first valid index for u in x direction
int StaggeredGrid::uIBegin() const
{
    return 0;
};

//! one after last valid index for u in x direction
int StaggeredGrid::uIEnd() const
{
    return nCells_[0];
};

//! first valid index for u in y direction
int StaggeredGrid::uJBegin() const
{
    return 1;
};

//! one after last valid index for u in y direction
int StaggeredGrid::uJEnd() const
{
    return nCells_[1];
};

//! get a reference to field variable u
const FieldVariable &StaggeredGrid::u() const
{
    return u_;
};

//! access value of u in element (i,j)
double StaggeredGrid::u(int i, int j) const
{
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    return u_(i, j - 1);
};

//! access value of u in element (i,j)
double &StaggeredGrid::u(int i, int j)
{
    return StaggeredGrid::u(i, j);
};

/*
 * velocity in y-direction v
 */

//! first valid index for v in x direction
int StaggeredGrid::vIBegin() const
{
    return 1;
};

//! one after last valid index for v in x direction
int StaggeredGrid::vIEnd() const
{
    return nCells_[0];
};

//! first valid index for v in y direction
int StaggeredGrid::vJBegin() const
{
    return 0;
};

//! one after last valid index for v in y direction
int StaggeredGrid::vJEnd() const
{
    return nCells_[1];
};

//! get a reference to field variable v
const FieldVariable &StaggeredGrid::v() const
{
    return v_;
};

//! access value of v in element (i,j)
double StaggeredGrid::v(int i, int j) const
{
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    return v_(i - 1, j);
};

//! access value of v in element (i,j)
double &StaggeredGrid::v(int i, int j)
{
    return StaggeredGrid::v(i, j);
};
