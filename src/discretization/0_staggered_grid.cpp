#include "discretization/0_staggered_grid.h"
#include <cassert>
#include <iostream>

//! constructor
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells,
                             std::array<double, 2> meshWidth) :
    nCells_(nCells),
    meshWidth_(meshWidth),
    f_({nCells[0] + 1, nCells[1]}, {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
    g_({nCells[0], nCells[1] + 1}, {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
    p_({nCells[0] + 2, nCells[1] + 2}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    rhs_({nCells[0] + 2, nCells[1] + 2}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    u_({nCells[0] + 1, nCells[1] + 2}, {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
    v_({nCells[0] + 2, nCells[1] + 1}, {meshWidth[0] / 2., meshWidth[1]}, meshWidth)
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
 * pressure variable p
 */

//! first valid index for p in x direction
int StaggeredGrid::pIBegin() const
{
    return 0;
};

//! one after last valid index for p in x direction
int StaggeredGrid::pIEnd() const
{
    return nCells_[0] + 1;
};

//! first valid index for p in y direction
int StaggeredGrid::pJBegin() const
{
    return 0;
};

//! one after last valid index for p in y direction
int StaggeredGrid::pJEnd() const
{
    return nCells_[1] + 1;
};
//! first valid Interior index for p in x direction
int StaggeredGrid::pInteriorIBegin() const
{
    return pIBegin()+1;
};

//! one after last valid Interior index for p in x direction
int StaggeredGrid::pInteriorIEnd() const
{
    return pIEnd()-1;
};

//! first valid Interior index for p in y direction
int StaggeredGrid::pInteriorJBegin() const
{
    return pJBegin()+1;
};

//! one after last valid Interior index for p in y direction
int StaggeredGrid::pInteriorJEnd() const
{
    return pJEnd()-1;
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
    return p_(i, j);
};

//! evaluate field variable p in an element (i,j)
double &StaggeredGrid::p(int i, int j)
{
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
    return p_(i, j);
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
    return 0;
};

//! one after last valid index for u in y direction
int StaggeredGrid::uJEnd() const
{
    return nCells_[1] + 1;
};
//! first valid field index for u in x direction
int StaggeredGrid::uInteriorIBegin() const
{

    return uIBegin()+1;
};

//! one after last valid Interior index for u in x direction
int StaggeredGrid::uInteriorIEnd() const
{

    return uIEnd()-1;

};

//! first valid Interior index for u in y direction
int StaggeredGrid::uInteriorJBegin() const
{

    return uJBegin()+1;

};

//! one after last valid Interior index for u in y direction
int StaggeredGrid::uInteriorJEnd() const
{

    return uJEnd()-1;

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
    return u_(i, j);
};

//! access value of u in element (i,j)
double &StaggeredGrid::u(int i, int j)
{
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    return u_(i, j);
};

/*
 * velocity in y-direction v
 */

//! first valid index for v in x direction
int StaggeredGrid::vIBegin() const
{
    return 0;
};

//! one after last valid index for v in x direction
int StaggeredGrid::vIEnd() const
{
    return nCells_[0] + 1;
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

//! first valid Interior index for v in x direction
int StaggeredGrid::vInteriorIBegin() const
{
    return vIBegin()+1;

};

//! one after last valid Interior index for v in x direction
int StaggeredGrid::vInteriorIEnd() const
{
    return vIEnd()-1;

};

//! first valid Interior index for v in y direction
int StaggeredGrid::vInteriorJBegin() const
{
    return vJBegin()+1;

};

//! one after last valid Interior index for v in y direction
int StaggeredGrid::vInteriorJEnd() const
{
    return vJEnd()-1;

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
    return v_(i, j);
};

//! access value of v in element (i,j)
double &StaggeredGrid::v(int i, int j)
{
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    return v_(i, j);
};

/*
 * right hand side rhs
 */

//! first valid index for rhs in x direction
int StaggeredGrid::rhsIBegin() const
{
    return 1;
};

//! one after last valid index for rhs in x direction
int StaggeredGrid::rhsIEnd() const
{
    return nCells_[0];
};

//! first valid index for rhs in y direction
int StaggeredGrid::rhsJBegin() const
{
    return 1;
};

//! one after last valid index for rhs in y direction
int StaggeredGrid::rhsJEnd() const
{
    return nCells_[1];
};

//! access value of rhs in element (i,j)
double StaggeredGrid::rhs(int i, int j) const
{
    assert((rhsIBegin() <= i) && (i <= rhsIEnd()));
    assert((rhsJBegin() <= j) && (j <= rhsJEnd()));
    return rhs_(i, j);
};

//! access value of rhs in element (i,j)
double &StaggeredGrid::rhs(int i, int j)
{
    assert((rhsIBegin() <= i) && (i <= rhsIEnd()));
    assert((rhsJBegin() <= j) && (j <= rhsJEnd()));
    return rhs_(i, j);
};

//! get reference to field variable rhs
const FieldVariable &StaggeredGrid::rhs() const
{
    return rhs_;
};

/*
 * preliminary velocity F
 */

//! access value of F in element (i,j)
double StaggeredGrid::f(int i, int j) const
{
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    return f_(i, j);
};

//! access value of F in element (i,j)
double &StaggeredGrid::f(int i, int j)
{
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    return f_(i, j);
};

//! get reference to field variable F
const FieldVariable &StaggeredGrid::f() const
{
    return f_;
};

/*
 * preliminary velocity G
 */

//! access value of G in element (i,j)
double StaggeredGrid::g(int i, int j) const
{
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    return g_(i, j);
};

//! access value of G in element (i,j)
double &StaggeredGrid::g(int i, int j)
{
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    return g_(i, j);
};

//! get reference to field variable G
const FieldVariable &StaggeredGrid::g() const
{
    return g_;
};
