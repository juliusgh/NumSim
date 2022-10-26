#include "discretization/0_staggered_grid.h"

//! constructor
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells,
                std::array<double, 2> meshWidth)
: nCells_(nCells), meshWidth_(meshWidth)
{
	//TODO: implement
    std::array<double, 2> origin = {0.0, 0.0};
    f_ = FieldVariable(nCells, origin, meshWidth);
    g_ = FieldVariable(nCells, origin, meshWidth);
    p_ = FieldVariable(nCells, origin, meshWidth);
    rhs_ = FieldVariable(nCells, origin, meshWidth);
    u_ = FieldVariable(nCells, origin, meshWidth);
    v_ = FieldVariable(nCells, origin, meshWidth);
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

//! evaluate value of F in an element (i,j)
double StaggeredGrid::f(int i, int j) const
{
	//TODO: implement
    return 0.0;
};

//! evaluate value of G in an element (i,j)
double StaggeredGrid::g(int i, int j) const
{
	//TODO: implement
    return 0.0;
};

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

//! get reference to field variable p
const FieldVariable & StaggeredGrid::p() const
{
	return p_;
};

//! evaluate field variable p in an element (i,j)
double StaggeredGrid::p(int i, int j) const
{
	//TODO: implement
    return 0.0;
};

//! evaluate field variable p in an element (i,j)
double & StaggeredGrid::p(int i, int j)
{
	//TODO: implement
    return 0.0;
};

//! first valid index for p in x direction
int StaggeredGrid::pIBegin() const
{
	//TODO: implement
    return 0;
};

//! one after last valid index for p in x direction
int StaggeredGrid::pIEnd() const
{
	//TODO: implement
    return 0;
};

//! first valid index for p in y direction
int StaggeredGrid::pJBegin() const
{
	//TODO: implement
    return 0;
};

//! one after last valid index for p in y direction
int StaggeredGrid::pJEnd() const
{
	//TODO: implement
    return 0;
};

//! access value of rhs in element (i,j)
double & StaggeredGrid::rhs(int i, int j)
{
	//TODO: implement
    return 0.0;
};

//! get a reference to field variable u
const FieldVariable & StaggeredGrid::u() const
{
	return u_;
};

//! access value of u in element (i,j)
double StaggeredGrid::u(int i, int j) const
{
	//TODO: implement
    return 0.0;
};

//! access value of u in element (i,j)
double & StaggeredGrid::u(int i, int j)
{
	//TODO: implement
    return 0.0;
};

//! first valid index for u in x direction
int StaggeredGrid::uIBegin() const
{
	//TODO: implement
    return 0;
};

//! one after last valid index for u in x direction
int StaggeredGrid::uIEnd() const
{
	//TODO: implement
    return 0;
};

//! first valid index for u in y direction
int StaggeredGrid::uJBegin() const
{
	//TODO: implement
    return 0;
};

//! one after last valid index for u in y direction
int StaggeredGrid::uJEnd() const
{
	//TODO: implement
    return 0;
};

//! get a reference to field variable v
const FieldVariable & StaggeredGrid::v() const
{
    return v_;
};

//! access value of v in element (i,j)
double StaggeredGrid::v(int i, int j) const
{
	//TODO: implement
    return 0.0;
};

//! access value of v in element (i,j)
double & StaggeredGrid::v(int i, int j)
{
	//TODO: implement
    return 0.0;
};

//! first valid index for v in x direction
int StaggeredGrid::vIBegin() const
{
	//TODO: implement
    return 0;
};

//! one after last valid index for v in x direction
int StaggeredGrid::vIEnd() const
{
	//TODO: implement
    return 0;
};

//! first valid index for v in y direction
int StaggeredGrid::vJBegin() const
{
	//TODO: implement
    return 0;
};

//! one after last valid index for v in y direction
int StaggeredGrid::vJEnd() const
{
	//TODO: implement
    return 0;
};
