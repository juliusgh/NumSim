#include <iostream>
#include "discretization/0_staggered_grid.h"

/**
 * Implement staggered grid, providing a variety of parameters
 * @param nCells: number of cells
 * @param meshWidth: cell width in all directions
 */
StaggeredGrid::StaggeredGrid(std::shared_ptr<Partitioning> partitioning,
                             std::array<double, 2> meshWidth) :
    partitioning_(partitioning),
    nCells_(partitioning->nCellsLocal()),
    meshWidth_(meshWidth),
    f_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
    g_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
    p_(pSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    rhs_(rhsSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
    u_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
    v_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
    t_(tSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth)
{

};

/**
 * mesh information
 */

/**
 * get the mesh width
 * @return mesh width
 */
const std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
};

/**
 * get number of cells
 * @return number of cells
 */
const std::array<int, 2> StaggeredGrid::nCells() const
{
    return nCells_;
};

/**
 * get mesh width in x-direction
 * @return mesh width in x-direction
 */
double StaggeredGrid::dx() const
{
    return meshWidth_[0];
};

/**
 * get mesh width in y-direction
 * @return mesh width in y-direction
 */
double StaggeredGrid::dy() const
{
    return meshWidth_[1];
};

/**
 * pressure variable p
 */

/**
 * get first valid index for p in x direction
 * @return first valid index for p in x direction
 */
int StaggeredGrid::pIBegin() const
{
    return -1;
};

/**
 * get one after last valid index for p in x direction
 * @return one after last valid index for p in x direction
 */
int StaggeredGrid::pIEnd() const
{
    return nCells_[0] + 1;
};

/**
 * get first valid index for p in y direction
 * @return first valid index for p in y direction
 */
int StaggeredGrid::pJBegin() const
{
    return -1;
};

/**
 * get one after last valid index for p in y direction
 * @return one after last valid index for p in y direction
 */
int StaggeredGrid::pJEnd() const
{
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable p
 * @return size of FieldVariable p
 */
std::array<int, 2> StaggeredGrid::pSize() const
{
    return {pIEnd() - pIBegin(), pJEnd() - pJBegin()};
}

/**
 * get first valid Interior index for p in x direction
 * @return first valid Interior index for p in x direction
 */
int StaggeredGrid::pInteriorIBegin() const
{
    return pIBegin() + 1;
};

/**
 * get one after last valid Interior index for p in x direction
 * @return one after last valid Interior index for p in x direction
 */
int StaggeredGrid::pInteriorIEnd() const
{
    return pIEnd() - 1;
};

/**
 * get first valid Interior index for p in y direction
 * @return first valid Interior index for p in y direction
 */
int StaggeredGrid::pInteriorJBegin() const
{
    return pJBegin() + 1;
};

/**
 * get first valid Interior index for p in y direction
 * @return first valid Interior index for p in y direction
 */
int StaggeredGrid::pInteriorJEnd() const
{
    return pJEnd() - 1;
};

/**
 * get reference to field variable p
 * @return reference to field variable p
 */
const FieldVariable &StaggeredGrid::p() const
{
    return p_;
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
double StaggeredGrid::p(int i, int j) const
{
    #ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
    #endif
    return p_(i - pIBegin(), j - pJBegin());
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
double &StaggeredGrid::p(int i, int j)
{
    #ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
    #endif
    return p_(i - pIBegin(), j - pJBegin());
};

/**
 * velocity in x-direction u
 */

/**
 * get first valid index for u in x direction
 * @return first valid index for u in x direction
 */
int StaggeredGrid::uIBegin() const
{
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        return -1;
    }
    return -2;
};

/**
 * get one after last valid index for u in x direction
 * @return one after last valid index for u in x direction
 */
int StaggeredGrid::uIEnd() const
{
    if (partitioning_->ownPartitionContainsRightBoundary()) {
        return nCells_[0];
    }
    return nCells_[0] + 1;
};

/**
 * get first valid index for u in y direction
 * @return first valid index for u in y direction
 */
int StaggeredGrid::uJBegin() const
{
    return -1;
};

/**
 * get one after last valid index for u in y direction
 * @return one after last valid index for u in y direction
 */
int StaggeredGrid::uJEnd() const
{
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable u
 * @return size of FieldVariable u
 */
std::array<int, 2> StaggeredGrid::uSize() const
{
    return {uIEnd() - uIBegin(), uJEnd() - uJBegin()};
}

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::uInteriorIBegin() const
{

    return uIBegin() + 1;
};

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::uInteriorIEnd() const
{

    return uIEnd() - 1;

};

/**
 * first valid Interior index for u in y direction
 * @return first inner grid value of u in j direction
 */
int StaggeredGrid::uInteriorJBegin() const
{

    return uJBegin() + 1;

};

/**
 * get one after last valid Interior index for u in y direction
 * @return one after last valid Interior index for u in y direction
 */
int StaggeredGrid::uInteriorJEnd() const
{

    return uJEnd() - 1;

};

/**
 * get a reference to field variable u
 * @return reference to field variable u
 */
const FieldVariable &StaggeredGrid::u() const
{
    return u_;
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double StaggeredGrid::u(int i, int j) const
{
    #ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    #endif
    return u_(i - uIBegin(), j - uJBegin());
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double &StaggeredGrid::u(int i, int j)
{
    #ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    #endif
    return u_(i - uIBegin(), j - uJBegin());
};

/**
 * velocity in y-direction v
 */

/**
 * get first valid index for v in x direction
 * @return first valid index for v in x direction
 */
int StaggeredGrid::vIBegin() const
{
    return -1;
};

/**
 * get one after last valid index for v in x direction
 * @return one after last valid index for v in x direction
 */
int StaggeredGrid::vIEnd() const
{
    return nCells_[0] + 1;
};

/**
 * get first valid index for v in y direction
 * @return first valid index for v in y direction
 */
int StaggeredGrid::vJBegin() const
{
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        return -1;
    }
    return -2;
};

/**
 * get one after last valid index for v in y direction
 * @return one after last valid index for v in y direction
 */
int StaggeredGrid::vJEnd() const
{
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        return nCells_[1];
    }
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable v
 * @return size of FieldVariable v
 */
std::array<int, 2> StaggeredGrid::vSize() const
{
    return {vIEnd() - vIBegin(), vJEnd() - vJBegin()};
}

/**
 * get first valid Interior index for v in x direction
 * @return first valid Interior index for v in x direction
 */
int StaggeredGrid::vInteriorIBegin() const
{
    return vIBegin() + 1;

};


/**
 * get one after last valid Interior index for v in x direction
 * @return one after last valid Interior index for v in x direction
 */
int StaggeredGrid::vInteriorIEnd() const
{
    return vIEnd() - 1;

};

/**
 * get first valid Interior index for v in y direction
 * @return first valid Interior index for v in y direction
 */
int StaggeredGrid::vInteriorJBegin() const
{
    return vJBegin() + 1;

};

/**
 * get one after last valid Interior index for v in y direction
 * @return one after last valid Interior index for v in y direction
 */
int StaggeredGrid::vInteriorJEnd() const
{
    return vJEnd() - 1;

};

/**
 * get a reference to field variable v
 * @return a reference to field variable v
 */
const FieldVariable &StaggeredGrid::v() const
{
    return v_;
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double StaggeredGrid::v(int i, int j) const
{
    #ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    #endif
    return v_(i - vIBegin(), j - vJBegin());
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double &StaggeredGrid::v(int i, int j)
{
    #ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    #endif
    return v_(i - vIBegin(), j - vJBegin());
};

/**
 * right hand side rhs
 */

/**
 * get first valid index for rhs in x direction
 * @return first valid index for rhs in x direction
 */
int StaggeredGrid::rhsIBegin() const
{
    return -1;
};

/**
 * get one after last valid index for rhs in x direction
 * @return one after last valid index for rhs in x direction
 */
int StaggeredGrid::rhsIEnd() const
{
    return nCells_[0] + 1;
};


/**
 * get first valid index for rhs in y direction
 * @return first valid index for rhs in y direction
 */
int StaggeredGrid::rhsJBegin() const
{
    return -1;
};

/**
 * get one after last valid index for rhs in y direction
 * @return one after last valid index for rhs in y direction
 */
int StaggeredGrid::rhsJEnd() const
{
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable rhs
 * @return size of FieldVariable rhs
 */
std::array<int, 2> StaggeredGrid::rhsSize() const
{
    return {rhsIEnd() - rhsIBegin(), rhsJEnd() - rhsJBegin()};
}

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::rhsInteriorIBegin() const
{

    return rhsIBegin() + 1;
};

/**
 * get one after last valid Interior index for u in x direction
 * @return one after last valid Interior index for u in x direction
 */
int StaggeredGrid::rhsInteriorIEnd() const
{

    return rhsIEnd() - 1;

};

/**
 * get first valid Interior index for u in y direction
 * @return first valid Interior index for u in y direction
 */
int StaggeredGrid::rhsInteriorJBegin() const
{

    return rhsJBegin() + 1;

};

/**
 * get one after last valid Interior index for u in y direction
 * @return one after last valid Interior index for u in y direction
 */
int StaggeredGrid::rhsInteriorJEnd() const
{

    return rhsJEnd() - 1;

};

/**
 * access value of rhs in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of rhs in element (i,j)
 */
double StaggeredGrid::rhs(int i, int j) const
{
    #ifndef NDEBUG
    assert((rhsIBegin() <= i) && (i <= rhsIEnd()));
    assert((rhsJBegin() <= j) && (j <= rhsJEnd()));
    #endif
    return rhs_(i - rhsIBegin(), j - rhsJBegin());
};

/**
 * access value of rhs in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of rhs in element (i,j)
 */
double &StaggeredGrid::rhs(int i, int j)
{
    #ifndef NDEBUG
    assert((rhsIBegin() <= i) && (i <= rhsIEnd()));
    assert((rhsJBegin() <= j) && (j <= rhsJEnd()));
    #endif
    return rhs_(i - rhsIBegin(), j - rhsJBegin());
};

/**
 * get reference to field variable rhs
 * @return reference to field variable rhs
 */
const FieldVariable &StaggeredGrid::rhs() const
{
    return rhs_;
};

/**
 * preliminary velocity F
 */

/**
 * access value of F in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of F in element (i,j)
 */
double StaggeredGrid::f(int i, int j) const
{
    #ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    #endif
    return f_(i - uIBegin(), j - uJBegin());
};

/**
 * access value of F in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of F in element (i,j)
 */
double &StaggeredGrid::f(int i, int j)
{
    #ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
    #endif
    return f_(i - uIBegin(), j - uJBegin());
};


/**
 * get reference to field variable F
 * @return reference to field variable F
 */
const FieldVariable &StaggeredGrid::f() const
{
    return f_;
};

/**
 * preliminary velocity G
 */


/**
 * access value of G in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of G in element (i,j)
 */
double StaggeredGrid::g(int i, int j) const
{
    #ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    #endif
    return g_(i - vIBegin(), j - vJBegin());
};

/**
 * access value of G in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of G in element (i,j)
 */
double &StaggeredGrid::g(int i, int j)
{
    #ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
    #endif
    return g_(i - vIBegin(), j - vJBegin());
};

/**
 * get reference to field variable G
 * @return reference to field variable G
 */
const FieldVariable &StaggeredGrid::g() const
{
    return g_;
};

/**
 * temperature variable t
 */

/**
 * get first valid index for t in x direction
 * @return first valid index for t in x direction
 */
int StaggeredGrid::tIBegin() const
{
    return -1;
};

/**
 * get one after last valid index for t in x direction
 * @return one after last valid index for t in x direction
 */
int StaggeredGrid::tIEnd() const
{
    return nCells_[0] + 1;
};

/**
 * get first valid index for t in y direction
 * @return first valid index for t in y direction
 */
int StaggeredGrid::tJBegin() const
{
    return -1;
};

/**
 * get one after last valid index for t in y direction
 * @return one after last valid index for t in y direction
 */
int StaggeredGrid::tJEnd() const
{
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable t
 * @return size of FieldVariable t
 */
std::array<int, 2> StaggeredGrid::tSize() const
{
    return {tIEnd() - tIBegin(), tJEnd() - tJBegin()};
}

/**
 * get first valid Interior index for t in x direction
 * @return first valid Interior index for t in x direction
 */
int StaggeredGrid::tInteriorIBegin() const
{
    return tIBegin() + 1;
};

/**
 * get one after last valid Interior index for t in x direction
 * @return one after last valid Interior index for t in x direction
 */
int StaggeredGrid::tInteriorIEnd() const
{
    return tIEnd() - 1;
};

/**
 * get first valid Interior index for t in y direction
 * @return first valid Interior index for t in y direction
 */
int StaggeredGrid::tInteriorJBegin() const
{
    return tJBegin() + 1;
};

/**
 * get first valid Interior index for t in y direction
 * @return first valid Interior index for t in y direction
 */
int StaggeredGrid::tInteriorJEnd() const
{
    return tJEnd() - 1;
};

/**
 * get reference to field variable t
 * @return reference to field variable t
 */
const FieldVariable &StaggeredGrid::t() const
{
    return t_;
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double StaggeredGrid::t(int i, int j) const
{
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()));
    assert((tJBegin() <= j) && (j <= tJEnd()));
#endif
    return t_(i - tIBegin(), j - tJBegin());
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double &StaggeredGrid::t(int i, int j)
{
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()));
    assert((tJBegin() <= j) && (j <= tJEnd()));
#endif
    return t_(i - tIBegin(), j - tJBegin());
};