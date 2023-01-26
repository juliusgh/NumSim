#include <iostream>
#include "discretization/0_staggered_grid.h"

/**
 * Implement staggered grid, providing a variety of parameters
 * @param nCells: number of cells
 * @param meshWidth: cell width in all directions
 */
StaggeredGrid::StaggeredGrid(std::shared_ptr<Partitioning> partitioning,
                             std::array<double, 2> meshWidth,
                             std::shared_ptr<Settings> settings) :
        partitioning_(partitioning),
        settings_(settings),
        nCells_(partitioning->nCellsLocal()),
        meshWidth_(meshWidth),
        marker_(pSize()),
        f_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
        g_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
        p_(pSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        rhs_(rhsSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        u_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
        v_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
        uLast_(uSize(), {meshWidth[0], meshWidth[1] / 2.}, meshWidth),
        vLast_(vSize(), {meshWidth[0] / 2., meshWidth[1]}, meshWidth),
        t_(tSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
        q_(tSize(), {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth) {
    // set markers
    // TODO: read markers from file
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            marker(i, j) = MARKER::FLUID;
        }
    }
    for (int i = pIBegin(); i < pIEnd(); i++) {
        // bottom
        marker(i, pJBegin()) = MARKER::NOSLIP;
        // top
        marker(i, pJEnd() - 1) = MARKER::INFLOW;
    }
    for (int j = pJBegin(); j < pJEnd(); j++) {
        // left
        marker(pIBegin(), j) = MARKER::NOSLIP;
        // right
        marker(pIEnd() - 1, j) = MARKER::NOSLIP;
    }
    std::cout << "finished setting markers" << std::endl;
};

/**
 * mesh information
 */

/**
 * get the mesh width
 * @return mesh width
 */
const std::array<double, 2> StaggeredGrid::meshWidth() const {
    return meshWidth_;
};

/**
 * get number of cells
 * @return number of cells
 */
const std::array<int, 2> StaggeredGrid::nCells() const {
    return nCells_;
};

/**
 * get mesh width in x-direction
 * @return mesh width in x-direction
 */
double StaggeredGrid::dx() const {
    return meshWidth_[0];
};

/**
 * get mesh width in y-direction
 * @return mesh width in y-direction
 */
double StaggeredGrid::dy() const {
    return meshWidth_[1];
};

/**
 * pressure variable p
 */

/**
 * get first valid index for p in x direction
 * @return first valid index for p in x direction
 */
int StaggeredGrid::pIBegin() const {
    return -1;
};

/**
 * get one after last valid index for p in x direction
 * @return one after last valid index for p in x direction
 */
int StaggeredGrid::pIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for p in y direction
 * @return first valid index for p in y direction
 */
int StaggeredGrid::pJBegin() const {
    return -1;
};

/**
 * get one after last valid index for p in y direction
 * @return one after last valid index for p in y direction
 */
int StaggeredGrid::pJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable p
 * @return size of FieldVariable p
 */
std::array<int, 2> StaggeredGrid::pSize() const {
    return {pIEnd() - pIBegin(), pJEnd() - pJBegin()};
}

/**
 * get first valid Interior index for p in x direction
 * @return first valid Interior index for p in x direction
 */
int StaggeredGrid::pInteriorIBegin() const {
    return pIBegin() + 1;
};

/**
 * get one after last valid Interior index for p in x direction
 * @return one after last valid Interior index for p in x direction
 */
int StaggeredGrid::pInteriorIEnd() const {
    return pIEnd() - 1;
};

/**
 * get first valid Interior index for p in y direction
 * @return first valid Interior index for p in y direction
 */
int StaggeredGrid::pInteriorJBegin() const {
    return pJBegin() + 1;
};

/**
 * get first valid Interior index for p in y direction
 * @return first valid Interior index for p in y direction
 */
int StaggeredGrid::pInteriorJEnd() const {
    return pJEnd() - 1;
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
StaggeredGrid::MARKER StaggeredGrid::marker(int i, int j) const {
#ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
#endif
    return (MARKER)marker_(i - pIBegin(), j - pJBegin());
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
StaggeredGrid::MARKER &StaggeredGrid::marker(int i, int j) {
#ifndef NDEBUG
    assert((pIBegin() <= i) && (i <= pIEnd()));
    assert((pJBegin() <= j) && (j <= pJEnd()));
#endif
    return (MARKER&)marker_(i - pIBegin(), j - pJBegin());
};

/**
 * get reference to field variable p
 * @return reference to field variable p
 */
const FieldVariable &StaggeredGrid::p() const {
    return p_;
};

/**
 * evaluate field variable p in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable p in an element (i,j)
 */
double StaggeredGrid::p(int i, int j) const {
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
double &StaggeredGrid::p(int i, int j) {
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
int StaggeredGrid::uIBegin() const {
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        return -1;
    }
    return -2;
};

/**
 * get one after last valid index for u in x direction
 * @return one after last valid index for u in x direction
 */
int StaggeredGrid::uIEnd() const {
    if (partitioning_->ownPartitionContainsRightBoundary()) {
        return nCells_[0];
    }
    return nCells_[0] + 1;
};

/**
 * get first valid index for u in y direction
 * @return first valid index for u in y direction
 */
int StaggeredGrid::uJBegin() const {
    return -1;
};

/**
 * get one after last valid index for u in y direction
 * @return one after last valid index for u in y direction
 */
int StaggeredGrid::uJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable u
 * @return size of FieldVariable u
 */
std::array<int, 2> StaggeredGrid::uSize() const {
    return {uIEnd() - uIBegin(), uJEnd() - uJBegin()};
}

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::uInteriorIBegin() const {

    return uIBegin() + 1;
};

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::uInteriorIEnd() const {

    return uIEnd() - 1;

};

/**
 * first valid Interior index for u in y direction
 * @return first inner grid value of u in j direction
 */
int StaggeredGrid::uInteriorJBegin() const {

    return uJBegin() + 1;

};

/**
 * get one after last valid Interior index for u in y direction
 * @return one after last valid Interior index for u in y direction
 */
int StaggeredGrid::uInteriorJEnd() const {

    return uJEnd() - 1;

};

/**
 * get a reference to field variable u
 * @return reference to field variable u
 */
const FieldVariable &StaggeredGrid::u() const {
    return u_;
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double StaggeredGrid::u(int i, int j) const {
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
double &StaggeredGrid::u(int i, int j) {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return u_(i - uIBegin(), j - uJBegin());
};

/**
 * get a reference to field variable u
 * @return reference to field variable u
 */
const FieldVariable &StaggeredGrid::uLast() const {
    return uLast_;
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double StaggeredGrid::uLast(int i, int j) const {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return uLast_(i - uIBegin(), j - uJBegin());
};

/**
 * access value of u in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of u in element (i,j)
 */
double &StaggeredGrid::uLast(int i, int j) {
#ifndef NDEBUG
    assert((uIBegin() <= i) && (i <= uIEnd()));
    assert((uJBegin() <= j) && (j <= uJEnd()));
#endif
    return uLast_(i - uIBegin(), j - uJBegin());
};

/**
 * velocity in y-direction v
 */

/**
 * get first valid index for v in x direction
 * @return first valid index for v in x direction
 */
int StaggeredGrid::vIBegin() const {
    return -1;
};

/**
 * get one after last valid index for v in x direction
 * @return one after last valid index for v in x direction
 */
int StaggeredGrid::vIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for v in y direction
 * @return first valid index for v in y direction
 */
int StaggeredGrid::vJBegin() const {
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        return -1;
    }
    return -2;
};

/**
 * get one after last valid index for v in y direction
 * @return one after last valid index for v in y direction
 */
int StaggeredGrid::vJEnd() const {
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        return nCells_[1];
    }
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable v
 * @return size of FieldVariable v
 */
std::array<int, 2> StaggeredGrid::vSize() const {
    return {vIEnd() - vIBegin(), vJEnd() - vJBegin()};
}

/**
 * get first valid Interior index for v in x direction
 * @return first valid Interior index for v in x direction
 */
int StaggeredGrid::vInteriorIBegin() const {
    return vIBegin() + 1;

};


/**
 * get one after last valid Interior index for v in x direction
 * @return one after last valid Interior index for v in x direction
 */
int StaggeredGrid::vInteriorIEnd() const {
    return vIEnd() - 1;

};

/**
 * get first valid Interior index for v in y direction
 * @return first valid Interior index for v in y direction
 */
int StaggeredGrid::vInteriorJBegin() const {
    return vJBegin() + 1;

};

/**
 * get one after last valid Interior index for v in y direction
 * @return one after last valid Interior index for v in y direction
 */
int StaggeredGrid::vInteriorJEnd() const {
    return vJEnd() - 1;

};

/**
 * get a reference to field variable v
 * @return a reference to field variable v
 */
const FieldVariable &StaggeredGrid::v() const {
    return v_;
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double StaggeredGrid::v(int i, int j) const {
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
double &StaggeredGrid::v(int i, int j) {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return v_(i - vIBegin(), j - vJBegin());
};

/**
 * get a reference to field variable v
 * @return a reference to field variable v
 */
const FieldVariable &StaggeredGrid::vLast() const {
    return vLast_;
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double StaggeredGrid::vLast(int i, int j) const {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return vLast_(i - vIBegin(), j - vJBegin());
};

/**
 * access value of v in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of v in element (i,j)
 */
double &StaggeredGrid::vLast(int i, int j) {
#ifndef NDEBUG
    assert((vIBegin() <= i) && (i <= vIEnd()));
    assert((vJBegin() <= j) && (j <= vJEnd()));
#endif
    return vLast_(i - vIBegin(), j - vJBegin());
};

/**
 * right hand side rhs
 */

/**
 * get first valid index for rhs in x direction
 * @return first valid index for rhs in x direction
 */
int StaggeredGrid::rhsIBegin() const {
    return -1;
};

/**
 * get one after last valid index for rhs in x direction
 * @return one after last valid index for rhs in x direction
 */
int StaggeredGrid::rhsIEnd() const {
    return nCells_[0] + 1;
};


/**
 * get first valid index for rhs in y direction
 * @return first valid index for rhs in y direction
 */
int StaggeredGrid::rhsJBegin() const {
    return -1;
};

/**
 * get one after last valid index for rhs in y direction
 * @return one after last valid index for rhs in y direction
 */
int StaggeredGrid::rhsJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable rhs
 * @return size of FieldVariable rhs
 */
std::array<int, 2> StaggeredGrid::rhsSize() const {
    return {rhsIEnd() - rhsIBegin(), rhsJEnd() - rhsJBegin()};
}

/**
 * get first valid field index for u in x direction
 * @return first valid field index for u in x direction
 */
int StaggeredGrid::rhsInteriorIBegin() const {

    return rhsIBegin() + 1;
};

/**
 * get one after last valid Interior index for u in x direction
 * @return one after last valid Interior index for u in x direction
 */
int StaggeredGrid::rhsInteriorIEnd() const {

    return rhsIEnd() - 1;

};

/**
 * get first valid Interior index for u in y direction
 * @return first valid Interior index for u in y direction
 */
int StaggeredGrid::rhsInteriorJBegin() const {

    return rhsJBegin() + 1;

};

/**
 * get one after last valid Interior index for u in y direction
 * @return one after last valid Interior index for u in y direction
 */
int StaggeredGrid::rhsInteriorJEnd() const {

    return rhsJEnd() - 1;

};

/**
 * access value of rhs in element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return value of rhs in element (i,j)
 */
double StaggeredGrid::rhs(int i, int j) const {
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
double &StaggeredGrid::rhs(int i, int j) {
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
const FieldVariable &StaggeredGrid::rhs() const {
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
double StaggeredGrid::f(int i, int j) const {
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
double &StaggeredGrid::f(int i, int j) {
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
const FieldVariable &StaggeredGrid::f() const {
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
double StaggeredGrid::g(int i, int j) const {
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
double &StaggeredGrid::g(int i, int j) {
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
const FieldVariable &StaggeredGrid::g() const {
    return g_;
};

/**
 * temperature variable t
 */

/**
 * get first valid index for t in x direction
 * @return first valid index for t in x direction
 */
int StaggeredGrid::tIBegin() const {
    return -1;
};

/**
 * get one after last valid index for t in x direction
 * @return one after last valid index for t in x direction
 */
int StaggeredGrid::tIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for t in y direction
 * @return first valid index for t in y direction
 */
int StaggeredGrid::tJBegin() const {
    return -1;
};

/**
 * get one after last valid index for t in y direction
 * @return one after last valid index for t in y direction
 */
int StaggeredGrid::tJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable t
 * @return size of FieldVariable t
 */
std::array<int, 2> StaggeredGrid::tSize() const {
    return {tIEnd() - tIBegin(), tJEnd() - tJBegin()};
}

/**
 * get first valid Interior index for t in x direction
 * @return first valid Interior index for t in x direction
 */
int StaggeredGrid::tInteriorIBegin() const {
    return tIBegin() + 1;
};

/**
 * get one after last valid Interior index for t in x direction
 * @return one after last valid Interior index for t in x direction
 */
int StaggeredGrid::tInteriorIEnd() const {
    return tIEnd() - 1;
};

/**
 * get first valid Interior index for t in y direction
 * @return first valid Interior index for t in y direction
 */
int StaggeredGrid::tInteriorJBegin() const {
    return tJBegin() + 1;
};

/**
 * get first valid Interior index for t in y direction
 * @return first valid Interior index for t in y direction
 */
int StaggeredGrid::tInteriorJEnd() const {
    return tJEnd() - 1;
};

/**
 * get reference to field variable t
 * @return reference to field variable t
 */
const FieldVariable &StaggeredGrid::t() const {
    return t_;
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double StaggeredGrid::t(int i, int j) const {
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()) && "temperature i failed in const");
    assert((tJBegin() <= j) && (j <= tJEnd()) && "temperature j failed in const");
#endif
    return t_(i - tIBegin(), j - tJBegin());
};

/**
 * evaluate field variable t in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable t in an element (i,j)
 */
double &StaggeredGrid::t(int i, int j) {
#ifndef NDEBUG
    assert((tIBegin() <= i) && (i <= tIEnd()) && "temperature i failed");
    assert((tJBegin() <= j) && (j <= tJEnd()) && "temperature j failed");
#endif
    return t_(i - tIBegin(), j - tJBegin());
};

/**
 * sources and sinks for the heat equation q
 */

/**
 * get first valid index for q in x direction
 * @return first valid index for q in x direction
 */
int StaggeredGrid::qIBegin() const {
    return -1;
};

/**
 * get one after last valid index for q in x direction
 * @return one after last valid index for q in x direction
 */
int StaggeredGrid::qIEnd() const {
    return nCells_[0] + 1;
};

/**
 * get first valid index for q in y direction
 * @return first valid index for q in y direction
 */
int StaggeredGrid::qJBegin() const {
    return -1;
};

/**
 * get one after last valid index for q in y direction
 * @return one after last valid index for q in y direction
 */
int StaggeredGrid::qJEnd() const {
    return nCells_[1] + 1;
};

/**
 * get size of FieldVariable t
 * @return size of FieldVariable t
 */
std::array<int, 2> StaggeredGrid::qSize() const {
    return {qIEnd() - qIBegin(), qJEnd() - qJBegin()};
}

/**
 * get first valid Interior index for q in x direction
 * @return first valid Interior index for q in x direction
 */
int StaggeredGrid::qInteriorIBegin() const {
    return qIBegin() + 1;
};

/**
 * get one after last valid Interior index for q in x direction
 * @return one after last valid Interior index for q in x direction
 */
int StaggeredGrid::qInteriorIEnd() const {
    return qIEnd() - 1;
};

/**
 * get first valid Interior index for q in y direction
 * @return first valid Interior index for q in y direction
 */
int StaggeredGrid::qInteriorJBegin() const {
    return qJBegin() + 1;
};

/**
 * get first valid Interior index for q in y direction
 * @return first valid Interior index for q in y direction
 */
int StaggeredGrid::qInteriorJEnd() const {
    return qJEnd() - 1;
};

/**
 * get reference to field variable t
 * @return reference to field variable t
 */
const FieldVariable &StaggeredGrid::q() const {
    return q_;
};

/**
 * evaluate field variable q in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable q in an element (i,j)
 */
double StaggeredGrid::q(int i, int j) const {
#ifndef NDEBUG
    assert((qIBegin() <= i) && (i <= qIEnd()) && "Q i failed in const");
    assert((qJBegin() <= j) && (j <= qJEnd()) && "Q j failed in const");
#endif
    return q_(i - qIBegin(), j - qJBegin());
};

/**
 * evaluate field variable q in an element (i,j)
 * @param i: position in x direction in discretized grid
 * @param j: position in y direction in discretized grid
 * @return field variable q in an element (i,j)
 */
double &StaggeredGrid::q(int i, int j) {
#ifndef NDEBUG
    assert((qIBegin() <= i) && (i <= qIEnd()) && "Q i failed");
    assert((qJBegin() <= j) && (j <= qJEnd()) && "Q j failed");
#endif
    return q_(i - qIBegin(), j - qJBegin());
};

void StaggeredGrid::applyBoundaryVelocities() {
    for (int i = pInteriorIBegin(); i < pInteriorIEnd(); i++) {
        for (int j = pInteriorJBegin(); j < pInteriorJEnd(); j++) {
            switch (marker(i, j)) {
                case FLUID:
                case FREE:
                    break;
                case OBSTACLE:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    break;
                case OBSTACLE_LEFT:
                    f(i, j) = u(i, j) = 0.0;
                    f(i - 1, j) = u(i - 1, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i - 1, j);
                    break;
                case OBSTACLE_RIGHT:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i + 1, j);
                    break;
                case OBSTACLE_TOP:
                    f(i, j) = u(i, j) = -u(i, j + 1);
                    g(i, j) = v(i, j) = 0.0;
                    break;
                case OBSTACLE_BOTTOM:
                    f(i, j) = u(i, j) = -u(i, j - 1);
                    g(i, j) = v(i, j) = 0.0;
                    g(i, j - 1) = v(i, j - 1) = 0.0;
                    break;
                case OBSTACLE_LEFT_TOP:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    f(i - 1, j) = u(i - 1, j) = 0.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = 0.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    f(i, j) = u(i, j) = -u(i,j - 1);
                    f(i - 1, j) = u(i - 1, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i - 1, j);
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    f(i, j) = u(i, j) = 0.0;
                    g(i, j) = v(i, j) = -v(i + 1, j);
                    g(i, j - 1) = v(i, j - 1) = 0.0;
                    break;
                default:
                    break;
            }
        }
    }

    // set boundary values for u and v at bottom and top side (lower priority)
    for (int i = pIBegin(); i < pIEnd(); i++) {
        // set boundary values at bottom side
        const int uOffs = uIBegin() - pIBegin();
        const int vOffs = vIBegin() - pIBegin();
        switch (marker(i, pJBegin())) {
            case NOSLIP:
                if (i < pIEnd() - 1) {
                    f(i + uOffs, uJBegin()) = u(i + uOffs, uJBegin()) = -u(i + uOffs, uInteriorJBegin());
                }
                g(i + vOffs, vJBegin()) = v(i + vOffs, vJBegin()) = 0.0;
                break;
            case INFLOW:
                if (i < pIEnd() - 1) {
                    f(i + uOffs, uJBegin()) =
                    u(i + uOffs, uJBegin()) = 2.0 * settings_->dirichletBcBottom[0]
                                              - u(i + uOffs, uInteriorJBegin());
                }
                g(i + vOffs, vJBegin()) =
                v(i + vOffs, vJBegin()) = settings_->dirichletBcBottom[1];
                break;
            case OUTFLOW:
                if (i < pIEnd() - 1) {
                    f(i + uOffs, uJBegin()) =
                    u(i + uOffs, uJBegin()) = u(i + uOffs, uInteriorJBegin());
                }
                g(i + vOffs, vJBegin()) =
                v(i + vOffs, vJBegin()) = v(i + vOffs, vInteriorJBegin());
                break;
            default:
                break;
        }

        // set boundary values for u at top side
        switch (marker(i, uJEnd() - 1)) {
            case NOSLIP:
                if (i < pIEnd() - 1) {
                    u(i, uJEnd() - 1) = -u(i, uInteriorJEnd() - 1);
                }
                v(i, vJEnd() - 1) = 0.0;
                break;
            case INFLOW:
                if (i < pIEnd() - 1) {
                    u(i, uJEnd() - 1) = 2.0 * settings_->dirichletBcTop[0] - u(i, uInteriorJEnd() - 1);
                }
                v(i, vJEnd() - 1) = settings_->dirichletBcTop[1];
                break;
            case OUTFLOW:
                if (i < pIEnd() - 1) {
                    u(i, uJEnd() - 1) = u(i, uInteriorJEnd() - 1);
                }
                v(i, vJEnd() - 1) = v(i, vInteriorJEnd() - 1);
                break;
            default:
                break;
        }
    }

    // set boundary values for u and v at left and right side (higher priority)
    for (int j = pJBegin(); j < pJEnd(); j++) {
        // set boundary values for u at left side
        switch (marker(pIBegin(), j)) {
            case NOSLIP:
                u(uIBegin(), j) = 0.0;
                if (j < pJEnd() - 1) {
                    v(vIBegin(), j) = -v(vInteriorIBegin(), j);
                }
                break;
            case INFLOW:
                u(uIBegin(), j) = settings_->dirichletBcLeft[0];
                if (j < pJEnd() - 1) {
                    v(vIBegin(), j) = 2.0 * settings_->dirichletBcLeft[1]
                                      - v(vInteriorIBegin(), j);
                }
                break;
            case OUTFLOW:
                u(uIBegin(), j) = u(uInteriorIBegin(), j);
                if (j < pJEnd() - 1) {
                    v(vIBegin(), j) = v(vInteriorIBegin(), j);
                }
                break;
            default:
                break;
        }
        // set boundary values for u at right side
        switch (marker(uIEnd() - 1, j)) {
            case NOSLIP:
                u(uIEnd() - 1, j) = 0.0;
                if (j < pJEnd() - 1) {
                    v(vIEnd() - 1, j) = -v(vInteriorIEnd() - 1, j);
                }
                break;
            case INFLOW:
                u(uIEnd() - 1, j) = settings_->dirichletBcRight[0];
                if (j < pJEnd() - 1) {
                    v(vIEnd() - 1, j) = settings_->dirichletBcRight[1]
                                        - u(vInteriorIEnd() - 1, j);
                }
                break;
            case OUTFLOW:
                u(uIEnd() - 1, j) = u(uInteriorIEnd() - 1, j);
                if (j < pJEnd() - 1) {
                    v(vIEnd() - 1, j) = v(vInteriorIEnd() - 1, j);
                }
                break;
            default:
                break;
        }
    }
};

void StaggeredGrid::applyBoundaryPressure() {
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            switch (marker(i, j)) {
                case FLUID:
                case FREE:
                    break;
                case OBSTACLE:
                    p(i, j) = p(i - 1, j);
                    break;
                case OBSTACLE_LEFT:
                    p(i, j) = p(i - 1, j);
                    break;
                case OBSTACLE_RIGHT:
                    p(i, j) = p(i + 1, j);
                    break;
                case OBSTACLE_TOP:
                    p(i, j) = p(i,j + 1);
                    break;
                case OBSTACLE_BOTTOM:
                    p(i, j) = p(i, j - 1);
                    break;
                case OBSTACLE_LEFT_TOP:
                    p(i, j) = (p(i - 1, j) + p(i,j + 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    p(i, j) = (p(i + 1, j) + p(i,j + 1)) / 2.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    p(i, j) = (p(i - 1, j) + p(i,j - 1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    p(i, j) = (p(i + 1, j) + p(i,j - 1)) / 2.0;
                    break;
                default:
                    break;
            }
        }
    }

    // set boundary values for p at bottom and top side (lower priority)
    for (int i = pIBegin(); i < pIEnd(); i++) {
        // set boundary values at bottom side
        switch (marker(i, pJBegin())) {
            case INFLOW:
            case NOSLIP:
                p(i, pJBegin()) = p(i, pInteriorJBegin());
                break;
            case OUTFLOW:
                p(i, uJBegin()) = -p(i, uInteriorJBegin());
                break;
            default:
                break;
        }
        // set boundary values for p at top side
        switch (marker(i, pJEnd() - 1)) {
            case NOSLIP:
            case INFLOW:
                p(i, pJEnd() - 1) = p(i, pInteriorJEnd() - 1);
                break;
            case OUTFLOW:
                p(i, uJEnd() - 1) = -p(i, uInteriorJEnd() - 1);
                break;
            default:
                break;
        }
    }

    // set boundary values for p at left and right side (higher priority)
    for (int j = pJBegin(); j < pJEnd(); j++) {
        // set boundary values for p at left side
        switch (marker(pIBegin(), j)) {
            case NOSLIP:
            case INFLOW:
                p(pIBegin(), j) = p(pInteriorIBegin(), j);
                break;
            case OUTFLOW:
                p(pIBegin(), j) = -p(uInteriorIBegin(), j);
                break;
            default:
                break;
        }
        // set boundary values for u at right side
        switch (marker(uIEnd() - 1, j)) {
            case NOSLIP:
            case INFLOW:
                p(pIEnd() - 1, j) = p(pInteriorIEnd() - 1, j);
                break;
            case OUTFLOW:
                p(pIEnd() - 1, j) = -p(uInteriorIEnd() - 1, j);
                break;
            default:
                break;
        }
    }
};

void StaggeredGrid::applyBoundaryTemperature() {
    for (int i = pIBegin(); i < pIEnd(); i++) {
        for (int j = pJBegin(); j < pJEnd(); j++) {
            switch (marker(i, j)) {
                case FLUID:
                    break;
                case FREE:
                    break;
                case OBSTACLE:
                    break;
                case OBSTACLE_LEFT:
                    t(i,j) = t(i-1,j);
                    break;
                case OBSTACLE_RIGHT:
                    t(i,j) = t(i+1,j);
                    break;
                case OBSTACLE_TOP:
                    t(i,j) = t(i,j+1);
                    break;
                case OBSTACLE_BOTTOM:
                    t(i,j) = t(i,j-1);
                    break;
                case OBSTACLE_LEFT_TOP:
                    t(i,j) = (t(i-1,j) + t(i,j+1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_TOP:
                    t(i,j) = (t(i+1,j) + t(i,j+1)) / 2.0;
                    break;
                case OBSTACLE_LEFT_BOTTOM:
                    t(i,j) = (t(i-1,j) + t(i,j-1)) / 2.0;
                    break;
                case OBSTACLE_RIGHT_BOTTOM:
                    t(i,j) = (t(i+1,j) + t(i,j-1)) / 2.0;
                    break;
            }
        }
    }

    // set boundary values for t at bottom and top side (lower priotity)
    for (int i = tIBegin(); i < tIEnd(); i++) {
        if (settings_->setFixedTempBottom) {
            t(i, tJBegin()) = 2.0 * settings_->tempBcBottom
                               - t(i, tInteriorJBegin());
            //std::cout << "t(i, discretization_->tJBegin()) = " << discretization_->t(i, discretization_->tJBegin()) << std::endl;
        } else {
            t(i, tJBegin()) = t(i, tInteriorJBegin())
                               - dy() * settings_->tempBcBottom;
        }
        if (settings_->setFixedTempTop) {
            t(i, tJEnd() - 1) = 2.0 * settings_->tempBcTop
                                 - t(i, tInteriorJEnd() - 1);
        } else {
            t(i, tJEnd() - 1) = t(i, tInteriorJEnd() - 1)
                                 - dy() * settings_->tempBcTop;
        }
    }

    // set boundary values for t at left and right side (higher priority)
    for (int j = tJBegin(); j < tJEnd(); j++) {
        if (settings_->setFixedTempLeft) {
            t(tIBegin(), j) = 2.0 * settings_->tempBcLeft
                               - t(tInteriorIBegin(), j);
        } else {
            t(tIBegin(), j) = t(tInteriorIBegin(), j)
                                 - dx() * settings_->tempBcLeft;
        }
        if (settings_->setFixedTempRight) {
            t(tIEnd() - 1, j) = 2.0 * settings_->tempBcRight
                                 - t(tInteriorIEnd() - 1, j);
        } else {
            t(tIEnd() - 1, j) = t(tInteriorIEnd() - 1, j)
                                 - dx() * settings_->tempBcRight;
        }
    }
};