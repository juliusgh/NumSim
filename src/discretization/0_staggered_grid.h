#pragma once

#include <cassert>
#include <memory>
#include "storage/field_variable.h"
#include "storage/marker2d.h"
#include "partitioning/partitioning.h"
#include "settings.h"

/**
 * Implement staggered grid, providing a variety of parameters
 */
class StaggeredGrid {
public:
    /**
    * constructor
    * @param nCells: number of cells
    * @param meshWidth: cell width in all directions
    */
    StaggeredGrid(const std::shared_ptr<Partitioning>& partitioning,
                  std::array<double, 2> meshWidth,
                  std::shared_ptr<Settings> settings);

    /**
     * get mesh width in x-direction
     */
    double dx() const;

    /**
     * get mesh width in y-direction
     */
    double dy() const;

    /**
     * evaluate value of F in an element (i,j)
     * @param i: position in x direction in discretized grid
     * @param j: position in y direction in discretized grid
     */
    double f(int i, int j) const;

    /**
     *  evaluate value of F in an element (i,j)
     * @param i: position in x direction in discretized grid
     * @param j: position in y direction in discretized grid
     * @return
     */
    double &f(int i, int j);

    /**
     * get a reference to field variable F
     */
    const FieldVariable &f() const;

    /**
     * evaluate value of G in an element (i,j)
    */
    double g(int i, int j) const;

    /**
     * get a reference to field variable G
    */
    const FieldVariable &g() const;

    /**
     * evaluate value of G in an element (i,j)
    */
    double &g(int i, int j);

    /**
     * get the mesh width
    */
    const std::array<double, 2> meshWidth() const;

    /**
     * get number of cells
    */
    const std::array<int, 2> nCells() const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    MARKER marker(int i, int j) const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    MARKER &marker(int i, int j);

    /**
     * get reference to field variable p
    */
    const FieldVariable &p() const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double p(int i, int j) const;

    /**
     * evaluate field variable p in an element (i,j)
    */
    double &p(int i, int j);

    /**
     * first valid index for p in x direction
    */
    int pIBegin() const;

    /**
     * one after last valid index for p in x direction
    */
    int pIEnd() const;

    /**
     * first valid index for p in y direction
    */
    int pJBegin() const;

    /**
     * one after last valid index for p in y direction
    */
    int pJEnd() const;

    /**
     * get size of FieldVariable p
    */
    std::array<int, 2> pSize() const;

    /**
     * first valid Interior index for p in x direction
   */
    int pInteriorIBegin() const;

    /**
     * one after last valid Interior index for p in x direction
    */
    int pInteriorIEnd() const;

    /**
     * first valid Interior index for p in y direction
    */
    int pInteriorJBegin() const;

    /**
     * one after last valid Interior index for p in y direction
    */
    int pInteriorJEnd() const;

    /**
     * get a reference to field variable u
    */
    const FieldVariable &u() const;

    /**
     * access value of u in element (i,j)
    */
    double u(int i, int j) const;

    /**
     * access value of u in element (i,j)
    */
    double &u(int i, int j);

    /**
     * get a reference to field variable u
    */
    const FieldVariable &uLast() const;

    /**
     * access value of u in element (i,j)
    */
    double uLast(int i, int j) const;

    /**
     * access value of u in element (i,j)
    */
    double &uLast(int i, int j);

    /**
     * first valid index for u in x direction
    */
    int uIBegin() const;

    /**
     * one after last valid index for u in x direction
    */
    int uIEnd() const;

    /**
     * first valid index for u in y direction
    */
    int uJBegin() const;

    /**
     * one after last valid index for u in y direction
    */
    int uJEnd() const;

    /**
     * get size of FieldVariable u
    */
    std::array<int, 2> uSize() const;

    /**
     * first valid field index for u in x direction
    */
    int uInteriorIBegin() const;

    /**
     * one after last valid Interior index for u in x direction
    */
    int uInteriorIEnd() const;

    /**
     * first valid Interior index for u in y direction
    */
    int uInteriorJBegin() const;

    /**
     * one after last valid Interior index for u in y direction
    */
    int uInteriorJEnd() const;

    /**
     * get a reference to Interior variable v
    */
    const FieldVariable &v() const;

    /**
     * access value of v in element (i,j)
    */
    double v(int i, int j) const;

    /**
     * access value of v in element (i,j)
    */
    double &v(int i, int j);

    /**
     * get a reference to Interior variable v
    */
    const FieldVariable &vLast() const;

    /**
     * access value of v in element (i,j)
    */
    double vLast(int i, int j) const;

    /**
     * access value of v in element (i,j)
    */
    double &vLast(int i, int j);

    /**
     * first valid index for v in x direction
    */
    int vIBegin() const;

    /**
     * one after last valid index for v in x direction
    */
    int vIEnd() const;

    /**
     * first valid index for v in y direction
    */
    int vJBegin() const;

    /**
     * one after last valid index for v in y direction
    */
    int vJEnd() const;

    /**
     * get size of FieldVariable v
    */
    std::array<int, 2> vSize() const;

    /**
     * first valid Interior index for v in x direction
    */
    int vInteriorIBegin() const;

    /**
     * one after last valid Interior index for v in x direction
    */
    int vInteriorIEnd() const;

    /**
     * first valid Interior index for v in y direction
    */
    int vInteriorJBegin() const;

    /**
     * one after last valid Interior index for v in y direction
    */
    int vInteriorJEnd() const;

    /**
     * first valid index for rhs in x direction
    */
    int rhsIBegin() const;

    /**
     * one after last valid index for rhs in x direction
    */
    int rhsIEnd() const;

    /**
     * first valid index for rhs in y direction
    */
    int rhsJBegin() const;

    /**
     * one after last valid index for rhs in y direction
    */
    int rhsJEnd() const;

    /**
     * first valid Interior index for v in x direction
    */
    int rhsInteriorIBegin() const;

    /**
     * one after last valid Interior index for v in x direction
    */
    int rhsInteriorIEnd() const;

    /**
     * first valid Interior index for v in y direction
    */
    int rhsInteriorJBegin() const;

    /**
     * one after last valid Interior index for v in y direction
    */
    int rhsInteriorJEnd() const;

    /**
     * get size of FieldVariable rhs
    */
    std::array<int, 2> rhsSize() const;

    /**
     * access value of rhs in element (i,j)
    */
    double rhs(int i, int j) const;

    /**
     * access value of rhs in element (i,j)
    */
    double &rhs(int i, int j);

    /**
    * get a reference to field variable rhs
    */
    const FieldVariable &rhs() const;

    /**
    * get reference to field variable t
    */
    const FieldVariable &t() const;

    /**
     * evaluate field variable t in an element (i,j)
    */
    double t(int i, int j) const;

    /**
     * evaluate field variable t in an element (i,j)
    */
    double &t(int i, int j);

    /**
     * first valid index for t in x direction
    */
    int tIBegin() const;

    /**
     * one after last valid index for t in x direction
    */
    int tIEnd() const;

    /**
     * first valid index for t in y direction
    */
    int tJBegin() const;

    /**
     * one after last valid index for t in y direction
    */
    int tJEnd() const;

    /**
     * get size of FieldVariable t
    */
    std::array<int, 2> tSize() const;

    /**
     * first valid Interior index for t in x direction
   */
    int tInteriorIBegin() const;

    /**
     * one after last valid Interior index for t in x direction
    */
    int tInteriorIEnd() const;

    /**
     * first valid Interior index for t in y direction
    */
    int tInteriorJBegin() const;

    /**
     * one after last valid Interior index for t in y direction
    */
    int tInteriorJEnd() const;

    /**
    * get reference to field variable q
    */
    const FieldVariable &q() const;

    /**
     * evaluate field variable q in an element (i,j)
    */
    double q(int i, int j) const;

    /**
     * evaluate field variable q in an element (i,j)
    */
    double &q(int i, int j);

    /**
     * first valid index for q in x direction
    */
    int qIBegin() const;

    /**
     * one after last valid index for q in x direction
    */
    int qIEnd() const;

    /**
     * first valid index for q in y direction
    */
    int qJBegin() const;

    /**
     * one after last valid index for q in y direction
    */
    int qJEnd() const;

    /**
     * get size of FieldVariable q
    */
    std::array<int, 2> qSize() const;

    /**
     * first valid Interior index for q in x direction
   */
    int qInteriorIBegin() const;

    /**
     * one after last valid Interior index for q in x direction
    */
    int qInteriorIEnd() const;

    /**
     * first valid Interior index for q in y direction
    */
    int qInteriorJBegin() const;

    /**
     * one after last valid Interior index for q in y direction
    */
    int qInteriorJEnd() const;

    void applyBoundaryVelocities();

    void applyBoundaryPressure();

    void applyBoundaryTemperature();

    void setObstacleValues();


protected:
    const std::array<double, 2> meshWidth_;
    const std::array<int, 2> nCells_;
    std::shared_ptr<Partitioning> partitioning_;
    FieldVariable f_;
    FieldVariable g_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable u_;
    FieldVariable v_;
    FieldVariable uLast_;
    FieldVariable vLast_;
    FieldVariable t_;
    FieldVariable q_;
    Marker2D marker_;
    shared_ptr<Settings> settings_;
};