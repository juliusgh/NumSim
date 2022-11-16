#pragma once

#include "storage/fieldvariable.h"

/**
 * Implement staggered grid, providing a variety of parameters
 */
class StaggeredGrid
{
public:
    /**
    * constructor
    * @param nCells: number of cells
    * @param meshWidth: cell width in all directions
    */
    StaggeredGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth);
    
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
    double & f(int i, int j);

    /**
     * get a reference to field variable F
     */
    const FieldVariable & f() const;

    /**
     * evaluate value of G in an element (i,j)
    */
     double g(int i, int j) const;

    /**
     * get a reference to field variable G
    */
     const FieldVariable & g() const;

    /**
     * evaluate value of G in an element (i,j)
    */
     double & g(int i, int j);

    /**
     * get the mesh width
    */
     const std::array<double, 2> meshWidth() const;

    /**
     * get number of cells
    */
     const std::array<int, 2> nCells() const;

    /**
     * get reference to field variable p
    */
     const FieldVariable & p() const;

    /**
     * evaluate field variable p in an element (i,j)
    */
     double p(int i, int j) const;

    /**
     * evaluate field variable p in an element (i,j)
    */
     double & p(int i, int j);
    
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
     const FieldVariable & u() const;
    
    /**
     * access value of u in element (i,j)
    */
     double u(int i, int j) const;
    
    /**
     * access value of u in element (i,j)
    */
     double & u(int i, int j);
    
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
     const FieldVariable & v() const;
    
    /**
     * access value of v in element (i,j)
    */
     double v(int i, int j) const;
    
    /**
     * access value of v in element (i,j)
    */
     double & v(int i, int j);

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
     double & rhs(int i, int j);

    /**
    * get a reference to field variable u
    */
     const FieldVariable & rhs() const;

protected:
    const std::array<double, 2> meshWidth_;
    const std::array<int, 2> nCells_;
    FieldVariable f_;
    FieldVariable g_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable u_;
    FieldVariable v_;
};