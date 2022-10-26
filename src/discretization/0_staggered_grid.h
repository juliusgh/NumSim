#pragma once

#include "storage/fieldvariable.h"

class StaggeredGrid
{
public:
    //! constructor
    StaggeredGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth);

    //! get mesh width in x-direction
    double dx() const;
    
    //! get mesh width in y-direction
    double dy() const;
    
    //! evaluate value of F in an element (i,j)
    double f(int i, int j) const;
    
    //! evaluate value of G in an element (i,j)
    double g(int i, int j) const;

    //! get the mesh width
    const std::array<double, 2> meshWidth() const;

    //! get number of cells
    const std::array<int, 2> nCells() const;

    //! get reference to field variable p
    const FieldVariable & p() const;

    //! evaluate field variable p in an element (i,j)
    double p(int i, int j) const;

    //! evaluate field variable p in an element (i,j)
    double & p(int i, int j);
    
    //! first valid index for p in x direction
    int pIBegin() const;
    
    //! one after last valid index for p in x direction
    int pIEnd() const;
    
    //! first valid index for p in y direction
    int pJBegin() const;
    
    //! one after last valid index for p in y direction
    int pJEnd() const;
    
    //! access value of rhs in element (i,j)
    double & rhs(int i, int j);
    
    //! get a reference to field variable u
    const FieldVariable & u() const;
    
    //! access value of u in element (i,j)
    double u(int i, int j) const;
    
    //! access value of u in element (i,j)
    double & u(int i, int j);
    
    //! first valid index for u in x direction
    int uIBegin() const;
    
    //! one after last valid index for u in x direction
    int uIEnd() const;
    
    //! first valid index for u in y direction
    int uJBegin() const;
    
    //! one after last valid index for u in y direction
    int uJEnd() const;
    
    //! get a reference to field variable v
    const FieldVariable & v() const;
    
    //! access value of v in element (i,j)
    double v(int i, int j) const;
    
    //! access value of v in element (i,j)
    double & v(int i, int j);
    
    //! first valid index for v in x direction
    int vIBegin() const;
    
    //! one after last valid index for v in x direction
    int vIEnd() const;
    
    //! first valid index for v in y direction
    int vJBegin() const;
    
    //! one after last valid index for v in y direction
    int vJEnd() const;

protected:
    FieldVariable f_;
    FieldVariable g_;
    const std::array<double, 2> meshWidth_;
    const std::array<int, 2> nCells_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable u_;
    FieldVariable v_;
};