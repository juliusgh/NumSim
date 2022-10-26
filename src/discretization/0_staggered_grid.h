#pragma once

#include "storage/fieldvariable.h"

class StaggeredGrid
{
public:
    //! constructor
    StaggeredGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth);

    //! get mesh width in x-direction
    double StaggeredGrid::dx() const;
    
    //! get mesh width in y-direction
    double StaggeredGrid::dy() const;
    
    //! evaluate value of F in an element (i,j)
    double StaggeredGrid::f(int i, int j) const;
    
    //! evaluate value of G in an element (i,j)
    double StaggeredGrid::g(int i, int j) const;

    //! get the mesh width
    const std::array<double, 2> StaggeredGrid::meshWidth() const;

    //! get number of cells
    const std::array<int, 2> StaggeredGrid::nCells() const;

    //! get reference to field variable p
    const FieldVariable & StaggeredGrid::p() const;

    //! evaluate field variable p in an element (i,j)
    double StaggeredGrid::p(int i, int j) const;

    //! evaluate field variable p in an element (i,j)
    double & StaggeredGrid::p(int i, int j);
    
    
protected:

};