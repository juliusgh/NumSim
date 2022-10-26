#pragma once

#include "discretization/1_discretization.h"

class DonorCell : public Discretization
{
public:

  //! use the constructor of the base class
  DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha);

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const;

  //! compute the 1st derivative ∂ v^2 / ∂x
  virtual double computeDv2Dy(int i, int j) const;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const;

//...
};