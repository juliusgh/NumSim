#pragma once

#include "discretization/1_discretization.h"

class CentralDifferences : public Discretization
{
public:

  //! use the constructor of the base class
  using Discretization::Discretization;

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const;

  //! compute the 1st derivative ∂ v^2 / ∂y
  virtual double computeDv2Dy(int i, int j) const;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const;

  //! compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j) const;

  
};