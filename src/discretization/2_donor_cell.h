#pragma once

#include "discretization/1_discretization.h"

/**
 *  calculate derivatives needed for pressure calculations using the donor cell approach
 */

class DonorCell : public Discretization
{
public:

  /**
   * use the constructor of the base class
   * @param nCells: number of inner cells
   * @param meshWidth: width of a cell in both directions
   * @param alpha: donor cell weight parameter
   */
  DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha);

  /**
   * compute the 1st derivative ∂ u^2 / ∂x
   * @param i: discretized position in x direcetion
   * @param j: discretiszed position in y direction
   * @return donor cell derivative approximation of the derivative stated above
   */
  virtual double computeDu2Dx(int i, int j) const;

  /**
  * compute the 1st derivative ∂ v^2 / ∂y
  * @param i: discretized position in x direcetion
  * @param j: discretiszed position in y direction
  * @return donor cell derivative approximation of the derivative stated above
  */
  virtual double computeDv2Dy(int i, int j) const;

  /**
  * compute the 1st derivative ∂ (uv) / ∂x
  * @param i: discretized position in x direcetion
  * @param j: discretiszed position in y direction
  * @return donor cell derivative approximation of the derivative stated above
  */
  virtual double computeDuvDx(int i, int j) const;

  /**
  * compute the 1st derivative ∂ (uv) / ∂y
  * @param i: discretized position in x direcetion
  * @param j: discretiszed position in y direction
  * @return donor cell derivative approximation of the derivative stated above
  */
  virtual double computeDuvDy(int i, int j) const;
  
private:
  double alpha_;
};