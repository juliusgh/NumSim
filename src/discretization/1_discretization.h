#pragma once

#include <array>
#include <memory>
#include "discretization/0_staggered_grid.h"
#include "partitioning/partitioning.h"

/**
 * Calculate derivatives needed for velocity calculations
 */

class Discretization :
        public StaggeredGrid {
public:


    /**
     * construct the object with given number of cells in x and y direction
     * @param nCells: number of inner cells
     * @param meshWidth: width of grid cell in one direction
     */
    Discretization(std::shared_ptr<Partitioning> partitioning, std::array<double, 2> meshWidth);

    /**
    * compute the 2nd derivative ∂^2 u / ∂x^2
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeD2uDx2(int i, int j) const;

    /**
    * compute the 2nd derivative ∂^2 u / ∂y^2
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeD2uDy2(int i, int j) const;

    /**
    * compute the 2nd derivative ∂^2 v / ∂x^2
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeD2vDx2(int i, int j) const;

    /**
  * compute the 2nd derivative ∂^2 v / ∂y^2
  * @param i: discretized position in x direcetion
  * @param j: discretized position in y direction
  * @return derivative approximation of the derivative stated above
  */
    virtual double computeD2vDy2(int i, int j) const;

    /**
  * compute the 1st derivative ∂p / ∂x
  * @param i: discretized position in x direcetion
  * @param j: discretized position in y direction
  * @return derivative approximation of the derivative stated above
  */
    virtual double computeDpDx(int i, int j) const;

    /**
   * compute the 1st derivative ∂p / ∂y
   * @param i: discretized position in x direcetion
   * @param j: discretized position in y direction
   * @return derivative approximation of the derivative stated above
   */
    virtual double computeDpDy(int i, int j) const;

    /**
    * compute the 1st derivative ∂ u^2 / ∂x
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeDu2Dx(int i, int j) const = 0;

    /**
    * compute the 1st derivative ∂ (uv) / ∂x
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeDuvDx(int i, int j) const = 0;

    /**
    * compute the 1st derivative ∂ (uv) / ∂y
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeDuvDy(int i, int j) const = 0;

    /**
    * compute the 1st derivative ∂ v^2 / ∂y
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeDv2Dy(int i, int j) const = 0;

    /**
    * compute the 2nd derivative ∂²t / ∂x²
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeD2tD2x(int i, int j) const;

    /**
    * compute the 2nd derivative ∂²t / ∂y²
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeD2tD2y(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (ut) / ∂x
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeDutDx(int i, int j) const = 0;

    /**
    * compute the 1st derivative ∂ (vt) / ∂y
    * @param i: discretized position in x direcetion
    * @param j: discretized position in y direction
    * @return derivative approximation of the derivative stated above
    */
    virtual double computeDvtDy(int i, int j) const = 0;

};