#pragma once

#include "discretization/1_discretization.h"

/**
 *  calculate derivatives needed for pressure calculations using the central differences approach
 */

class CentralDifferences : public Discretization {
public:
    CentralDifferences(std::shared_ptr<Partitioning> partitioning,
                       std::array<double, 2> meshWidth,
                       Settings settings);

    /**
    * compute the 1st derivative ∂ u^2 / ∂x
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return central differences derivative approximation of the derivative stated above
    */
    virtual double computeDu2Dx(int i, int j) const;

    /**
    * compute the 1st derivative ∂ v^2 / ∂y
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return central differences derivative approximation of the derivative stated above
    */
    virtual double computeDv2Dy(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (uv) / ∂x
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return central differences derivative approximation of the derivative stated above
    */
    virtual double computeDuvDx(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (uv) / ∂y
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return central differences derivative approximation of the derivative stated above
    */
    virtual double computeDuvDy(int i, int j) const;

    /**
     * compute the 1st derivative ∂ (ut) / ∂x
     * @param i: discretized position in x direction
     * @param j: discretized position in y direction
     * @return central differences derivative approximation of the derivative stated above
     */
    virtual double computeDutDx(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (vt) / ∂y
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return central differences derivative approximation of the derivative stated above
    */
    virtual double computeDvtDy(int i, int j) const;

};