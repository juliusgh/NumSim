#pragma once

#include "discretization/1_discretization.h"

/**
 * calculate derivatives needed for pressure calculations using the donor cell approach
 */

class DonorCell : public Discretization {
public:

    /**
    * use the constructor of the base class
    * @param partitioning: encapsulate functionality corresponding to subdomain handling
    * @param meshWidth: width of a cell in both directions
    * @param alpha: donor cell weight parameter
    * @param gamma: donor cell weight parameter for temperature
    * @param settings: information about the settings received from parameter file
    */
    DonorCell(std::shared_ptr<Partitioning> partitioning,
              std::array<double, 2> meshWidth,
              double alpha,
              double gamma,
              Settings settings);


    /**
    * compute the 1st derivative ∂ u^2 / ∂x
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return donor cell derivative approximation of the derivative stated above
    */
    virtual double computeDu2Dx(int i, int j) const;

    /**
    * compute the 1st derivative ∂ v^2 / ∂y
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return donor cell derivative approximation of the derivative stated above
    */
    virtual double computeDv2Dy(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (uv) / ∂x
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return donor cell derivative approximation of the derivative stated above
    */
    virtual double computeDuvDx(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (uv) / ∂y
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return donor cell derivative approximation of the derivative stated above
    */
    virtual double computeDuvDy(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (ut) / ∂x
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return donor cell derivative approximation of the derivative stated above
    */
    virtual double computeDutDx(int i, int j) const;

    /**
    * compute the 1st derivative ∂ (vt) / ∂y
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return donor cell derivative approximation of the derivative stated above
    */
    virtual double computeDvtDy(int i, int j) const;

private:
    double alpha_;
    double gamma_;
};