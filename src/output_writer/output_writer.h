#pragma once

#include "discretization/1_discretization.h"

/** Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 */

class OutputWriter {
public:
    //! constructor
    OutputWriter(std::shared_ptr <Discretization> discretization);

    //! write current velocities to file, filename is output_<count>.vti
    virtual void writeFile(double currentTime) = 0;

protected:
    std::shared_ptr <Discretization> discretization_;
    int fileNo_ = 0;
};
