#pragma once

#include "discretization/1_discretization.h"
#include "partitioning/partitioning.h"

#include <memory>

/** Inteface class for writing simulation data output.
 */
class OutputWriter
{
public:
  //! constructor
  OutputWriter(std::shared_ptr<Discretization> discretization, const Partitioning &partitioning);

  //! write current velocities to file, filename is output_<count>.vti
  virtual void writeFile(double currentTime) = 0;

protected:

  std::shared_ptr<Discretization> discretization_;  //< a shared pointer to the discretization which contains all data that will be written to the file
  const Partitioning partitioning_;                 //< the partitioning object that knowns about the domain decomposition, only significant when executing in parallel
  int fileNo_;   //< a counter that increments for every file, this number is part of the file name of output files
};
