#pragma once

#include "output_writer/output_writer_text.h"
#include "partitioning/partitioning.h"

/** Write *.txt files that are useful for debugging.
 *  All values are written to the file as they are stored in the field variables,
 *  no interpolation takes place.
 */
class OutputWriterTextParallel : 
  public OutputWriterText
{
public:
    //! constructor
    OutputWriterTextParallel(std::shared_ptr<Discretization> discretization, std::shared_ptr<Partitioning> partitioning);

    //! write current velocities to file, filename is output_<count>.<rankNo>.txt
  void writeFile(double currentTime);

  std::shared_ptr<Partitioning> partitioning_;                 //< the partitioning object that knowns about the domain decomposition, only significant when executing in parallel

  //! write only current values of pressure to file, filename is pressure_<count>.<rankNo>.txt
  void writePressureFile();
};
