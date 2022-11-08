#pragma once

#include "output_writer/output_writer.h"

/** Write *.txt files that are useful for debugging.
 *  All values are written to the file as they are stored in the field variables,
 *  no interpolation takes place.
 */
class OutputWriterText : 
  public OutputWriter
{
public:
  //! use constructor of base class
  using OutputWriter::OutputWriter;

  //! write current velocities to file, filename is output_<count>.txt
  void writeFile(double currentTime);

  //! write only current values of pressure to file, filename is pressure_<count>.txt
  void writePressureFile();

};
