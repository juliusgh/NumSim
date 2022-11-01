#pragma once

#include "output_writer/output_writer.h"

class OutputWriterText : OutputWriter {
public:
    //! constructor
    using OutputWriter::OutputWriter;

    //! write current velocities to file, filename is output_<count>.vti
    void writeFile(double currentTime);

    //! write only current values of pressure to file, filename is pressure_<count>.txt
    void writePressureFile();
};
