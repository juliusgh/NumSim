#pragma once

#include "output_writer/output_writer.h"
#include "discretization/1_discretization.h"

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

//! Write a file out/output_<fileNo>.vti to be visualized in ParaView.
//! It contains 10x10 nodes with an artifical pressure field.
//! This method is only for demonstration purpose and does nothing useful.
//! However, we will provide similar files, e.g. "output_writer_paraview.h", to be used in the submission code.


class OutputWriterParaview : OutputWriter {
public:
    //! constructor
    OutputWriterParaview(std::shared_ptr <Discretization> discretization);

    //! write current velocities to file, filename is output_<count>.vti
    void writeFile(double currentTime);

private:
    vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;   //< vtk writer to write ImageData
};
