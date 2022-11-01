#include "output_writer/output_writer_paraview.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <memory>

OutputWriterParaview::OutputWriterParaview(std::shared_ptr<Discretization> discretization) :
        OutputWriter(discretization)
{
    vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}

void OutputWriterParaview::writeFile(double currentTime)
{

}