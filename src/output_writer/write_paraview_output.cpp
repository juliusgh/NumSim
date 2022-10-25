#include "output_writer/write_paraview_output.h"

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <cstdlib>
#include <iostream>

void writeParaviewOutput(int fileNo)
{
  // create "out" subdirectory if it does not yet exist
  int returnValue = system("mkdir -p out");
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl;

  // Create a vtkWriter
  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNo << "." << vtkWriter->GetDefaultFileExtension();

  std::cout << "Write file \"" << fileName.str() << "\"." << std::endl;

  // assign the new file name to the output vtkWriter
  vtkWriter->SetFileName(fileName.str().c_str());

  // initialize data set that will be output to the file
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  const double dx = 1;
  const double dy = 1;
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(10, 10, 1);

  // add pressure field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayPressure->SetName("pressure");

  // loop over the nodes of the mesh and assign an artifical value that changes with fileNo
  for (int j = 0; j < 10; j++)
  {
    for (int i = 0; i < 10; i++)
    {
      int index = j*10 + i;
      arrayPressure->SetValue(index, i+j-fileNo*(i-j));
    }
  }

  // Add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayPressure);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  vtkWriter->SetInputData(dataSet);

  //vtkWriter->SetDataModeToAscii();     // comment this in to get ascii text files: those can be checked in an editor
  vtkWriter->SetDataModeToBinary();      // set file mode to binary files: smaller file sizes

  // finally write out the data
  vtkWriter->Write();
}

