#pragma once

#include "output_writer/output_writer.h"
#include "discretization/1_discretization.h"

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <memory>

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the computational domain.
 *  All values are given for the nodes of the mesh, i.e., the corners of each cell.
 *  This means, values will be interpolated because the values are stored at positions given by the staggered grid.
 */
class OutputWriterParaviewParallel : 
  public OutputWriter
{
public:
  //! constructor
  OutputWriterParaviewParallel(std::shared_ptr<Discretization> discretization, const Partitioning &partitioning);

  //! write current velocities to file, filename is output_<count>.vti
  void writeFile(double currentTime);

private:

  //! gather u,v and p values from all ranks to rank 0 and store them in the global field variables
  void gatherData();

  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;   //< vtk writer to write ImageData

  std::array<int,2> nCellsGlobal_;   //< global number of cells
  std::array<int,2> nPointsGlobal_;  //< global number of points

  FieldVariable u_;    // field variable for u with global size, contains only the local values, other entries are 0
  FieldVariable v_;    // field variable for v with global size, contains only the local values, other entries are 0
  FieldVariable p_;    // field variable for p with global size, contains only the local values, other entries are 0

  FieldVariable uGlobal_;    // on rank 0: field variable for u that gathers values from all ranks, on other ranks: nullptr
  FieldVariable vGlobal_;    // on rank 0: field variable for v that gathers values from all ranks, on other ranks: nullptr
  FieldVariable pGlobal_;    // on rank 0: field variable for p that gathers values from all ranks, on other ranks: nullptr
};
