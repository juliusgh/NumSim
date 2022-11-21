#include "output_writer/output_writer_paraview_parallel.h"

#include "storage/field_variable.h"
#include "discretization/1_discretization.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <mpi.h>

OutputWriterParaviewParallel::OutputWriterParaviewParallel(std::shared_ptr<Discretization> discretization, const Partitioning &partitioning) :
   OutputWriter(discretization, partitioning),

  nCellsGlobal_(partitioning_.nCellsGlobal()),
  nPointsGlobal_ {nCellsGlobal_[0]+1, nCellsGlobal_[1]+1},    // we have one point more than cells in every coordinate direction
  
  // create field variables for resulting values, only for local data as send buffer
  u_(nPointsGlobal_, std::array<double,2>{0.,0.}, discretization_->meshWidth()),
  v_(nPointsGlobal_, std::array<double,2>{0.,0.}, discretization_->meshWidth()),
  p_(nPointsGlobal_, std::array<double,2>{0.,0.}, discretization_->meshWidth()),
  
  // create field variables for resulting values, after MPI communication
  uGlobal_(nPointsGlobal_, std::array<double,2>{0.,0.}, discretization_->meshWidth()),
  vGlobal_(nPointsGlobal_, std::array<double,2>{0.,0.}, discretization_->meshWidth()),
  pGlobal_(nPointsGlobal_, std::array<double,2>{0.,0.}, discretization_->meshWidth())
{
  // Create a vtkWriter_
  vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}

void OutputWriterParaviewParallel::gatherData()
{
 // std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth

  int nPointsGlobalTotal = nPointsGlobal_[0] * nPointsGlobal_[1];
  const double dx = discretization_->meshWidth()[0];
  const double dy = discretization_->meshWidth()[1];

  // set values in own subdomain, other values are left at zero

  // determine data range {0,…,iEnd-1} x {0,…,jEnd-1}
  std::array<int,2> nCells = discretization_->nCells();
  int jEnd = nCells[1];
  int iEnd = nCells[0];

  // add right-most points at ranks with right boundary
  if (partitioning_.ownPartitionContainsRightBoundary())
    iEnd += 1;

  // add right-most points at ranks with top boundary
  if (partitioning_.ownPartitionContainsTopBoundary())
    jEnd += 1;

  std::array<int,2> nodeOffset = partitioning_.nodeOffset();

  u_.setToZero();
  v_.setToZero();
  p_.setToZero();

  for (int j = 0; j < jEnd; j++)
  {
    for (int i = 0; i < iEnd; i++)
    {
      const double x = i*dx;
      const double y = j*dy;

      // get global indices
      int iGlobal = nodeOffset[0] + i;
      int jGlobal = nodeOffset[1] + j;

      u_(iGlobal,jGlobal) = discretization_->u().interpolateAt(x,y);
      v_(iGlobal,jGlobal) = discretization_->v().interpolateAt(x,y);
      p_(iGlobal,jGlobal) = discretization_->p().interpolateAt(x,y);
    }
  }

  // sum up values from all ranks, not set values are zero
  MPI_Reduce(u_.data(), uGlobal_.data(), nPointsGlobalTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(v_.data(), vGlobal_.data(), nPointsGlobalTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(p_.data(), pGlobal_.data(), nPointsGlobalTotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

}

void OutputWriterParaviewParallel::writeFile(double currentTime)
{
  // communicate all data to rank 0
  gatherData();

  // only continue to write the file on rank 0
  if (partitioning_.ownRankNo() != 0)
  {
    return;
  }

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNo_ << "." << vtkWriter_->GetDefaultFileExtension();
  
  // increment file no.
  fileNo_++;

  // assign the new file name to the output vtkWriter_
  vtkWriter_->SetFileName(fileName.str().c_str());
  
  // initialize data set that will be output to the file
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  const double dx = discretization_->meshWidth()[0];
  const double dy = discretization_->meshWidth()[1];
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(nCellsGlobal_[0]+1, nCellsGlobal_[1]+1, 1);  // we want to have points at each corner of each cell
  

  // add pressure field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  arrayPressure->SetName("pressure");

  // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  int index = 0;   // index for the vtk data structure, will be incremented in the inner loop
  for (int j = 0; j < nCellsGlobal_[1]+1; j++)
  {
    for (int i = 0; i < nCellsGlobal_[0]+1; i++, index++)
    {
      arrayPressure->SetValue(index, pGlobal_(i,j));
    }
  }

  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayPressure);
  
  // add velocity field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayVelocity = vtkDoubleArray::New();

  // here we have two components (u,v), but ParaView will only allow vector glyphs if we have an ℝ^3 vector, 
  // therefore we use a 3-dimensional vector and set the 3rd component to zero
  arrayVelocity->SetNumberOfComponents(3);

  // set the number of values
  arrayVelocity->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  arrayVelocity->SetName("velocity");

  // loop over the mesh where p is defined and assign the values in the vtk data structure
  index = 0;   // index for the vtk data structure
  for (int j = 0; j < nCellsGlobal_[1]+1; j++)
  {
    const double y = j*dy;

    for (int i = 0; i < nCellsGlobal_[0]+1; i++, index++)
    {
      const double x = i*dx;

      std::array<double,3> velocityVector;
      velocityVector[0] = uGlobal_(i,j);
      velocityVector[1] = vGlobal_(i,j);
      velocityVector[2] = 0.0;    // z-direction is 0

      arrayVelocity->SetTuple(index, velocityVector.data());
    }
  }
  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayVelocity);
  
  // add current time 
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, currentTime);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();
  
  // Write the data
  vtkWriter_->SetInputData(dataSet);
  
  //vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii text files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary();      // set file mode to binary files: smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}