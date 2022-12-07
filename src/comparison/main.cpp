//#include <vtkXMLUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkImageData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <cmath>

//#include <experimental/filesystem>
//using namespace std::experimental::filesystem;

#include "boost/filesystem.hpp"
#include <iostream>
using namespace boost::filesystem;

// the data for one timestep
struct MeshData
{
  std::vector<double> u;
  std::vector<double> v;
};

// the data of one directory, i.e. for multiple timesteps
struct SimulationData
{
  std::vector<double> times;      //< time points of the parsed values
  std::vector<MeshData> values;   //< the values for a single time step, each, corresponding to times
};

void parseData(std::string filename, MeshData &meshData, double &t)
{
  // Create a reader
  vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
 
  reader->SetFileName(filename.c_str());
  reader->Update();
  

  if (reader->GetNumberOfPointArrays() == 0)
  {
    std::cout << "Error: File \"" << filename << "\" contains no data." << std::endl;
  }

  vtkSmartPointer<vtkDataSet> dataSet = reader->GetOutputAsDataSet(0);
  vtkSmartPointer<vtkFieldData> fieldData = dataSet->GetFieldData();
  
  /*
  std::cout << "file \"" << filename << "\" contains " << reader->GetNumberOfPointArrays() << " point arrays, "
    << reader->GetNumberOfCellArrays() << "," << reader->GetNumberOfColumnArrays() << std::endl;
  std::cout << " field has " << fieldData->GetNumberOfArrays() << " arrays.";
  */
  
  // get stored time
  if (fieldData->GetNumberOfArrays() > 0)
  {
    vtkSmartPointer<vtkDataArray> array = fieldData->GetArray(0);
    double *value = array->GetTuple(0);
    t = *value;
  }
  else 
  {
    std::cout << "Warning: File \"" << filename << "\" has no time." << std::endl;
  }
  
  // get stored velocity values
  vtkSmartPointer<vtkPointData> pointData = dataSet->GetPointData();
  for (int arrayNo = 0; arrayNo < pointData->GetNumberOfArrays(); arrayNo++)
  {
    vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast(pointData->GetAbstractArray(arrayNo));
    int nEntries = array->GetNumberOfTuples();
    //int nComponents = array->GetNumberOfComponents();
    //std::cout << "  " << array->GetName() << ", nEntries: " << nEntries << ", nComponents: " << nComponents << std::endl;
    
    if (std::string(array->GetName()) != "velocity")
      continue;
    
    meshData.u.resize(nEntries);
    meshData.v.resize(nEntries);
    
    // loop over values
    for (int entryNo = 0; entryNo < nEntries; entryNo++)
    {
      std::array<double,3> values;
      array->GetTuple(entryNo, values.data());
      meshData.u[entryNo] = values[0];
      meshData.v[entryNo] = values[1];
    }
  }
}

void parseDirectory(std::string directory, SimulationData &simulationData)
{
  // get files in directory
  std::vector<std::string> filenames;
  
  directory_iterator end_iter; // default construction yields past-the-end
  for (directory_iterator iter(directory); iter != end_iter; ++iter)
  {
    if (is_directory(iter->status()))
      continue;
    
    if (is_regular_file(iter->status()))
    {
      std::string filename(iter->path().string());
      
      if (filename.find(".vti") != std::string::npos)
      {
        filenames.push_back(filename);
      }
    }
  }
  
  // sort files by filename
  std::sort(filenames.begin(), filenames.end());
      
  //std::cout << "Parse data in \"" << directory << "\"." << std::endl;
  
  // loop over filenames and parse files
  for (std::vector<std::string>::iterator iter = filenames.begin(); iter != filenames.end(); iter++)
  {
    std::string filename = *iter;
    //std::cout << "  File " << filename << std::endl;
    double t;
    simulationData.values.emplace_back();
    parseData(filename, simulationData.values.back(), t);
    simulationData.times.push_back(t);
  }
}

double compareData(SimulationData &testData, SimulationData &referenceData)
{
  const bool output = false;
  
  // loop over entries in test data
  double differenceNormData = 0;
  int nFiles = testData.times.size();
  for (int i = 0; i < nFiles; i++)
  {
    double t = testData.times[i];
    MeshData &testMeshData = testData.values[i];
    
    if (output)
      std::cout << "test t = " << t << std::endl;
    
    // get corresponding data sets of reference Data
    int firstReferenceIndex = 0;
    int secondReferenceIndex = 0;
    for (int j = 0; j < (int)referenceData.times.size(); j++)
    {
      if (referenceData.times[j] > t)
      {
        secondReferenceIndex = j;
        break;
      }
      firstReferenceIndex = j;
    }
    
    double alpha = 0.0;
    
    if (secondReferenceIndex == 0)
    {
      secondReferenceIndex = firstReferenceIndex;
    }
    else
    {
      alpha = (t - referenceData.times[firstReferenceIndex]) / (referenceData.times[secondReferenceIndex] - referenceData.times[firstReferenceIndex]);
    }
    
    if (output)
    {
      std::cout << "   corresponding reference datasets: " << std::endl
        << "     t=" << referenceData.times[firstReferenceIndex] << " at index " << firstReferenceIndex << std::endl
        << "     t=" << referenceData.times[secondReferenceIndex] << " at index " << secondReferenceIndex << ", alpha: " << alpha << std::endl;
    }
    
    // loop over values
    double differenceNormDataset = 0;
    int nPoints = testMeshData.u.size();
    for (int j = 0; j < nPoints; j++)
    {
      double testU = testMeshData.u[j];
      double testV = testMeshData.v[j];
      
      double referenceU = (1.-alpha) * referenceData.values[firstReferenceIndex].u[j] + alpha * referenceData.values[secondReferenceIndex].u[j];
      double referenceV = (1.-alpha) * referenceData.values[firstReferenceIndex].v[j] + alpha * referenceData.values[secondReferenceIndex].v[j];
      
      double differenceNorm = sqrt((referenceU - testU)*(referenceU - testU) + (referenceV - testV)*(referenceV - testV));
      
      //std::cout << "      j=" << j << "/" << nPoints << ", error: " << differenceNorm << std::endl;
      differenceNormDataset += differenceNorm;
    }
    differenceNormDataset /= nPoints;
    differenceNormData += differenceNormDataset;
  }
  differenceNormData /= nFiles;
  std::cout << nFiles << " files with " << testData.values[0].u.size() << " entries each." << std::endl 
    << "average 2-norm error per velocity vector: " << differenceNormData << std::endl;
  
  if (fabs(differenceNormData) > 1e-4)
  {
    std::cout << "Error is higher than tolerance of 1e-4!" << std::endl;
  }

  return differenceNormData;
}

int main(int argc, char *argv[])
{
  // parse command line arguments
  if (argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <directory to test> <directory with reference data>" << std::endl;
    exit(-1);
  }
  
  std::string directoryTest = argv[1];
  std::string directoryReference = argv[2];

  SimulationData testData;
  SimulationData referenceData;
  
  // parse data in two directories
  parseDirectory(directoryTest, testData);
  parseDirectory(directoryReference, referenceData);
  
  if (testData.times.empty() || referenceData.times.empty())
    return -1;

  double endTimeTest = testData.times.back();
  double endTimeReference = referenceData.times.back();

  if (fabs(endTimeTest - endTimeReference) / endTimeReference > 0.1)
  {
    std::cout << "End time does not match! Test: " << endTimeTest << ", Reference: " << endTimeReference << std::endl;
  }
  
  // compare data
  compareData(testData, referenceData);
  
  return 0;
}
       
