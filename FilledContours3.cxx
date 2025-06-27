#include <vtkActor.h>
#include <vtkAlgorithm.h>
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyData.h>

#include <json/json.h> // For JSON parsing

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

// 读取json成功

// Function to read JSON data
vtkSmartPointer<vtkPolyData> ReadJSONData(const std::string& filename, double& zMin, double& zMax)
{
  // Read JSON file
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return nullptr;
  }

  Json::Value root;
  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  bool ok = Json::parseFromStream(builder, file, &root, &errs);
  if (!ok) {
    std::cerr << "Error parsing JSON: " << errs << std::endl;
    return nullptr;
  }

  // Extract data from JSON
  const Json::Value& data = root["data"];
  const Json::Value& xArray = data["x"];
  const Json::Value& yArray = data["y"];
  const Json::Value& vArray = data["v"];
  
  if (data.isMember("zmin") && data.isMember("zmax")) {
    zMin = data["zmin"].asDouble();
    zMax = data["zmax"].asDouble();
  } else {
    zMin = std::numeric_limits<double>::max();
    zMax = std::numeric_limits<double>::lowest();
  }

  int numX = xArray.size();
  int numY = yArray.size();
  
  // Create points for structured grid
  vtkNew<vtkPoints> points;
  vtkNew<vtkDoubleArray> scalarArray;
  scalarArray->SetName("Elevation");
  
  // Initialize with appropriate size
  points->SetNumberOfPoints(numX * numY);
  scalarArray->SetNumberOfValues(numX * numY);
  
  // Fill the points and scalar values
  int pointId = 0;
  for (int j = 0; j < numY; j++) {
    double y = yArray[j].asDouble();
    for (int i = 0; i < numX; i++) {
      double x = xArray[i].asDouble();
      double z = vArray[j][i].asDouble();
      
      // Update min/max if needed
      if (z < zMin) zMin = z;
      if (z > zMax) zMax = z;
      
      points->SetPoint(pointId, x, y, 0.0); // Create 2D points (z=0)
      scalarArray->SetValue(pointId, z);
      pointId++;
    }
  }
  
  // Create a Delaunay triangulation
  vtkNew<vtkDelaunay2D> delaunay;
  
  // Create polydata to hold the points
  vtkNew<vtkPolyData> polydata;
  polydata->SetPoints(points);
  
  // Perform Delaunay triangulation
  delaunay->SetInputData(polydata);
  delaunay->Update();
  
  // Get the output and add scalar values
  vtkSmartPointer<vtkPolyData> output = delaunay->GetOutput();
  output->GetPointData()->SetScalars(scalarArray);
  
  // Clean the output to ensure proper topology
  vtkNew<vtkCleanPolyData> cleaner;
  cleaner->SetInputData(output);
  cleaner->Update();
  
  return cleaner->GetOutput();
}

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr
        << "Usage: " << argv[0]
        << " InputJSONFile(.json) NumberOfContours [SmoothingLevel=0]" << std::endl
        << "Example: FilledContours.exe testData.json 10 3" << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkNamedColors> colors;

  // Read the JSON file
  double scalarRange[2];
  vtkSmartPointer<vtkPolyData> inputPolyData = ReadJSONData(argv[1], scalarRange[0], scalarRange[1]);
  
  if (!inputPolyData) {
    std::cerr << "Failed to read JSON data" << std::endl;
    return EXIT_FAILURE;
  }

  int numberOfContours = atoi(argv[2]);
  
  // Smoothing level, default is 0 (no smoothing)
  int smoothingLevel = 0;
  if (argc > 3)
  {
    smoothingLevel = atoi(argv[3]);
  }

  double delta = (scalarRange[1] - scalarRange[0]) /
      static_cast<double>(numberOfContours - 1);

  vtkNew<vtkAppendPolyData> appendFilledContours;

  // Keep the clippers alive.
  std::vector<vtkSmartPointer<vtkClipPolyData>> clippersLo;
  std::vector<vtkSmartPointer<vtkClipPolyData>> clippersHi;

  for (int i = 0; i < numberOfContours; i++)
  {
    double valueLo = scalarRange[0] + static_cast<double>(i) * delta;
    double valueHi = scalarRange[0] + static_cast<double>(i + 1) * delta;
    
    clippersLo.push_back(vtkSmartPointer<vtkClipPolyData>::New());
    clippersLo[i]->SetValue(valueLo);
    
    if (i == 0)
    {
      clippersLo[i]->SetInputData(inputPolyData);
    }
    else
    {
      clippersLo[i]->SetInputConnection(clippersHi[i - 1]->GetOutputPort(1));
    }
    
    clippersLo[i]->InsideOutOff();
    clippersLo[i]->Update();

    clippersHi.push_back(vtkSmartPointer<vtkClipPolyData>::New());
    clippersHi[i]->SetValue(valueHi);
    clippersHi[i]->SetInputConnection(clippersLo[i]->GetOutputPort());
    clippersHi[i]->GenerateClippedOutputOn();
    clippersHi[i]->InsideOutOn();
    clippersHi[i]->Update();
    
    if (clippersHi[i]->GetOutput()->GetNumberOfCells() == 0)
    {
      continue;
    }

    vtkNew<vtkFloatArray> cd;
    cd->SetNumberOfComponents(1);
    cd->SetNumberOfTuples(clippersHi[i]->GetOutput()->GetNumberOfCells());
    cd->FillComponent(0, valueLo);

    clippersHi[i]->GetOutput()->GetCellData()->SetScalars(cd);
    appendFilledContours->AddInputConnection(clippersHi[i]->GetOutputPort());
  }

  vtkNew<vtkCleanPolyData> filledContours;
  filledContours->SetInputConnection(appendFilledContours->GetOutputPort());
  
  // Store the final data source, which may be original or smoothed
  vtkSmartPointer<vtkAlgorithm> finalDataSource;
  
  // First convert quads to triangles
  vtkNew<vtkTriangleFilter> triangleFilter;
  triangleFilter->SetInputConnection(filledContours->GetOutputPort());
  triangleFilter->Update();
  
  // Apply smoothing to the contour surface
  if (smoothingLevel > 0)
  {
    // Use simple vtkSmoothPolyDataFilter which is more robust with non-manifold data
    vtkNew<vtkSmoothPolyDataFilter> smoothFilter;
    smoothFilter->SetInputConnection(triangleFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(smoothingLevel * 10);
    smoothFilter->SetRelaxationFactor(0.2);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    
    finalDataSource = smoothFilter;
  }
  else
  {
    finalDataSource = triangleFilter;
  }

  // Set up the contour surface mapper
  vtkNew<vtkPolyDataMapper> contourMapper;
  contourMapper->SetInputConnection(finalDataSource->GetOutputPort());
  
  vtkNew<vtkLookupTable> lut;
  lut->SetNumberOfTableValues(numberOfContours + 1);
  lut->Build();
  
  contourMapper->SetScalarRange(scalarRange[0], scalarRange[1]);
  contourMapper->SetScalarModeToUseCellData();
  contourMapper->SetLookupTable(lut);

  vtkNew<vtkActor> contourActor;
  contourActor->SetMapper(contourMapper);
  contourActor->GetProperty()->SetInterpolationToFlat();

  // Generate contour lines from the same data source
  vtkNew<vtkContourFilter> contours;
  contours->SetInputConnection(finalDataSource->GetOutputPort());
  contours->GenerateValues(numberOfContours, scalarRange[0], scalarRange[1]);

  vtkNew<vtkPolyDataMapper> contourLineMapperer;
  contourLineMapperer->SetInputConnection(contours->GetOutputPort());
  contourLineMapperer->SetScalarRange(scalarRange[0], scalarRange[1]);
  contourLineMapperer->ScalarVisibilityOff();

  vtkNew<vtkActor> contourLineActor;
  contourLineActor->SetMapper(contourLineMapperer);
  contourLineActor->GetProperty()->SetLineWidth(2);

  // The usual renderer, render window and interactor.
  vtkNew<vtkRenderer> ren1;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  renWin->AddRenderer(ren1);
  renWin->SetWindowName("FilledContours");

  iren->SetRenderWindow(renWin);

  // Add the actors
  ren1->AddActor(contourActor);
  ren1->AddActor(contourLineActor);
  ren1->SetBackground(colors->GetColor3d("MidnightBlue").GetData());

  // Begin interaction.
  renWin->Render();
  iren->Start();

  return EXIT_SUCCESS;
}