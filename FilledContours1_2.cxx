#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkContourFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkFeatureEdges.h>
#include <vtkSmoothPolyDataFilter.h>

#include <json/json.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// This program loads JSON data with null values and creates filled contours visualization
// The null value boundaries are used as masks for contour lines and filled contours
// 原始的demo 加载json文件 json中可以包含空值

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr
        << "Usage: " << argv[0]
        << " InputJSONFile(.json) NumberOfContours e.g testData1.json 10"
        << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkNamedColors> colors;

  // Read the JSON file
  std::ifstream jsonFile(argv[1]);
  if (!jsonFile.is_open()) {
    std::cerr << "Failed to open JSON file: " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  Json::Value root;
  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  if (!Json::parseFromStream(builder, jsonFile, &root, &errs)) {
    std::cerr << "Error parsing JSON: " << errs << std::endl;
    return EXIT_FAILURE;
  }

  // Extract data from JSON
  Json::Value& data = root["data"];
  Json::Value& x = data["x"];
  Json::Value& y = data["y"];
  Json::Value& v = data["v"];
  double zmin = data.get("zmin", 0.0).asDouble();
  double zmax = data.get("zmax", 1.0).asDouble();

  // Create points and values
  vtkNew<vtkPoints> points;
  vtkNew<vtkFloatArray> scalars;
  scalars->SetName("Elevation");

  // Create a mask for valid points
  std::vector<bool> validPoint;

  // Process the data
  for (int j = 0; j < y.size(); ++j) {
    double yCoord = y[j].asDouble();
    for (int i = 0; i < x.size(); ++i) {
      double xCoord = x[i].asDouble();
      double value = 0.0;
      
      bool isValid = !v[j][i].isNull();
      validPoint.push_back(isValid);
      
      if (isValid) {
        value = v[j][i].asDouble();
      } else {
        // Skip null values for now
        continue;
      }
      
      // Add point and scalar value
      points->InsertNextPoint(xCoord, yCoord, 0.0);
      scalars->InsertNextValue(value);
    }
  }

  // Create polydata
  vtkNew<vtkPolyData> polydata;
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(scalars);

  // Use Delaunay triangulation to create a surface from points
  vtkNew<vtkDelaunay2D> delaunay;
  delaunay->SetInputData(polydata);
  delaunay->SetTolerance(0.001);
  delaunay->Update();

  // Clip the triangulation to remove triangles in null regions
  vtkNew<vtkAppendPolyData> appendFilledContours;

  int numberOfContours = atoi(argv[2]);
  
  double scalarRange[2];
  scalarRange[0] = zmin;
  scalarRange[1] = zmax;

  double delta = (scalarRange[1] - scalarRange[0]) /
      static_cast<double>(numberOfContours - 1);

  // Keep the clippers alive
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
      clippersLo[i]->SetInputConnection(delaunay->GetOutputPort());
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

  // Clean the data
  vtkNew<vtkCleanPolyData> filledContours;
  filledContours->SetInputConnection(appendFilledContours->GetOutputPort());

  // Create a lookup table for coloring
  vtkNew<vtkLookupTable> lut;
  lut->SetNumberOfTableValues(numberOfContours + 1);
  lut->Build();

  // Set custom color for each contour band
  for (int i = 0; i < numberOfContours; i++) {
    double t = static_cast<double>(i) / numberOfContours;
    double r = 0.0, g = 0.0, b = 1.0;
    
    // Simple blue to red gradient
    r = t;
    b = 1.0 - t;
    
    lut->SetTableValue(i, r, g, b);
  }

  // Create the mapper for filled contours
  vtkNew<vtkPolyDataMapper> contourMapper;
  contourMapper->SetInputConnection(filledContours->GetOutputPort());
  contourMapper->SetScalarRange(scalarRange[0], scalarRange[1]);
  contourMapper->SetScalarModeToUseCellData();
  contourMapper->SetLookupTable(lut);

  vtkNew<vtkActor> contourActor;
  contourActor->SetMapper(contourMapper);
  contourActor->GetProperty()->SetInterpolationToFlat();

  // Create contour lines
  vtkNew<vtkContourFilter> contours;
  contours->SetInputConnection(filledContours->GetOutputPort());
  contours->GenerateValues(numberOfContours, scalarRange[0], scalarRange[1]);

  // Smooth the contours
  vtkNew<vtkSmoothPolyDataFilter> smoothContours;
  smoothContours->SetInputConnection(contours->GetOutputPort());
  smoothContours->SetNumberOfIterations(2);
  smoothContours->SetRelaxationFactor(0.2);
  smoothContours->Update();

  vtkNew<vtkPolyDataMapper> contourLineMapperer;
  contourLineMapperer->SetInputConnection(smoothContours->GetOutputPort());
  contourLineMapperer->SetScalarRange(scalarRange[0], scalarRange[1]);
  contourLineMapperer->ScalarVisibilityOff();

  vtkNew<vtkActor> contourLineActor;
  contourLineActor->SetMapper(contourLineMapperer);
  contourLineActor->GetProperty()->SetLineWidth(2);
  contourLineActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

  // Create boundary edges around valid data regions
  vtkNew<vtkFeatureEdges> boundaryEdges;
  boundaryEdges->SetInputConnection(delaunay->GetOutputPort());
  boundaryEdges->BoundaryEdgesOn();
  boundaryEdges->FeatureEdgesOff();
  boundaryEdges->NonManifoldEdgesOff();
  boundaryEdges->ManifoldEdgesOff();
  boundaryEdges->Update();

  vtkNew<vtkPolyDataMapper> boundaryMapper;
  boundaryMapper->SetInputConnection(boundaryEdges->GetOutputPort());
  
  vtkNew<vtkActor> boundaryActor;
  boundaryActor->SetMapper(boundaryMapper);
  boundaryActor->GetProperty()->SetLineWidth(3);
  boundaryActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

  // The usual renderer, render window and interactor
  vtkNew<vtkRenderer> ren1;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  renWin->AddRenderer(ren1);
  renWin->SetWindowName("FilledContours");

  iren->SetRenderWindow(renWin);

  // Add the actors
  ren1->AddActor(contourActor);
  ren1->AddActor(contourLineActor);
  ren1->AddActor(boundaryActor);
  ren1->SetBackground(colors->GetColor3d("MidnightBlue").GetData());

  // Begin interaction
  renWin->Render();
  iren->Start();

  return EXIT_SUCCESS;
}
