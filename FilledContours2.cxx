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
#include <vtkXMLPolyDataReader.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkButterflySubdivisionFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkLinearSubdivisionFilter.h>

#include <iostream>
#include <string>
#include <vector>



 // 使用了平滑参数  但是感觉作用不大
int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr
        << "Usage: " << argv[0]
        << " InputPolyDataFile(.vtp) NumberOfContours [SmoothingLevel=0]" << std::endl
        << "Example: FilledContours.vtp 10 3" << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkNamedColors> colors;

  // Read the file
  vtkNew<vtkXMLPolyDataReader> reader;

  reader->SetFileName(argv[1]);
  reader->Update(); // Update so that we can get the scalar range.

  double scalarRange[2];
  reader->GetOutput()->GetPointData()->GetScalars()->GetRange(scalarRange);

  vtkNew<vtkAppendPolyData> appendFilledContours;

  int numberOfContours = atoi(argv[2]);
  
  // Smoothing level, default is 0 (no smoothing)
  int smoothingLevel = 0;
  if (argc > 3)
  {
    smoothingLevel = atoi(argv[3]);
  }

  double delta = (scalarRange[1] - scalarRange[0]) /
      static_cast<double>(numberOfContours - 1);

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
      clippersLo[i]->SetInputConnection(reader->GetOutputPort());
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
    // Use linear subdivision for smoother results
    vtkNew<vtkLinearSubdivisionFilter> subdivisionFilter;
    subdivisionFilter->SetInputConnection(triangleFilter->GetOutputPort());
    subdivisionFilter->SetNumberOfSubdivisions(1);
    subdivisionFilter->Update();
    
    // Use windowed sinc filter for final smoothing with preserved features
    vtkNew<vtkWindowedSincPolyDataFilter> smoothFilter;
    smoothFilter->SetInputConnection(subdivisionFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(smoothingLevel * 5); // Adjust iterations based on smoothing level
    smoothFilter->SetPassBand(0.1); // Lower values produce smoother results
    smoothFilter->SetFeatureEdgeSmoothing(false);
    smoothFilter->SetBoundarySmoothing(true);
    smoothFilter->SetNormalizeCoordinates(true);
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