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
#include <vtkPolygon.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>

#include <json/json.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

// This program loads JSON data with null values and creates filled contours visualization
// The null value boundaries are used as masks for contour lines and filled contours
// 原始的demo 加载json文件 然后导出成 geojson 格式的数据 

// Function to export contour lines as GeoJSON
void ExportContoursToGeoJSON(vtkPolyData* contours, const std::string& filename, double zMin, double zMax) {
    // Create GeoJSON root object
    Json::Value geoJSON;
    geoJSON["type"] = "FeatureCollection";
    geoJSON["features"] = Json::Value(Json::arrayValue);

    // Get contour values from point data
    vtkDataArray* scalarData = contours->GetPointData()->GetScalars();
    
    // Get cell data (for line segments)
    vtkCellArray* lines = contours->GetLines();
    
    if (!lines || lines->GetNumberOfCells() == 0) {
        std::cout << "No contour lines found to export." << std::endl;
        return;
    }
    
    // Maps to group lines by contour value
    std::map<double, std::vector<std::vector<double*>>> contourMap;
    
    // For each cell (line)
    vtkIdType numCells = contours->GetNumberOfCells();
    for (vtkIdType cellId = 0; cellId < numCells; cellId++) {
        vtkCell* cell = contours->GetCell(cellId);
        
        if (cell->GetCellType() != VTK_LINE) {
            continue; // Skip non-line cells
        }
        
        // Get the first point of the line to determine contour value
        vtkIdType pointId = cell->GetPointId(0);
        double contourValue = scalarData->GetTuple1(pointId);
        
        // Create a line segment
        std::vector<double*> lineCoords;
        for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); i++) {
            vtkIdType pid = cell->GetPointId(i);
            double* point = new double[3];
            contours->GetPoint(pid, point);
            lineCoords.push_back(point);
        }
        
        // Add to the contour map
        contourMap[contourValue].push_back(lineCoords);
    }
    
    // For each contour value, create a GeoJSON feature
    for (auto& contourPair : contourMap) {
        double contourValue = contourPair.first;
        auto& contourLines = contourPair.second;
        
        // Normalize value for color computation
        double normalizedValue = (contourValue - zMin) / (zMax - zMin);
        
        // Create feature for this contour value
        Json::Value feature;
        feature["type"] = "Feature";
        
        // Add properties
        feature["properties"] = Json::Value(Json::objectValue);
        feature["properties"]["contourValue"] = contourValue;
        feature["properties"]["normalizedValue"] = normalizedValue;
        
        // Create MultiLineString geometry
        feature["geometry"] = Json::Value(Json::objectValue);
        feature["geometry"]["type"] = "MultiLineString";
        feature["geometry"]["coordinates"] = Json::Value(Json::arrayValue);
        
        // Add all line segments for this contour value
        for (auto& line : contourLines) {
            Json::Value lineCoords(Json::arrayValue);
            for (double* point : line) {
                Json::Value coord(Json::arrayValue);
                coord.append(point[0]); // x
                coord.append(point[1]); // y
                lineCoords.append(coord);
                delete[] point; // Clean up memory
            }
            feature["geometry"]["coordinates"].append(lineCoords);
        }
        
        // Add feature to collection
        geoJSON["features"].append(feature);
    }
    
    // Write to file
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        Json::StyledWriter writer;
        outFile << writer.write(geoJSON);
        outFile.close();
        std::cout << "Contour lines exported to " << filename << std::endl;
    } else {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
    }
}

// Function to export filled contours as GeoJSON
void ExportFilledContoursToGeoJSON(vtkPolyData* filledContours, const std::string& filename, double zMin, double zMax) {
    // First, ensure we have polygons by triangulating if needed
    vtkNew<vtkTriangleFilter> triangleFilter;
    triangleFilter->SetInputData(filledContours);
    triangleFilter->Update();
    
    // Create GeoJSON root object
    Json::Value geoJSON;
    geoJSON["type"] = "FeatureCollection";
    geoJSON["features"] = Json::Value(Json::arrayValue);
    
    // Get scalar data (contour values)
    vtkDataArray* scalarData = triangleFilter->GetOutput()->GetCellData()->GetScalars();
    
    if (!scalarData) {
        std::cout << "No scalar data found for filled contours." << std::endl;
        return;
    }
    
    // Group cells by scalar value (contour band)
    std::map<double, vtkSmartPointer<vtkPolyData>> contourBands;
    
    vtkIdType numCells = triangleFilter->GetOutput()->GetNumberOfCells();
    for (vtkIdType i = 0; i < numCells; i++) {
        double contourValue = scalarData->GetTuple1(i);
        
        // Create a new polydata for this contour value if it doesn't exist
        if (contourBands.find(contourValue) == contourBands.end()) {
            contourBands[contourValue] = vtkSmartPointer<vtkPolyData>::New();
            contourBands[contourValue]->SetPoints(triangleFilter->GetOutput()->GetPoints());
            contourBands[contourValue]->Allocate(numCells); // Preallocate
        }
        
        // Add this cell to the polydata for this contour value
        vtkCell* cell = triangleFilter->GetOutput()->GetCell(i);
        vtkIdList* idList = cell->GetPointIds();
        contourBands[contourValue]->InsertNextCell(cell->GetCellType(), idList);
    }
    
    // Process each contour band
    for (auto& bandPair : contourBands) {
        double contourValue = bandPair.first;
        vtkPolyData* band = bandPair.second;
        
        // Use connectivity filter to separate disconnected regions
        vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
        connectivityFilter->SetInputData(band);
        connectivityFilter->SetExtractionModeToAllRegions();
        connectivityFilter->Update();
        
        // Normalize value for properties
        double normalizedValue = (contourValue - zMin) / (zMax - zMin);
        
        // Create feature for this contour band
        Json::Value feature;
        feature["type"] = "Feature";
        
        // Add properties
        feature["properties"] = Json::Value(Json::objectValue);
        feature["properties"]["contourValue"] = contourValue;
        feature["properties"]["normalizedValue"] = normalizedValue;
        
        // Create MultiPolygon geometry
        feature["geometry"] = Json::Value(Json::objectValue);
        feature["geometry"]["type"] = "MultiPolygon";
        feature["geometry"]["coordinates"] = Json::Value(Json::arrayValue);
        
        // Extract boundaries of each region
        vtkNew<vtkFeatureEdges> boundaryEdges;
        boundaryEdges->SetInputData(connectivityFilter->GetOutput());
        boundaryEdges->BoundaryEdgesOn();
        boundaryEdges->FeatureEdgesOff();
        boundaryEdges->NonManifoldEdgesOff();
        boundaryEdges->ManifoldEdgesOff();
        boundaryEdges->Update();
        
        // Group boundary edges into continuous loops
        vtkPolyData* boundaries = boundaryEdges->GetOutput();
        vtkCellArray* lines = boundaries->GetLines();
        
        if (!lines || lines->GetNumberOfCells() == 0) {
            continue; // Skip if no boundaries
        }
        
        // Extract and organize boundaries into polygon rings
        std::vector<std::vector<std::vector<double>>> polygons;
        std::vector<std::vector<double>> currentRing;
        
        vtkNew<vtkIdList> pointIds;
        lines->InitTraversal();
        
        // Set to track processed cells
        std::set<vtkIdType> processedCells;
        
        // Process each cell
        for (vtkIdType cellId = 0; cellId < boundaries->GetNumberOfCells(); cellId++) {
            if (processedCells.find(cellId) != processedCells.end()) {
                continue; // Skip processed cells
            }
            
            // Start a new ring
            currentRing.clear();
            
            // Get the current cell's points
            boundaries->GetCellPoints(cellId, pointIds);
            
            // Start point
            double point[3];
            boundaries->GetPoint(pointIds->GetId(0), point);
            currentRing.push_back({point[0], point[1]});
            
            // End point of first segment
            boundaries->GetPoint(pointIds->GetId(1), point);
            currentRing.push_back({point[0], point[1]});
            
            processedCells.insert(cellId);
            
            // Keep connecting segments until we form a closed loop
            bool closedLoop = false;
            vtkIdType startPointId = pointIds->GetId(0);
            vtkIdType currentEndPointId = pointIds->GetId(1);
            
            while (!closedLoop) {
                bool foundNextSegment = false;
                
                // Look for a segment that connects to the current end point
                for (vtkIdType nextCellId = 0; nextCellId < boundaries->GetNumberOfCells(); nextCellId++) {
                    if (processedCells.find(nextCellId) != processedCells.end()) {
                        continue; // Skip processed cells
                    }
                    
                    vtkNew<vtkIdList> nextPointIds;
                    boundaries->GetCellPoints(nextCellId, nextPointIds);
                    
                    if (nextPointIds->GetId(0) == currentEndPointId) {
                        // This segment connects to our current end
                        boundaries->GetPoint(nextPointIds->GetId(1), point);
                        currentRing.push_back({point[0], point[1]});
                        
                        currentEndPointId = nextPointIds->GetId(1);
                        processedCells.insert(nextCellId);
                        foundNextSegment = true;
                        
                        // Check if we've closed the loop
                        if (currentEndPointId == startPointId) {
                            closedLoop = true;
                        }
                        
                        break;
                    }
                    else if (nextPointIds->GetId(1) == currentEndPointId) {
                        // This segment connects to our current end (in reverse)
                        boundaries->GetPoint(nextPointIds->GetId(0), point);
                        currentRing.push_back({point[0], point[1]});
                        
                        currentEndPointId = nextPointIds->GetId(0);
                        processedCells.insert(nextCellId);
                        foundNextSegment = true;
                        
                        // Check if we've closed the loop
                        if (currentEndPointId == startPointId) {
                            closedLoop = true;
                        }
                        
                        break;
                    }
                }
                
                if (!foundNextSegment) {
                    // If we can't find a connecting segment, we have an open loop - handle as needed
                    break;
                }
            }
            
            // Close the ring if needed by duplicating the first point
            if (closedLoop && (currentRing.front()[0] != currentRing.back()[0] || 
                              currentRing.front()[1] != currentRing.back()[1])) {
                currentRing.push_back(currentRing.front());
            }
            
            // If we have a valid ring, add it to polygons
            if (currentRing.size() >= 4) { // Minimum 3 points plus closing point
                std::vector<std::vector<double>> ringCopy = currentRing;
                polygons.push_back(ringCopy);
            }
        }
        
        // Convert polygon rings to GeoJSON format
        // First polygon is the outer boundary, any others are holes
        if (!polygons.empty()) {
            // Calculate polygon areas to determine outer vs. inner rings
            std::vector<double> areas(polygons.size());
            for (size_t i = 0; i < polygons.size(); i++) {
                double area = 0.0;
                for (size_t j = 0; j < polygons[i].size() - 1; j++) {
                    area += (polygons[i][j][0] * polygons[i][j+1][1] - 
                             polygons[i][j+1][0] * polygons[i][j][1]);
                }
                areas[i] = area / 2.0;
            }
            
            // Group rings into polygons (outer ring + holes)
            std::vector<size_t> outerRings;
            for (size_t i = 0; i < areas.size(); i++) {
                if (areas[i] > 0) {
                    outerRings.push_back(i);
                }
            }
            
            // For each outer ring, find its holes
            for (size_t outerIdx : outerRings) {
                Json::Value polygonCoords(Json::arrayValue);
                
                // Add the outer ring
                Json::Value outerRingCoords(Json::arrayValue);
                for (const auto& point : polygons[outerIdx]) {
                    Json::Value coord(Json::arrayValue);
                    coord.append(point[0]);
                    coord.append(point[1]);
                    outerRingCoords.append(coord);
                }
                polygonCoords.append(outerRingCoords);
                
                // Add holes (inner rings)
                for (size_t i = 0; i < polygons.size(); i++) {
                    if (i != outerIdx && areas[i] < 0) {
                        // This is a hole - add it
                        Json::Value innerRingCoords(Json::arrayValue);
                        for (const auto& point : polygons[i]) {
                            Json::Value coord(Json::arrayValue);
                            coord.append(point[0]);
                            coord.append(point[1]);
                            innerRingCoords.append(coord);
                        }
                        polygonCoords.append(innerRingCoords);
                    }
                }
                
                // Add this polygon to the MultiPolygon
                feature["geometry"]["coordinates"].append(polygonCoords);
            }
        }
        
        // Add feature to the collection
        geoJSON["features"].append(feature);
    }
    
    // Write to file
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        Json::StyledWriter writer;
        outFile << writer.write(geoJSON);
        outFile.close();
        std::cout << "Filled contours exported to " << filename << std::endl;
    } else {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
    }
}

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
  filledContours->Update();

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
  contours->Update();

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

  // Export contour lines to GeoJSON
  std::string contourLinesGeoJSON = "contour_lines.geojson";
  ExportContoursToGeoJSON(smoothContours->GetOutput(), contourLinesGeoJSON, scalarRange[0], scalarRange[1]);
  
  // Export filled contours to GeoJSON
  std::string filledContoursGeoJSON = "filled_contours.geojson";
  ExportFilledContoursToGeoJSON(filledContours->GetOutput(), filledContoursGeoJSON, scalarRange[0], scalarRange[1]);

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
