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
#include <sstream>
#include <iomanip>

// This program loads JSON data with null values and creates filled contours visualization
// The null value boundaries are used as masks for contour lines and filled contours
// 原始的demo 加载json文件 然后导出成 geojson 格式的数据    新增输入的参数 colorScale的参数 用于设置颜色  (修改了colorScale的参数的格式 ：格式如colorScale_1.json)

// Function to convert hex color to RGB
void hexToRGB(const std::string& hexColor, double rgb[3]) {
    // Remove # if present
    std::string hex = hexColor;
    if (hex[0] == '#') {
        hex = hex.substr(1);
    }
    
    // Convert hex to RGB
    int r, g, b;
    std::stringstream ss;
    ss << std::hex << hex.substr(0, 2);
    ss >> r;
    ss.clear();
    ss << std::hex << hex.substr(2, 2);
    ss >> g;
    ss.clear();
    ss << std::hex << hex.substr(4, 2);
    ss >> b;
    
    // Normalize to 0-1 range
    rgb[0] = r / 255.0;
    rgb[1] = g / 255.0;
    rgb[2] = b / 255.0;
}

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
  std::cout << "Program started. Arguments count: " << argc << std::endl;
  
  if (argc < 3)
  {
    std::cerr
        << "Usage: " << argv[0]
        << " InputJSONFile(.json) ColorScaleFile(.json) e.g testData.json colorScale.json"
        << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input JSON file: " << argv[1] << std::endl;
  std::cout << "Color scale file: " << argv[2] << std::endl;

  vtkNew<vtkNamedColors> colors;

  // Read the JSON data file
  std::cout << "Attempting to open data JSON file..." << std::endl;
  std::ifstream jsonFile(argv[1]);
  if (!jsonFile.is_open()) {
    std::cerr << "Failed to open JSON file: " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Data JSON file opened successfully" << std::endl;

  Json::Value root;
  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  std::cout << "Parsing data JSON..." << std::endl;
  if (!Json::parseFromStream(builder, jsonFile, &root, &errs)) {
    std::cerr << "Error parsing JSON: " << errs << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Data JSON parsed successfully" << std::endl;

  // Read the color scale JSON file
  std::cout << "Attempting to open color scale JSON file..." << std::endl;
  std::ifstream colorScaleFile(argv[2]);
  if (!colorScaleFile.is_open()) {
    std::cerr << "Failed to open color scale JSON file: " << argv[2] << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Color scale JSON file opened successfully" << std::endl;

  Json::Value colorScaleRoot;
  std::cout << "Parsing color scale JSON..." << std::endl;
  if (!Json::parseFromStream(builder, colorScaleFile, &colorScaleRoot, &errs)) {
    std::cerr << "Error parsing color scale JSON: " << errs << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Color scale JSON parsed successfully" << std::endl;

  // Extract color scale data - colorScale_1.json is a direct array, not wrapped in a "colorScale" field
  Json::Value& colorScale = colorScaleRoot;
  if (colorScale.isNull() || !colorScale.isArray()) {
    std::cerr << "Invalid color scale format. Expected an array of color scale entries." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Found color scale with " << colorScale.size() << " entries" << std::endl;

  // Extract data from JSON
  std::cout << "Extracting data from JSON..." << std::endl;
  Json::Value& data = root["data"];
  if (data.isNull()) {
    std::cerr << "Error: No 'data' field found in JSON" << std::endl;
    return EXIT_FAILURE;
  }
  
  Json::Value& x = data["x"];
  if (x.isNull() || !x.isArray()) {
    std::cerr << "Error: Invalid or missing 'x' array in data" << std::endl;
    return EXIT_FAILURE;
  }
  
  Json::Value& y = data["y"];
  if (y.isNull() || !y.isArray()) {
    std::cerr << "Error: Invalid or missing 'y' array in data" << std::endl;
    return EXIT_FAILURE;
  }
  
  Json::Value& v = data["v"];
  if (v.isNull() || !v.isArray()) {
    std::cerr << "Error: Invalid or missing 'v' array in data" << std::endl;
    return EXIT_FAILURE;
  }
  
  double zmin = data.get("zmin", 0.0).asDouble();
  double zmax = data.get("zmax", 1.0).asDouble();
  std::cout << "Data extracted: x size=" << x.size() << ", y size=" << y.size() 
            << ", zmin=" << zmin << ", zmax=" << zmax << std::endl;

  // Create points and values
  std::cout << "Creating points and values..." << std::endl;
  vtkNew<vtkPoints> points;
  vtkNew<vtkFloatArray> scalars;
  scalars->SetName("Elevation");

  // Create a mask for valid points
  std::vector<bool> validPoint;
  int validPointCount = 0;

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
        validPointCount++;
        // Add point and scalar value
        points->InsertNextPoint(xCoord, yCoord, 0.0);
        scalars->InsertNextValue(value);
      }
    }
  }
  std::cout << "Created " << validPointCount << " valid points out of " << x.size() * y.size() << " total points" << std::endl;

  // Create polydata
  std::cout << "Creating polydata..." << std::endl;
  vtkNew<vtkPolyData> polydata;
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(scalars);

  // Use Delaunay triangulation to create a surface from points
  std::cout << "Performing Delaunay triangulation..." << std::endl;
  vtkNew<vtkDelaunay2D> delaunay;
  delaunay->SetInputData(polydata);
  delaunay->SetTolerance(0.001);
  delaunay->Update();
  std::cout << "Triangulation complete. Number of cells: " << delaunay->GetOutput()->GetNumberOfCells() << std::endl;

  // Create contour levels from colorScale
  std::cout << "Creating contour levels..." << std::endl;
  int numberOfContours = colorScale.size();
  std::vector<double> contourLevels;
  
  for (int i = 0; i < numberOfContours; i++) {
    contourLevels.push_back(colorScale[i]["level"].asDouble());
    std::cout << "  Level " << i << ": " << contourLevels[i] << std::endl;
  }

  // Clip the triangulation to remove triangles in null regions
  std::cout << "Creating filled contours..." << std::endl;
  vtkNew<vtkAppendPolyData> appendFilledContours;

  // Keep the clippers alive
  std::vector<vtkSmartPointer<vtkClipPolyData>> clippers;
  
  // Instead of clipping between adjacent contour levels, handle the entire data range
  // First, handle data below the first contour level
  if (!contourLevels.empty()) {
    // Handle values below the first contour level
    vtkNew<vtkClipPolyData> clipLow;
    clipLow->SetValue(contourLevels[0]);
    clipLow->SetInputConnection(delaunay->GetOutputPort());
    clipLow->InsideOutOn(); // Get points BELOW the first level
    clipLow->Update();
    
    if (clipLow->GetOutput()->GetNumberOfCells() > 0) {
      std::cout << "  Processing values below first contour level (" << contourLevels[0] << ")" << std::endl;
      
      vtkNew<vtkFloatArray> cd;
      cd->SetNumberOfComponents(1);
      cd->SetNumberOfTuples(clipLow->GetOutput()->GetNumberOfCells());
      cd->FillComponent(0, contourLevels[0] - 1.0); // Use a value lower than first level
      
      clipLow->GetOutput()->GetCellData()->SetScalars(cd);
      appendFilledContours->AddInputConnection(clipLow->GetOutputPort());
      std::cout << "    Added band with " << clipLow->GetOutput()->GetNumberOfCells() << " cells" << std::endl;
    }
    
    // Now handle each interval between contour levels
    for (int i = 0; i < numberOfContours - 1; i++) {
      double valueLo = contourLevels[i];
      double valueHi = contourLevels[i + 1];
      
      std::cout << "  Processing contour band " << i << " [" << valueLo << " - " << valueHi << "]" << std::endl;
      
      vtkNew<vtkClipPolyData> clipLo;
      clipLo->SetValue(valueLo);
      clipLo->SetInputConnection(delaunay->GetOutputPort());
      clipLo->InsideOutOff(); // Get points ABOVE the lower level
      clipLo->Update();
      
      vtkNew<vtkClipPolyData> clipHi;
      clipHi->SetValue(valueHi);
      clipHi->SetInputConnection(clipLo->GetOutputPort());
      clipHi->InsideOutOn(); // Get points BELOW the higher level
      clipHi->Update();
      
      if (clipHi->GetOutput()->GetNumberOfCells() == 0) {
        std::cout << "    No cells in this contour band, skipping" << std::endl;
        continue;
      }
      
      vtkNew<vtkFloatArray> cd;
      cd->SetNumberOfComponents(1);
      cd->SetNumberOfTuples(clipHi->GetOutput()->GetNumberOfCells());
      cd->FillComponent(0, valueHi); // Use the higher level value for this band
      
      clipHi->GetOutput()->GetCellData()->SetScalars(cd);
      appendFilledContours->AddInputConnection(clipHi->GetOutputPort());
      clippers.push_back(clipHi); // Keep the clipper alive
      std::cout << "    Added band with " << clipHi->GetOutput()->GetNumberOfCells() << " cells" << std::endl;
    }
    
    // Finally, handle values above the highest contour level
    vtkNew<vtkClipPolyData> clipHigh;
    clipHigh->SetValue(contourLevels.back());
    clipHigh->SetInputConnection(delaunay->GetOutputPort());
    clipHigh->InsideOutOff(); // Get points ABOVE the highest level
    clipHigh->Update();
    
    if (clipHigh->GetOutput()->GetNumberOfCells() > 0) {
      std::cout << "  Processing values above last contour level (" << contourLevels.back() << ")" << std::endl;
      
      vtkNew<vtkFloatArray> cd;
      cd->SetNumberOfComponents(1);
      cd->SetNumberOfTuples(clipHigh->GetOutput()->GetNumberOfCells());
      cd->FillComponent(0, contourLevels.back() + 1.0); // Use a value higher than last level
      
      clipHigh->GetOutput()->GetCellData()->SetScalars(cd);
      appendFilledContours->AddInputConnection(clipHigh->GetOutputPort());
      std::cout << "    Added band with " << clipHigh->GetOutput()->GetNumberOfCells() << " cells" << std::endl;
    }
  }

  // Clean the data
  std::cout << "Cleaning polydata..." << std::endl;
  vtkNew<vtkCleanPolyData> filledContours;
  filledContours->SetInputConnection(appendFilledContours->GetOutputPort());
  filledContours->Update();
  std::cout << "Filled contours created with " << filledContours->GetOutput()->GetNumberOfCells() << " cells" << std::endl;

  // Create a lookup table for coloring based on colorScale
  std::cout << "Creating color lookup table..." << std::endl;
  vtkNew<vtkLookupTable> lut;
  lut->SetNumberOfTableValues(numberOfContours + 2); // +2 for below lowest and above highest
  lut->SetTableRange(contourLevels.front() - 2.0, contourLevels.back() + 2.0); // Extend range
  lut->Build();

  // Set custom color for each contour band from colorScale
  // First value: for data below the lowest level
  std::string firstHexColor = colorScale[0]["color"].asString();
  double firstRgb[3];
  hexToRGB(firstHexColor, firstRgb);
  lut->SetTableValue(0, firstRgb[0], firstRgb[1], firstRgb[2]);
  std::cout << "  Color for values below " << contourLevels.front() << ": " << firstHexColor << std::endl;
  
  // Middle values: for each contour level
  for (int i = 0; i < numberOfContours; i++) {
    std::string hexColor = colorScale[i]["color"].asString();
    double rgb[3];
    hexToRGB(hexColor, rgb);
    lut->SetTableValue(i + 1, rgb[0], rgb[1], rgb[2]);
    std::cout << "  Color for level " << contourLevels[i] << ": " << hexColor << " -> RGB(" << rgb[0] << "," << rgb[1] << "," << rgb[2] << ")" << std::endl;
  }
  
  // Last value: for data above the highest level
  std::string lastHexColor = colorScale[numberOfContours - 1]["color"].asString();
  double lastRgb[3];
  hexToRGB(lastHexColor, lastRgb);
  lut->SetTableValue(numberOfContours + 1, lastRgb[0], lastRgb[1], lastRgb[2]);
  std::cout << "  Color for values above " << contourLevels.back() << ": " << lastHexColor << std::endl;

  // Create the mapper for filled contours
  std::cout << "Creating filled contour mapper..." << std::endl;
  vtkNew<vtkPolyDataMapper> contourMapper;
  contourMapper->SetInputConnection(filledContours->GetOutputPort());
  contourMapper->SetScalarRange(contourLevels.front() - 2.0, contourLevels.back() + 2.0); // Extend range
  contourMapper->SetScalarModeToUseCellData();
  contourMapper->SetLookupTable(lut);

  vtkNew<vtkActor> contourActor;
  contourActor->SetMapper(contourMapper);
  contourActor->GetProperty()->SetInterpolationToFlat();

  // Create contour lines
  std::cout << "Creating contour lines..." << std::endl;
  vtkNew<vtkContourFilter> contours;
  contours->SetInputConnection(delaunay->GetOutputPort()); // Use original data to ensure complete contour lines
  
  // Set contour values from colorScale
  for (int i = 0; i < numberOfContours; i++) {
    contours->SetValue(i, contourLevels[i]);
  }
  contours->Update();
  std::cout << "Contour lines created with " << contours->GetOutput()->GetNumberOfCells() << " cells" << std::endl;

  // Smooth the contours
  std::cout << "Smoothing contour lines..." << std::endl;
  vtkNew<vtkSmoothPolyDataFilter> smoothContours;
  smoothContours->SetInputConnection(contours->GetOutputPort());
  smoothContours->SetNumberOfIterations(2);
  smoothContours->SetRelaxationFactor(0.2);
  try {
    smoothContours->Update();
    std::cout << "Smoothing complete. Output has " << smoothContours->GetOutput()->GetNumberOfCells() << " cells" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception during smoothing: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception during smoothing" << std::endl;
  }

  std::cout << "Creating contour line mapper..." << std::endl;
  vtkNew<vtkPolyDataMapper> contourLineMapperer;
  try {
    contourLineMapperer->SetInputConnection(smoothContours->GetOutputPort());
    contourLineMapperer->SetScalarRange(contourLevels.front(), contourLevels.back());
    contourLineMapperer->ScalarVisibilityOff();
    std::cout << "Contour line mapper created successfully" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception creating contour line mapper: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception creating contour line mapper" << std::endl;
  }

  std::cout << "Creating contour line actor..." << std::endl;
  vtkNew<vtkActor> contourLineActor;
  try {
    contourLineActor->SetMapper(contourLineMapperer);
    contourLineActor->GetProperty()->SetLineWidth(2);
    
    // Use lineColor from colorScale for contour lines
    std::string lineHexColor = colorScale[0]["lineColor"].asString();
    double lineRgb[3];
    hexToRGB(lineHexColor, lineRgb);
    contourLineActor->GetProperty()->SetColor(lineRgb);
    std::cout << "Contour line color: " << lineHexColor << " -> RGB(" << lineRgb[0] << "," << lineRgb[1] << "," << lineRgb[2] << ")" << std::endl;
    std::cout << "Contour line actor created successfully" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception creating contour line actor: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception creating contour line actor" << std::endl;
  }

  // Create boundary edges around valid data regions
  std::cout << "Creating boundary edges..." << std::endl;
  vtkNew<vtkFeatureEdges> boundaryEdges;
  try {
    // Use the filtered points (remove NULL_VALUE points) to create boundary
    boundaryEdges->SetInputConnection(filledContours->GetOutputPort());
    boundaryEdges->BoundaryEdgesOn();
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
    boundaryEdges->Update();
    std::cout << "Boundary edges created with " << boundaryEdges->GetOutput()->GetNumberOfCells() << " cells" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception creating boundary edges: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception creating boundary edges" << std::endl;
  }

  std::cout << "Creating boundary mapper..." << std::endl;
  vtkNew<vtkPolyDataMapper> boundaryMapper;
  try {
    boundaryMapper->SetInputConnection(boundaryEdges->GetOutputPort());
    std::cout << "Boundary mapper created successfully" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception creating boundary mapper: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception creating boundary mapper" << std::endl;
  }
  
  std::cout << "Creating boundary actor..." << std::endl;
  vtkNew<vtkActor> boundaryActor;
  try {
    boundaryActor->SetMapper(boundaryMapper);
    boundaryActor->GetProperty()->SetLineWidth(3);
    boundaryActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
    std::cout << "Boundary actor created successfully" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception creating boundary actor: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception creating boundary actor" << std::endl;
  }

  // Export contour lines to GeoJSON
  std::cout << "Exporting contour lines to GeoJSON..." << std::endl;
  try {
    std::string contourLinesGeoJSON = "contour_lines.geojson";
    ExportContoursToGeoJSON(smoothContours->GetOutput(), contourLinesGeoJSON, zmin, zmax);
    std::cout << "Contour lines export complete" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception exporting contour lines: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception exporting contour lines" << std::endl;
  }
  
  // Export filled contours to GeoJSON
  std::cout << "Exporting filled contours to GeoJSON..." << std::endl;
  try {
    std::string filledContoursGeoJSON = "filled_contours.geojson";
    ExportFilledContoursToGeoJSON(filledContours->GetOutput(), filledContoursGeoJSON, zmin, zmax);
    std::cout << "Filled contours export complete" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception exporting filled contours: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception exporting filled contours" << std::endl;
  }

  // The usual renderer, render window and interactor
  std::cout << "Setting up renderer and window..." << std::endl;
  vtkNew<vtkRenderer> ren1;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  try {
    renWin->AddRenderer(ren1);
    renWin->SetWindowName("FilledContours");
    // Set window size to make it more visible
    renWin->SetSize(800, 600);
    // 强制窗口在前台显示
    renWin->SetFullScreen(0);
    renWin->SetBorders(1);
    renWin->SetPosition(100, 100);
    renWin->SetOffScreenRendering(0);
    std::cout << "Render window setup complete" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception setting up render window: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception setting up render window" << std::endl;
  }

  try {
    iren->SetRenderWindow(renWin);
    // 设置交互器样式
    iren->Initialize();
    std::cout << "Interactor setup complete" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception setting up interactor: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception setting up interactor" << std::endl;
  }

  // Add the actors
  std::cout << "Adding actors to renderer..." << std::endl;
  try {
    ren1->AddActor(contourActor);
    ren1->AddActor(contourLineActor);
    ren1->AddActor(boundaryActor);
    ren1->SetBackground(colors->GetColor3d("MidnightBlue").GetData());
    std::cout << "Actors added successfully" << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception adding actors: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception adding actors" << std::endl;
  }

  // Begin interaction
  std::cout << "Starting rendering..." << std::endl;
  std::cout << "IMPORTANT: If no window appears, check if it's behind other windows or minimized." << std::endl;
  
  try {
    std::cout << "Calling renWin->Render()..." << std::endl;
    renWin->Render();
    std::cout << "Render complete. Starting interaction..." << std::endl;
    
    std::cout << "Calling iren->Start()..." << std::endl;
    iren->Start();
    std::cout << "Interaction complete. Program ending." << std::endl;
  }
  catch (std::exception& e) {
    std::cerr << "Exception during rendering or interaction: " << e.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception during rendering or interaction" << std::endl;
  }

  return EXIT_SUCCESS;
}
