#ifndef FILLED_CONTOURS_H
#define FILLED_CONTOURS_H

#include <string>

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif

/**
 * Process data from JSON files and generate contour and filled contour GeoJSON files
 * 
 * @param dataJsonPath Path to the input JSON data file with x, y, v arrays
 * @param colorScalePath Path to the color scale JSON file
 * @param contourLinesOutputPath Path where contour lines GeoJSON will be saved
 * @param filledContoursOutputPath Path where filled contours GeoJSON will be saved
 * @return 0 on success, non-zero on failure
 */
int processContours(
    const std::string& dataJsonPath,
    const std::string& colorScalePath,
    const std::string& contourLinesOutputPath = "contour_lines.geojson",
    const std::string& filledContoursOutputPath = "filled_contours.geojson"
);

/**
 * Version that takes JSON content as strings instead of files
 * 
 * @param dataJsonContent JSON data content string with x, y, v arrays
 * @param colorScaleContent Color scale JSON content string
 * @param contourLinesOutputPath Path where contour lines GeoJSON will be saved
 * @param filledContoursOutputPath Path where filled contours GeoJSON will be saved
 * @return 0 on success, non-zero on failure
 */
int processContoursFromStrings(
    const std::string& dataJsonContent,
    const std::string& colorScaleContent,
    const std::string& contourLinesOutputPath = "contour_lines.geojson",
    const std::string& filledContoursOutputPath = "filled_contours.geojson"
);

/**
 * Version that takes JSON content as strings and returns GeoJSON results as strings
 * 
 * @param dataJsonContent JSON data content string with x, y, v arrays
 * @param colorScaleContent Color scale JSON content string
 * @param contourLinesOutput String to store the contour lines GeoJSON output
 * @param filledContoursOutput String to store the filled contours GeoJSON output
 * @return 0 on success, non-zero on failure
 */
int processContoursToStrings(
    const std::string& dataJsonContent,
    const std::string& colorScaleContent,
    std::string& contourLinesOutput,
    std::string& filledContoursOutput
);

#endif // FILLED_CONTOURS_H 