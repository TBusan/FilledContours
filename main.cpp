#include "FilledContours.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <data_json_file> <color_scale_json_file>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string dataJsonPath = argv[1];
    std::string colorScalePath = argv[2];
    
    std::cout << "Processing contours with:" << std::endl;
    std::cout << "  Data JSON: " << dataJsonPath << std::endl;
    std::cout << "  Color Scale: " << colorScalePath << std::endl;

    int result = processContours(dataJsonPath, colorScalePath);
    
    if (result == EXIT_SUCCESS) {
        std::cout << "Successfully processed contours and saved GeoJSON files." << std::endl;
    } else {
        std::cerr << "Failed to process contours." << std::endl;
    }
    
    return result;
} 