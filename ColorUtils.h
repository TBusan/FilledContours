#ifndef COLOR_UTILS_H
#define COLOR_UTILS_H

#include <string>

namespace ColorUtils {

// Helper function for hex character conversion
inline int hexCharToInt(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return 0; // Default for invalid characters
}

// Function to convert hex color to RGB
inline void hexToRGB(const std::string& hexStr, double rgb[3]) {
    // Remove # if present
    std::string hex = hexStr;
    if (hex.length() > 0 && hex[0] == '#') {
        hex = hex.substr(1);
    }
    
    // Ensure we have at least 6 characters
    if (hex.length() < 6) {
        rgb[0] = rgb[1] = rgb[2] = 0.0;
        return;
    }
    
    // Convert hex to RGB manually without using std::hex
    int r = (hexCharToInt(hex[0]) << 4) + hexCharToInt(hex[1]);
    int g = (hexCharToInt(hex[2]) << 4) + hexCharToInt(hex[3]);
    int b = (hexCharToInt(hex[4]) << 4) + hexCharToInt(hex[5]);
    
    // Normalize to 0-1 range
    rgb[0] = r / 255.0;
    rgb[1] = g / 255.0;
    rgb[2] = b / 255.0;
}

} // namespace ColorUtils

#endif // COLOR_UTILS_H 