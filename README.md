# Filled Contours WebAssembly Module

This project converts a VTK-based filled contours generator to WebAssembly for use in web applications. It takes JSON data with x, y, v arrays and a color scale JSON to generate contour lines and filled contours as GeoJSON.

## Prerequisites

To build this project, you need:

1. [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html) (version 3.0.0 or later)
2. [CMake](https://cmake.org/download/) (version 3.12 or later)
3. [VTK](https://vtk.org/download/) built with Emscripten (see instructions below)

## Building VTK for WebAssembly

First, you need to build VTK for WebAssembly using Emscripten:

```bash
# Clone VTK repository
git clone https://github.com/Kitware/VTK.git
cd VTK
mkdir build-wasm
cd build-wasm

# Configure VTK for WebAssembly
emcmake cmake .. \
  -G "Ninja" \
  -DBUILD_SHARED_LIBS:BOOL=OFF \
  -DVTK_ENABLE_LOGGING:BOOL=OFF \
  -DVTK_WRAP_JAVASCRIPT:BOOL=ON \
  -DVTK_MODULE_ENABLE_VTK_hdf5:STRING=NO \
  -DVTK_MODULE_ENABLE_VTK_RenderingContextOpenGL2:STRING=DONT_WANT \
  -DVTK_MODULE_ENABLE_VTK_RenderingCellGrid:STRING=NO \
  -DVTK_MODULE_ENABLE_VTK_sqlite:STRING=NO

# Build VTK
ninja

# Install VTK (optional, use a destination directory if needed)
ninja install
```

## Building the Filled Contours Module

After building VTK for WebAssembly, you can build this project:

```bash
mkdir build-wasm
cd build-wasm

# Configure the project
emcmake cmake .. \
  -G "Ninja" \
  -DVTK_DIR=/path/to/vtk/build-wasm

# Build the project
ninja
```

The build process will generate:
- `FilledContours.js`: The JavaScript loader for the WebAssembly module
- `FilledContours.wasm`: The WebAssembly binary

## Using the Module in a Web Application

1. Copy the `FilledContours.js` and `FilledContours.wasm` files to your web project
2. Include the JavaScript file in your HTML:
   ```html
   <script src="FilledContours.js"></script>
   ```
3. Load and use the module in your JavaScript:
   ```javascript
   // Initialize the module
   const wasmModule = await FilledContoursModule();
   
   // Call the processContoursToStrings function
   let contourLinesOutput = '';
   let filledContoursOutput = '';
   
   const result = wasmModule.processContoursToStrings(
     dataJsonContent,  // Your data JSON string
     colorScaleContent, // Your color scale JSON string
     contourLinesOutput,
     filledContoursOutput
   );
   
   if (result === 0) {
     console.log("Contour lines:", contourLinesOutput);
     console.log("Filled contours:", filledContoursOutput);
   }
   ```

## Demo Application

A demo HTML application is included in this project:
1. Copy `index.html`, `FilledContours.js`, and `FilledContours.wasm` to a web server
2. Open `index.html` in a browser
3. Input your data and color scale JSON
4. Click "Process Contours" to generate the GeoJSON outputs

## Input Data Format

### Data JSON:
```json
{
  "data": {
    "x": [0, 1, 2, ...],           // X coordinates array
    "y": [0, 1, 2, ...],           // Y coordinates array
    "v": [                         // 2D array of values at each (x,y) point
      [0.1, 0.2, 0.3, ...],        // null values are allowed for missing data
      [0.2, null, 0.4, ...],
      ...
    ],
    "zmin": 0.0,                   // Minimum value range (optional)
    "zmax": 1.0                    // Maximum value range (optional)
  }
}
```

### Color Scale JSON:
```json
[
  {"level": 0.2, "color": "#0000FF", "lineColor": "#000000"},
  {"level": 0.4, "color": "#00FFFF", "lineColor": "#000000"},
  {"level": 0.6, "color": "#00FF00", "lineColor": "#000000"},
  {"level": 0.8, "color": "#FFFF00", "lineColor": "#000000"},
  {"level": 1.0, "color": "#FF0000", "lineColor": "#000000"}
]
```

## Output Format

The module outputs GeoJSON for both contour lines and filled contours:

### Contour Lines GeoJSON:
- Format: MultiLineString features
- Each feature represents a specific contour level
- Properties include contourValue and normalizedValue

### Filled Contours GeoJSON:
- Format: MultiPolygon features
- Each feature represents a band between contour levels
- Properties include contourValue and normalizedValue

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This project uses:
- [VTK](https://vtk.org/) - The Visualization Toolkit
- [Emscripten](https://emscripten.org/) - WebAssembly compiler toolchain
- [jsoncpp](https://github.com/open-source-parsers/jsoncpp) - JSON parser for C++ 