# FilledContours WebAssembly Setup

This document summarizes the setup created for compiling FilledContours1.cxx to WebAssembly.

## Files Created

1. **CMake Configuration**
   - `CMakeLists_wasm.txt`: CMake configuration for WebAssembly compilation
   
2. **HTML and Web UI**
   - `shell.html`: Template HTML file used by Emscripten for the compiled output
   - `web/index.html`: Landing page that loads the WebAssembly application
   - `web/style.css`: CSS styles for the web interface
   - `web/filledcontours-api.js`: JavaScript API wrapper for the WebAssembly module
   - `web/example.html`: Example usage of the JavaScript API

3. **Build and Run Scripts**
   - `build_wasm.sh`: Shell script for building WebAssembly on Unix-like systems
   - `build_wasm.bat`: Batch script for building WebAssembly on Windows
   - `run_server.py`: Python script to run a local web server for testing

4. **Documentation**
   - `README_WASM.md`: Detailed instructions for building and using the WebAssembly version
   - `WASM_README.md`: This summary file

## Build Process Overview

1. **Prerequisites**:
   - Emscripten SDK
   - VTK compiled with Emscripten
   - CMake 3.15 or higher

2. **Building**:
   - Update `VTK_DIR` in `CMakeLists_wasm.txt` to point to your Emscripten-compiled VTK
   - Run `build_wasm.sh` (Unix) or `build_wasm.bat` (Windows) to build the WebAssembly version
   - Output files are placed in the `build-wasm/bin/` directory

3. **Testing**:
   - Run `python run_server.py` to start a local web server
   - Open a browser and navigate to http://localhost:8000

## WebAssembly Output Files

After compilation, the following key files are produced:

- `FilledContours_wasm.js`: JavaScript glue code
- `FilledContours_wasm.wasm`: WebAssembly binary module
- `FilledContours_wasm.data`: Preloaded files (JSON data)
- `FilledContours_wasm.html`: HTML wrapper

## Integration with Web Applications

The FilledContours WebAssembly module can be integrated into web applications using the provided JavaScript API:

```javascript
// Initialize the API
const filledContours = new FilledContoursAPI({
  moduleUrl: 'FilledContours_wasm.js',
  onReady: () => {
    console.log('Module ready!');
  }
});

// Run visualization
filledContours.runVisualization('testData.json', 'colorScale.json');
```

See `web/example.html` for a complete usage example. 