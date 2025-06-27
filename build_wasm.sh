#!/bin/bash
# Script to build FilledContours as WebAssembly

# Check if Emscripten is activated
if [ -z "$EMSCRIPTEN" ]; then
    echo "Error: Emscripten environment not detected"
    echo "Please activate the Emscripten environment first by running:"
    echo "source /path/to/emsdk/emsdk_env.sh"
    exit 1
fi

# Create build directory if it doesn't exist
if [ ! -d "build-wasm" ]; then
    mkdir build-wasm
fi

# Enter build directory
cd build-wasm

# Configure with CMake
echo "Configuring with CMake..."
emcmake cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=$EMSCRIPTEN/cmake/Modules/Platform/Emscripten.cmake \
    -DCMAKE_CROSSCOMPILING_EMULATOR=node \
    -DCMAKE_INSTALL_PREFIX=./install \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -C ../CMakeLists_wasm.txt

if [ $? -ne 0 ]; then
    echo "CMake configuration failed"
    exit 1
fi

# Build
echo "Building..."
emmake make -j4

if [ $? -ne 0 ]; then
    echo "Build failed"
    exit 1
fi

echo "Build successful!"
echo "Output files are in build-wasm/bin/"
echo "To run a local web server:"
echo "cd build-wasm/bin && python -m http.server 8000"
echo "Then navigate to http://localhost:8000/ in your browser" 