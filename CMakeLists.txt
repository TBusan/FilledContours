cmake_minimum_required(VERSION 3.12)

project(FilledContours)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)  # Use a higher C++ standard for better compatibility
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for emscripten
if(EMSCRIPTEN)
  message(STATUS "Building with Emscripten for WebAssembly")
  
  # Emscripten specific compiler and linker flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s ALLOW_MEMORY_GROWTH=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap'] -s EXPORT_NAME=FilledContoursModule -s MODULARIZE=1")
  
  # If we want to enable threading (optional)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -sPROXY_TO_PTHREAD=1 -sOFFSCREENCANVAS_SUPPORT=1")
endif()

# Find VTK
find_package(VTK REQUIRED)
message(STATUS "Found VTK version: ${VTK_VERSION}")

# Optionally find standalone jsoncpp if not bundled with VTK
find_package(jsoncpp QUIET)
if(NOT jsoncpp_FOUND)
  message(STATUS "Using VTK bundled jsoncpp")
endif()

# Include directories
include_directories(${VTK_INCLUDE_DIRS})

# Ensure all VTK system headers are found
if(MSVC)
  # For Visual Studio, add system include directories
  include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
  add_compile_options(/EHsc /MP /W3 /utf-8)
else()
  # For gcc/clang, mark VTK includes as system headers
  include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
endif()

# Headers
set(HEADERS
    ColorUtils.h
)

# Sources
set(SOURCES 
    FilledContours1.cxx
)

# Create the executable/library
if(EMSCRIPTEN)
  # For WebAssembly, create a JavaScript module
  add_executable(FilledContours1 ${SOURCES} ${HEADERS})
  set_target_properties(FilledContours1 PROPERTIES SUFFIX ".js")
else()
  # For native builds, create a executable
  add_executable(FilledContours1 ${SOURCES} ${HEADERS})
endif()

# Link with VTK libraries
target_link_libraries(FilledContours1 ${VTK_LIBRARIES})

# VTK module initialization
if(VTK_VERSION VERSION_GREATER_EQUAL "8.90.0")
  vtk_module_autoinit(
    TARGETS FilledContours1
    MODULES ${VTK_LIBRARIES}
  )
endif()

# Installation
install(TARGETS FilledContours1
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
