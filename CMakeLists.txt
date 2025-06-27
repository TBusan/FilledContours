cmake_minimum_required(VERSION 3.12)

project(FilledContours)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
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

# Find jsoncpp (usually included with VTK)
find_package(jsoncpp QUIET)
if(NOT jsoncpp_FOUND)
  message(STATUS "Using VTK bundled jsoncpp")
endif()

# Include directories
include_directories(${VTK_INCLUDE_DIRS})

# Sources
set(SOURCES 
    FilledContours.cxx
)

set(HEADERS
    FilledContours.h
)

# Create the library
if(EMSCRIPTEN)
  # For WebAssembly, create a JavaScript module
  add_executable(FilledContours ${SOURCES} ${HEADERS})
  set_target_properties(FilledContours PROPERTIES SUFFIX ".js")
else()
  # For native builds, create a library
  add_library(FilledContours SHARED ${SOURCES} ${HEADERS})
endif()

# Link with VTK libraries
target_link_libraries(FilledContours ${VTK_LIBRARIES})

# For non-emscripten builds, create a test executable
if(NOT EMSCRIPTEN)
  add_executable(TestFilledContours main.cpp)
  target_link_libraries(TestFilledContours FilledContours ${VTK_LIBRARIES})
endif()

# Installation
install(TARGETS FilledContours
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

if(NOT EMSCRIPTEN)
    install(FILES ${HEADERS} DESTINATION include)
endif()
