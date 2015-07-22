project(LISACode)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package (GSL REQUIRED)
find_package (FFTW3 REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS})

add_subdirectory(ToolBox)
add_subdirectory(Noise)
add_subdirectory(GW)
add_subdirectory(Orbits)
add_subdirectory(ArmResponse)
add_subdirectory(USO)
add_subdirectory(OptBenPhaMet)
add_subdirectory(Detector)
add_subdirectory(TDI)
add_subdirectory(Main)
#add_subdirectory(Documentation/ManualTex)

