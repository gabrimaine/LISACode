project(USO)
INCLUDE_DIRECTORIES(Include ../Noise/Include/ ../ToolBox/Include/)
#file(GLOB_RECURSE USO_SOURCES "Src/*.c*")
add_library (lisacode2uso Src/LISACODE-USO.cpp Src/LISACODE-USONoiseDeriv.cpp)
target_link_libraries(lisacode2uso lisacode2noise lisacode2toolbox)

add_executable(LC2TestUSO "Test/LISACODE-TestUSO.cpp")
target_link_libraries(LC2TestUSO lisacode2uso lisacode2noise lisacode2toolbox  ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
