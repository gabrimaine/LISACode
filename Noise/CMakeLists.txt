project(Noise)
INCLUDE_DIRECTORIES(Include ../ToolBox/Include/)
add_library (lisacode2noise	Src/LISACODE-Noise.cpp 	Src/LISACODE-NoiseWhite.cpp 	Src/LISACODE-NoiseFilter.cpp    Src/LISACODE-NoiseFile.cpp)
target_link_libraries(lisacode2noise lisacode2toolbox)

add_executable(LC2TestNoise "Test/LISACODE-TestNoises.cpp")
target_link_libraries(LC2TestNoise lisacode2noise lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
