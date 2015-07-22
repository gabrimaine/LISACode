project(GW)
INCLUDE_DIRECTORIES(Include ../Noise/Include/ ../ToolBox/Include/)
add_library (lisacode2gw 	Src/LISACODE-GW.cpp Src/LISACODE-GWFile.cpp Src/LISACODE-GWStochastic.cpp Src/LISACODE-GWGalBin.cpp Src/LISACODE-GWSpinBBHHHarm1.cpp Src/LISACODE-GWSpinBBHNR1.cpp Src/LISACODE-GWCosmicString.cpp)
target_link_libraries(lisacode2gw lisacode2noise lisacode2toolbox )

add_executable(LC2GWbary "Test/LISACODE-TestGWs.cpp")
target_link_libraries(LC2GWbary lisacode2gw lisacode2noise lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
