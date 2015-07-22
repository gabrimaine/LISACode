project(TDI)
INCLUDE_DIRECTORIES(Include  ../Detector/Include ../OptBenPhaMet/Include ../USO/Include ../ArmResponse/Include ../Orbits/Include/ ../GW/Include/ ../Noise/Include/ ../ToolBox/Include/)
add_library (lisacode2tdi Src/LISACODE-TDIInt.cpp Src/LISACODE-TDIIntStd2002.cpp Src/LISACODE-TDIGen.cpp)
target_link_libraries(lisacode2tdi lisacode2detector lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox)

add_executable(LC2TestTDI "Test/LISACODE-TestTDI.cpp")
target_link_libraries(LC2TestTDI lisacode2tdi lisacode2detector lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
