project(OptBenPhaMet)
INCLUDE_DIRECTORIES(Include ../USO/Include ../ArmResponse/Include ../Orbits/Include ../GW/Include ../Noise/Include ../ToolBox/Include)
add_library (lisacode2obpm Src/LISACODE-OBPM.cpp Src/LISACODE-OBPMStd2002.cpp	Src/LISACODE-OBPMStd2010.cpp Src/LISACODE-OBPMOcta1.cpp)
target_link_libraries(lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox)

add_executable(LC2TestOBPM "Test/LISACODE-TestOBPM.cpp")
target_link_libraries(LC2TestOBPM lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
