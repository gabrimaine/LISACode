project(Detector)
INCLUDE_DIRECTORIES(Include ../OptBenPhaMet/Include ../USO/Include ../ArmResponse/Include ../Orbits/Include ../GW/Include ../Noise/Include ../ToolBox/Include)

add_library (lisacode2detector 	Src/LISACODE-Detector.cpp)
target_link_libraries(lisacode2detector lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox)

add_executable(LC2TestDetector "Test/LISACODE-TestDetector.cpp")
target_link_libraries(LC2TestDetector lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
