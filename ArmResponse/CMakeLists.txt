project(ArmResponse)
INCLUDE_DIRECTORIES(Include ../Orbits/Include/ ../GW/Include/ ../Noise/Include/ ../ToolBox/Include/)
add_library (lisacode2armresp Src/LISACODE-ArmResp.cpp)
target_link_libraries(lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox)

add_executable(LC2GWArmResp "Test/LISACODE-TestArmResp.cpp")
target_link_libraries(LC2GWArmResp lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
