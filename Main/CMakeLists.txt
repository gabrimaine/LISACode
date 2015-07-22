project(Main)
INCLUDE_DIRECTORIES(Include ../TDI/Include ../Detector/Include ../OptBenPhaMet/Include ../USO/Include ../ArmResponse/Include ../Orbits/Include/ ../GW/Include/ ../Noise/Include/ ../ToolBox/Include/)
add_library (lisacode2main Src/LISACODE-LISACode.cpp Src/LISACODE-ModSig.cpp)
target_link_libraries(lisacode2main lisacode2tdi lisacode2detector lisacode2obpm lisacode2uso lisacode2armresp lisacode2orbits lisacode2gw lisacode2noise lisacode2toolbox)

add_executable(LISACode2 "Exe/LISACODE-Run.cpp")
target_link_libraries(LISACode2 lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2Sensitivty "Exe/LISACODE-Sensitivity.cpp")
target_link_libraries(LC2Sensitivty lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2PSD "Exe/LISACODE-PSD.cpp")
target_link_libraries(LC2PSD lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2SNR "Exe/LISACODE-SNR.cpp")
target_link_libraries(LC2SNR lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2FIM "Exe/LISACODE-FIM.cpp")
target_link_libraries(LC2FIM lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2SNRFIMMulti "Exe/LISACODE-SNRFIMMulti.cpp")
target_link_libraries(LC2SNRFIMMulti lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2MCMC "Exe/LISACODE-MCMC.cpp")
target_link_libraries(LC2MCMC lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2MeanFiles "Exe/LISACODE-MeanFiles.cpp")
target_link_libraries(LC2MeanFiles lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

add_executable(LC2ManipSensi "Exe/LISACODE-ManipSensitivity.cpp")
target_link_libraries(LC2ManipSensi lisacode2main lisacode2armresp lisacode2detector lisacode2gw  lisacode2noise lisacode2obpm lisacode2orbits lisacode2tdi lisacode2toolbox lisacode2uso ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})

install(TARGETS LC2FIM LC2SNR LC2PSD LC2Sensitivty LISACode2 LC2SNRFIMMulti LC2MeanFiles LC2MCMC LC2SNRFIMMulti LC2ManipSensi DESTINATION /usr/local/bin)

install(FILES Exe/makeTDI-lisacode2.py DESTINATION /usr/local/bin)

#ADD_CUSTOM_COMMAND(TARGET Exe/makeTDI-lisacode2.py POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#        "${CMAKE_CURRENT_BINARY_DIR}/Exe/makeTDI-lisacode2.py"
#        "/usr/local/bin")
