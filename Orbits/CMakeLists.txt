project(Orbits)
INCLUDE_DIRECTORIES(Include ../ToolBox/Include/)
#file(GLOB_RECURSE Orbits_SOURCES "Src/*.c*")
add_library (lisacode2orbits 	Src/LISACODE-Orbits.cpp Src/LISACODE-OrbitsAnaLISA.cpp Src/LISACODE-OrbitsAnaMLDC.cpp Src/LISACODE-OrbitsData.cpp Src/LISACODE-OrbitsAnaOctahedron.cpp)
target_link_libraries(lisacode2orbits lisacode2toolbox)

add_executable(LC2Orbits "Test/LISACODE-TestOrbits.cpp")
target_link_libraries(LC2Orbits lisacode2orbits lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
