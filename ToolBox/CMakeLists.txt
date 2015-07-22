project(ToolBox)
INCLUDE_DIRECTORIES(Include)

#file(GLOB_RECURSE ToolBox_SOURCES "Src/*.c*")
add_library (lisacode2toolbox 	Src/ezxml.c Src/randlib/randlib.c Src/randlib/com.c Src/randlib/linpack.c Src/LISACODE-Tools.cpp Src/LISACODE-Serie2.cpp Src/LISACODE-Vector.cpp Src/LISACODE-DataFileWrite.cpp Src/LISACODE-DataFileRead.cpp Src/LISACODE-EllipticFilter.cpp Src/LISACODE-Filter.cpp Src/LISACODE-Matrix.cpp)

add_executable(LC2TestTools "Test/LISACODE-TestTools.cpp")
target_link_libraries(LC2TestTools lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
add_executable(LC2xml2ascii "Test/LISACODE-xml2ascii.cpp")
target_link_libraries(LC2xml2ascii lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
add_executable(LC2zDL "Test/LISACODE-zDL.cpp")
target_link_libraries(LC2zDL lisacode2toolbox ${FFTW3_LIBRARIES} ${GSL_LIBRARIES})
