CONFIGURE_FILE( Base.hh.in ${CMAKE_SOURCE_DIR}/tests/Base.hh )

cmake_minimum_required(VERSION 3.2.2)

#For testing:
include_directories(${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/include)
ADD_EXECUTABLE( test_jspec test_jspec.cc )
ADD_EXECUTABLE( test_eic test_eic.cc )
ADD_EXECUTABLE( test_functions test_functions.cc )
ADD_EXECUTABLE( test_beam test_beam.cc )
ADD_EXECUTABLE( test_ecooling test_ecooling.cc )

TARGET_LINK_LIBRARIES( test_jspec libjspec gsl gslcblas muparser )
TARGET_LINK_LIBRARIES( test_eic libjspec gsl gslcblas muparser )
TARGET_LINK_LIBRARIES( test_functions libjspec muparser )
TARGET_LINK_LIBRARIES( test_beam libjspec )
TARGET_LINK_LIBRARIES( test_ecooling libjspec gsl gslcblas )

ADD_TEST( ecooling test_ecooling )
ADD_TEST( functions test_functions )
ADD_TEST( beam test_beam )
ADD_TEST( forces_eic test_eic )
ADD_TEST( jspec test_jspec )

##Using as an example https://github.com/Pseudomanifold/Aleph/blob/master/tests/CMakeLists.txt
## as described in http://bastian.rieck.me/blog/posts/2017/simple_unit_tests/
