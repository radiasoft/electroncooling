cmake_minimum_required(VERSION 3.2.2)
project (jspec)
set (CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 11)

#set (CMAKE_CXX_FLAGS "-Wall")
set (CMAKE_CXX_FLAGS "-fopenmp")

#For the shared library:
set (PROJECT_LINK_LIBS muparser gslcblas gsl)
link_directories(lib)

include_directories(include)

file(GLOB SOURCES "src/*.cc")

#Build the library, which we will need for testing
# We don't want the main.cc to be a part of the test suite
list(FILTER SOURCES EXCLUDE REGEX ".*main.cc$")
add_library(lib${PROJECT_NAME} ${SOURCES})

#Build the executable from main.cc
add_executable(${PROJECT_NAME} ${SOURCES} "src/main.cc")

target_link_libraries(${PROJECT_NAME} ${PROJECT_LINK_LIBS})


#For testing:
ENABLE_TESTING()
ADD_SUBDIRECTORY( tests )
