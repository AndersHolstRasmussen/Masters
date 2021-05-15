cmake_minimum_required(VERSION 2.8)
project(I257)

# Tell cmake where to look for cmake modules.
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")

#Lets write modern C++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11")


#----------------------------------------------------------------------------
# Dependencies
#

# Variable for holding dependant libraries
set(LIBRARIES "")
SET(CMAKE_C_LINK_EXECUTABLE ${CMAKE_CXX_LINK_EXECUTABLE})
# set(LIBRARIES "-lboost_filesystem -lboost_system")


# Include AUSAlib. Look in cmake/FindAUSALIB.cmake for how the specify location
find_package(AUSALIB REQUIRED)
include_directories(${AUSALIB_INCLUDES})
list(APPEND LIBRARIES ${AUSALIB_LIBRARIES})

# Include ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
list(APPEND LIBRARIES ${ROOT_LIBRARIES} Spectrum MathMore)

find_package(CURL)
if (CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIRS})
    list(APPEND LIBRARIES ${CURL_LIBRARIES})
    add_definitions("-DAUSA_HAS_LIBCURL")
endif(CURL_FOUND)



link_libraries(${LIBRARIES})
include_directories(external/jsoncons/include)

add_executable(Li8 Li8Analysis.cpp)
add_executable(plt plots.cpp)
add_executable(my_eff Efficiency.cpp)
add_executable(eff eff.cpp)
