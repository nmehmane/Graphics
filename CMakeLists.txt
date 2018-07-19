CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(SSDR-2012 CXX C)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

INCLUDE_DIRECTORIES(external/)

FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})

MESSAGE("
*****************************************************************************
* REQUIRED EXTERNAL LIBRARIES
*   Eigen   
*   Boost  
*   SVG++ : header only/ headers included in the project
*****************************************************************************
")

#############################################################################
INCLUDE_DIRECTORIES(src)

AUX_SOURCE_DIRECTORY(src SRC_FILES)

ADD_EXECUTABLE(svg_parser ${SRC_FILES})
TARGET_LINK_LIBRARIES(svg_parser Eigen3::Eigen)
