CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(SSDR-2012 CXX C)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_BUILD_TYPE Debug)

SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

INCLUDE_DIRECTORIES(external/)

FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})

find_package(Cairo)
include_directories(${CAIRO_INCLUDE_DIRS})

find_package(osqp REQUIRED)

find_package(OsqpEigen REQUIRED)

## GCO
set(gco_srcs
    external/gco/GCoptimization.cpp
    external/gco/GCoptimization.h
    external/gco/LinkedBlockList.cpp
    external/gco/LinkedBlockList.h
    external/gco/block.h
    external/gco/energy.h
    external/gco/graph.h
    )
# INCLUDE_DIRECTORIES(external/gco)

add_library(gco ${gco_srcs})

## GRAPHICS_GEM
set(gem_srcs
    external/Graphics_Gems/GraphicsGems.cpp
    external/Graphics_Gems/GraphicsGems.h
    external/Graphics_Gems/NearestPoint.cpp
    external/Graphics_Gems/NearestPoint.h
    )

# INCLUDE_DIRECTORIES(external/Graphics_Gems)
add_library(gems ${gem_srcs})

MESSAGE("
*****************************************************************************
* REQUIRED EXTERNAL LIBRARIES
*   Eigen   
*   Boost  
*   SVG++ : header only/ headers included in the project
*   Cairo
*   OSQP
*   OSQP-Eigen
*****************************************************************************
")

#############################################################################
INCLUDE_DIRECTORIES(src)

AUX_SOURCE_DIRECTORY(src SRC_FILES)

ADD_EXECUTABLE(svg_parser ${SRC_FILES})
target_link_libraries(svg_parser Eigen3::Eigen)
target_link_libraries(svg_parser gco)
target_link_libraries(svg_parser gems)
target_link_libraries(svg_parser osqp::osqp)
target_link_libraries(svg_parser ${CAIRO_LIBRARIES})
target_link_libraries(svg_parser OsqpEigen::OsqpEigen)
