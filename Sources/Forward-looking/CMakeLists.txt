

cmake_minimum_required(VERSION 2.8.4)
project(ForwardLooking)

include(FindCPLEX.cmake)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/)
set(CMAKE_BINARY_DIR  ${PROJECT_SOURCE_DIR}/)
set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/)
set(PROJECT_BINARY_DIR  ${PROJECT_SOURCE_DIR}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


FIND_PACKAGE(CPLEX REQUIRED)
INCLUDE_DIRECTORIES(${CPLEX_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -fexceptions -fPIC -m64 -Wshadow ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 ${CMAKE_CXX_FLAGS_RELEASE}")

add_executable(ForwardLooking main.cpp Util.h Item.h Objeto.h InstanciaPorPeriodo.h Delta.h scpp_assert.cpp scpp_assert.h scpp_vector.h)

# a project that use cplex must define the macro IL_STD
set_property(
  TARGET ForwardLooking
  PROPERTY COMPILE_DEFINITIONS IL_STD
  PROPERTY COMPILE_DEFINITIONS NDEBUG
  )

TARGET_LINK_LIBRARIES( ForwardLooking pthread dl ${CPLEX_LIBRARIES} ${DRILLHOLES_LIBRARIES})
