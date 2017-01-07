project(ode-ext)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_definitions("-W -Wall -std=c++14")

set(CMAKE_BUILD_TYPE Release)

set(BOOST_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../boost)
set(KV_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../kv)
set(ODE_EXT_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../ode-ext)

include_directories(${BOOST_INCLUDE} ${KV_INCLUDE} ${ODE_EXT_INCLUDE})

add_executable(${PROJECT_NAME} ${SRC_LIST} main.cpp ${ODE_EXT_INCLUDE}/psa-ext.hpp ${ODE_EXT_INCLUDE}/ode-maffine-ext.hpp)
