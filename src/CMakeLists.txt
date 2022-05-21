cmake_minimum_required(VERSION 3.6)
project(new_from_scrath)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp config.cpp config.h distribution.cpp distribution.h export2D.cpp export2D.h position.cpp position.h sphere.cpp sphere.h states.cpp states.h system.cpp system.h timer.cpp timer.h)
add_executable(new_from_scrath ${SOURCE_FILES})