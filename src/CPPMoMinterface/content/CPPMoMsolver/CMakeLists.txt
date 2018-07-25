cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

project(mom)

file(GLOB SOURCES "src/*.cpp")

add_executable(mom ${SOURCES})
