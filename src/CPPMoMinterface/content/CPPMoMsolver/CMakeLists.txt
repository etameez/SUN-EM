cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)
project(mom)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

file(GLOB SOURCES "src/*.cpp")

add_executable(mom ${SOURCES})
target_link_libraries(mom ${MPI_LIBRARIES})
