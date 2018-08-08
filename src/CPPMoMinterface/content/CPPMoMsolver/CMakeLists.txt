cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

project(mom)

include_directories("lib/eigen/")

option(MPI "Compile with MPI" OFF)
if(MPI)
    MESSAGE("COMPILING WITH MPI")
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    set(SOURCES src/edge.cpp
                src/mom_data.cpp
                src/mom_file_reader.cpp
                src/mom_solver.cpp
                src/node.cpp
                src/quadrature.cpp
                src/timer.cpp
                src/triangle.cpp
                src/vrhs.cpp
                src/mpi/main_mpi.cpp
                src/mpi/mom_solver_mpi.cpp)
    add_executable(mom_mpi ${SOURCES})
    target_link_libraries(mom_mpi ${MPI_LIBRARIES})
else()
    MESSAGE("COMPILING WITHOUT MPI")
    file(GLOB SOURCES "src/*.cpp")
    add_executable(mom ${SOURCES})
endif(MPI)
unset(MPI CACHE)
