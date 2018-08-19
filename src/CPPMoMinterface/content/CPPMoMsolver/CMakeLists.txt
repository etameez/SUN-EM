cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

project(mom CXX C Fortran)

include_directories("lib/eigen/")

option(MPI "Compile with MPI" OFF)
option(ATLAS "Use optimized BLAS" OFF)
option(PARALLEL "Explicit Parallel Atlas" OFF)

if(MPI)
    MESSAGE("COMPILING WITH MPI")

    find_package(MPI REQUIRED)
    find_package(OpenMP REQUIRED)
    #find_package(ATLAS REQUIRED)

    if(ATLAS)
        message("COMPILING WITH ATLAS")
        find_package(LAPACK REQUIRED)
        find_package(BLAS REQUIRED)
        SET(SCALAPACK_DIR "/usr/lib/scalapack/build")
        find_package(SCALAPACK REQUIRED)
        SET(SCALAPACK_LIBRARIES "/usr/lib/scalapack/build/lib/libscalapack.so")
    elseif(PARALLEL)
        message("COMPILING PARALLEL ATLAS")
        SET(LAPACK_LIBRARIES "/usr/lib/libptlapack.a")
        SET(BLAS_LIBRARIES "/usr/lib/libtatlas.so")
        SET(SCALAPACK_LIBRARIES "/usr/lib/scalapack/build/lib/libscalapack.so")
    endif()

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
                src/mom_file_writer.cpp
                src/mpi/main_mpi.cpp
                src/mpi/mom_solver_mpi.cpp)
    add_executable(mom_mpi ${SOURCES})

    if(ATLAS OR PARALLEL)
        target_link_libraries(mom_mpi ${MPI_LIBRARIES}
                                    OpenMP::OpenMP_CXX
                                    ${SCALAPACK_LIBRARIES}
                                    ${BLAS_LIBRARIES}
                                    ${LAPACK_LIBRARIES})
    else()
        target_link_libraries(mom_mpi ${MPI_LIBRARIES}
                                    OpenMP::OpenMP_CXX
                                    -L/${CMAKE_CURRENT_SOURCE_DIR}/lib/scalapack
                                    -lscalapack
                                    -lrefblas
                                    -ltmg
                                    -lreflapack
                                    -L/usr/bin
                                    -lgfortran)
    endif()

else()
    MESSAGE("COMPILING WITHOUT MPI")
    file(GLOB SOURCES "src/*.cpp")
    add_executable(mom ${SOURCES})
    target_link_libraries(mom   -L/${CMAKE_CURRENT_SOURCE_DIR}/lib/scalapack
                                -lreflapack
                                -lrefblas
                                -ltmg 
                                -L/usr/bin
                                -lgfortran)
endif(MPI)
unset(MPI CACHE)
