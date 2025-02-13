cmake_minimum_required(VERSION 3.14)
cmake_policy(SET CMP0057 NEW)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

set(DEAL_II_DIR "D:/Documents/Programmation/CPP/libs/libdealii/lib/cmake/deal.II")
set(TBB_PATH "C:/Program Files\ (x86)/Intel/oneAPI/tbb/latest")
set(TBB_CMAKE_DIR "C:/Program Files\ (x86)/Intel/oneAPI/tbb/latest/lib/cmake/tbb")
set(MUMPS_DIR "D:/Documents/Programmation/CPP/libs/libmumps/cmake")

find_package(deal.II 9.6.2 REQUIRED HINTS ${DEAL_II_DIR})

if(NOT ${deal.II_FOUND})
    message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()
if(NOT ${TBB_FOUND})
    message(FATAL_ERROR "\n"
    "*** Could not locate a  version of oneTBB. ***\n\n"
    "You may want to either pass a flag -DTBB_DIR=/path/to/oneTBB to cmake\n"
    "or set an environment variable \"TBB_DIR\" that contains this path."
    )
endif()
deal_ii_initialize_cached_variables()

project(eddy_solver LANGUAGES CXX)

find_package(TBB HINTS TBB_CMAKE_DIR)

add_executable(exe_eddy_solveur ./src/main.cpp)
target_sources(exe_eddy_solveur PUBLIC ./src/eddysolver.cpp ./src/physical_parameters.cpp ./src/zparameters.cpp ./src/source_parameters.cpp)
# target_sources(exe_eddy_solveur PUBLIC ./src/eddysolver.cpp ./src/physical_parameters.cpp ./src/mumps_solver.cpp)
target_include_directories(exe_eddy_solveur PUBLIC ./include ../../../libs/libmumps/include)

deal_ii_setup_target(exe_eddy_solveur)

# TARGET_LINK_LIBRARIES(exe_eddy_solveur TBB::tbb ${METIS_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
TARGET_LINK_LIBRARIES(exe_eddy_solveur TBB::tbb)
# TARGET_LINK_LIBRARIES(exe_eddy_solveur TBB::tbb ../../../../libs/libmumps/lib/dmumps ../../../../libs/libmumps/lib/mumps_common ../../../../libs/libmumps/lib/pord ../../../../libs/libmumps/lib/mpiseq)

if(WIN32)
    set_property(TARGET exe_eddy_solveur PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")
    # add_custom_target(link_dll ALL COMMAND ${CMAKE_COMMAND} -E create_hardlink "${TBB_PATH}/bin/tbb12.dll" "${eddy_solver_BINARY_DIR}/tbb12.dll")
    # file(CREATE_LINK "${TBB_PATH}/bin/tbb12.dll" "${eddy_solver_BINARY_DIR}/tbb12.dll")
    add_custom_command(TARGET exe_eddy_solveur POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy "${TBB_PATH}/bin/tbb12.dll" ${eddy_solver_BINARY_DIR})
endif()

# if(WIN32)
#     add_custom_command(TARGET exe_eddy_solveur POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy D:/Documents/Programmation/CPP/libs/libopenblas/bin/libopenblas.dll ${eddy_solver_BINARY_DIR})
# endif()