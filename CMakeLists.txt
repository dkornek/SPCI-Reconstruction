# Include ROOT to be compiled
# target_compile_options(SPCI-Reconstruction PUBLIC -pthread -std=c++11
#     -m64 -I/opt/root/root-6.12.04/include -L/opt/root/root-6.12.04/lib
#     -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript
#     -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
#     -Wall -Wextra -Wfloat-conversion -fopenmp)

##### DEPENDENCIES #####
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

# find_package(Boost 1.60.0)
# target_link_libraries(SPCI-Reconstruction PUBLIC ${Boost_LIBRARIES})
#
# # list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# find_package(ROOT 6.12.04 EXACT)
# target_link_libraries(SPCI-Reconstruction PUBLIC ${ROOT_LIBRARIES})

#
# #----------------------------------------------------------------------------
# # Locate sources and headers for this project
# include_directories(${PROJECT_SOURCE_DIR}/src)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#
# ##### DEPENDENCIES #####
#
# #---Create  a main program using the library
# add_executable(SPCI-Reconstruction ${sources} ${headers})
# target_link_libraries(SPCI-Reconstruction ${ROOT_LIBRARIES} ${Boost_LIBRARIES})


cmake_minimum_required(VERSION 3.7)
project(SPCI-Reconstruction)

# ##### FILE GLOBBING #####
file(GLOB_RECURSE sources src/*.cpp src/*.h)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)
# find_package(ROOT 6.18.00 EXACT)
find_package(Boost 1.60.0)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# include_directories(${PROJECT_SOURCE_DIR}/src)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#---Create  a main program using the library
# add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# ##### PACKAGING #####
install(TARGETS SPCI-Reconstruction DESTINATION SPCI-Reconstruction_destination)
set(CPACK_PACKAGE_NAME "SPCI-Reconstruction")
set(CPACK_PACKAGE_VERSION "1.0.0")

# Build 1 executable file
set(CPACK_MONOLITHIC_INSTALL 1)

include(CPack)
