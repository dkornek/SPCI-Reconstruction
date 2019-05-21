cmake_minimum_required(VERSION 3.7)
project(SPCI-Reconstruction)

##### FILE GLOBBING #####
file(GLOB_RECURSE sources src/*.cpp src/*.h)

##### TARGET DEFINITIONS #####
# Add data to the target
add_executable(SPCI-Reconstruction ${sources})

# Include ROOT to be compiled
target_compile_options(SPCI-Reconstruction PUBLIC -pthread -std=c++11
    -m64 -I/opt/root/root-6.12.04/include -L/opt/root/root-6.12.04/lib
    -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript
    -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
    -Wall -Wextra -Wfloat-conversion -fopenmp)

##### DEPENDENCIES #####
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Boost 1.60.0)
target_link_libraries(SPCI-Reconstruction PUBLIC ${Boost_LIBRARIES})

# list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT 6.12.04 EXACT)
target_link_libraries(SPCI-Reconstruction PUBLIC ${ROOT_LIBRARIES})

##### PACKAGING #####
install(TARGETS SPCI-Reconstruction DESTINATION SPCI-Reconstruction_destination)
set(CPACK_PACKAGE_NAME "SPCI-Reconstruction")
set(CPACK_PACKAGE_VERSION "1.0.0")

# Build 1 executable file
set(CPACK_MONOLITHIC_INSTALL 1)

include(CPack)
