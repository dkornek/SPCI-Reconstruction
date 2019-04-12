cmake_minimum_required(VERSION 3.7)
project(SPCI-Reconstruction)

##### FILE GLOBBING #####
file(GLOB_RECURSE sources src/*.cpp src/*.h)

##### TARGET DEFINITIONS #####
# Add data to the target
add_executable(SPCI-Reconstruction ${sources})

# Include ROOT to be compiled
# execute_process(COMMAND root-config --cflags --glibs OUTPUT_VARIABLE root_flags)
# target_compile_options(SPCI-Reconstruction PUBLIC ${root_flags} -Wall -Wextra -Wfloat-conversion)  not working :(
target_compile_options(SPCI-Reconstruction PUBLIC -pthread -std=c++11
    -m64 -I/opt/root/root-6.12.04/include -L/opt/root/root-6.12.04/lib
    -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript
    -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
    -Wall -Wextra -Wfloat-conversion)

##### DEPENDENCIES #####
find_package(Boost 1.60.0)

target_link_libraries(SPCI-Reconstruction PUBLIC ${Boost_LIBRARIES})

##### PACKAGING #####
install(TARGETS SPCI-Reconstruction DESTINATION SPCI-Reconstruction_destination)
set(CPACK_PACKAGE_NAME "SPCI-Reconstruction")
set(CPACK_PACKAGE_VERSION "1.0.0")

# Build 1 executable file
set(CPACK_MONOLITHIC_INSTALL 1)

include(CPack)
