cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(stoch_sims)

include_directories(/usr/local/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

# Dependencies
include(cmake/CPM.cmake)
CPMAddPackage(
    GITHUB_REPOSITORY "smrfeld/LatticeGillespieCpp"
    NAME lattgillespie
    GIT_TAG 2.0
    )
CPMAddPackage(
    GITHUB_REPOSITORY "smrfeld/DynamicBoltzmann"
    NAME bmla
    GIT_TAG 4.5
    )

############################################################
# main
############################################################

add_executable(main src/main.cpp)

# Include dirs for bmla, lattgillespie
target_include_directories(main PRIVATE ${lattgillespie_SOURCE_DIR}/include/)

target_link_libraries(main PUBLIC lattgillespie)

############################################################
# learn_ic
############################################################

add_executable(calculate_moments src/calculate_moments.cpp)

# Include dirs for bmla
target_include_directories(calculate_moments PRIVATE ${bmla_SOURCE_DIR}/include/)

target_link_libraries(calculate_moments PUBLIC bmla)

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)