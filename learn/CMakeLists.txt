cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(learn)

include_directories(/usr/local/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

# Dependencies
include(cmake/CPM.cmake)
CPMFindPackage(
    GITHUB_REPOSITORY "smrfeld/Q3-C1-Finite-Elements"
    NAME q3c1
    GIT_TAG 3.0
    )
CPMFindPackage(
    GITHUB_REPOSITORY "smrfeld/DynamicBoltzmann"
    NAME dblz
    GIT_TAG 4.5
    )

############################################################
# diagnose_sampling
############################################################

add_executable(diagnose_sampling src/diagnose_sampling.cpp src/dbm_static.hpp)

# Include dirs for q3c1, dblz
target_include_directories(diagnose_sampling PRIVATE ${q3c1_SOURCE_DIR}/include/)
target_include_directories(diagnose_sampling PRIVATE ${dblz_SOURCE_DIR}/include/)

target_link_libraries(diagnose_sampling PUBLIC dblz)
target_link_libraries(diagnose_sampling PUBLIC q3c1)

############################################################
# learn_centered
############################################################

add_executable(learn_centered src/learn_centered.cpp src/dbm_centered.hpp)

# Include dirs for q3c1, dblz
target_include_directories(learn_centered PRIVATE ${q3c1_SOURCE_DIR}/include/)
target_include_directories(learn_centered PRIVATE ${dblz_SOURCE_DIR}/include/)

target_link_libraries(learn_centered PUBLIC dblz)
target_link_libraries(learn_centered PUBLIC q3c1)

############################################################
# sample_traj
############################################################

add_executable(sample_traj src/sample_traj.cpp src/dbm_static.hpp)

# Include dirs for q3c1, dblz
target_include_directories(sample_traj PRIVATE ${q3c1_SOURCE_DIR}/include/)
target_include_directories(sample_traj PRIVATE ${dblz_SOURCE_DIR}/include/)

target_link_libraries(sample_traj PUBLIC dblz)
target_link_libraries(sample_traj PUBLIC q3c1)

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)
