cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(learn)

include_directories(/usr/local/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

find_library(DBLZ_LIB dblz HINTS /usr/local/lib/)
find_library(Q3C1_LIB q3c1 HINTS /usr/local/lib/)

############################################################
# diagnose_sampling
############################################################

add_executable(diagnose_sampling src/diagnose_sampling.cpp src/dbm_static.hpp)

target_link_libraries(diagnose_sampling PUBLIC ${DBLZ_LIB})
target_link_libraries(diagnose_sampling PUBLIC ${Q3C1_LIB})

############################################################
# learn_centered
############################################################

add_executable(learn_centered src/learn_centered.cpp src/dbm_centered.hpp)

target_link_libraries(learn_centered PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_centered PUBLIC ${Q3C1_LIB})

############################################################
# sample_traj
############################################################

add_executable(sample_traj src/sample_traj.cpp src/dbm_static.hpp)

target_link_libraries(sample_traj PUBLIC ${DBLZ_LIB})
target_link_libraries(sample_traj PUBLIC ${Q3C1_LIB})

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)
