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
# learn_obs
############################################################

add_executable(learn_obs src/learn_obs.cpp src/dbm_obs.hpp)

target_link_libraries(learn_obs PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_obs PUBLIC ${Q3C1_LIB})

############################################################
# learn_params_rbm
############################################################

add_executable(learn_params_rbm src/learn_params_rbm.cpp src/dbm_params_rbm.hpp)

target_link_libraries(learn_params_rbm PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_params_rbm PUBLIC ${Q3C1_LIB})

############################################################
# learn_params_rbm_centered
############################################################

add_executable(learn_params_rbm_centered src/learn_params_rbm_centered.cpp src/dbm_params_rbm_centered.hpp)

target_link_libraries(learn_params_rbm_centered PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_params_rbm_centered PUBLIC ${Q3C1_LIB})

############################################################
# learn_params_dbm_centered
############################################################

add_executable(learn_params_dbm_centered src/learn_params_dbm_centered.cpp src/dbm_params_dbm_centered.hpp)

target_link_libraries(learn_params_dbm_centered PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_params_dbm_centered PUBLIC ${Q3C1_LIB})

############################################################
# learn_params_alternating_binary
############################################################

add_executable(learn_params_alternating_binary src/learn_params_alternating_binary.cpp src/dbm_params_alternating_binary.hpp)

target_link_libraries(learn_params_alternating_binary PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_params_alternating_binary PUBLIC ${Q3C1_LIB})

############################################################
# cd_sampling
############################################################

add_executable(cd_sampling src/cd_sampling.cpp src/dbm_ic.hpp)

target_link_libraries(cd_sampling PUBLIC ${DBLZ_LIB})
target_link_libraries(cd_sampling PUBLIC ${Q3C1_LIB})

############################################################
# sample_traj
############################################################

add_executable(sample_traj src/sample_traj.cpp src/dbm_ic.hpp)

target_link_libraries(sample_traj PUBLIC ${DBLZ_LIB})
target_link_libraries(sample_traj PUBLIC ${Q3C1_LIB})

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)
