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
# diagnose_dbm
############################################################

add_executable(diagnose_dbm src/diagnose_dbm.cpp src/dbm_static.hpp)

target_link_libraries(diagnose_dbm PUBLIC ${DBLZ_LIB})
target_link_libraries(diagnose_dbm PUBLIC ${Q3C1_LIB})

############################################################
# diagnose_rbm
############################################################

add_executable(diagnose_rbm src/diagnose_rbm.cpp src/rbm_static.hpp)

target_link_libraries(diagnose_rbm PUBLIC ${DBLZ_LIB})
target_link_libraries(diagnose_rbm PUBLIC ${Q3C1_LIB})

############################################################
# learn_alternating_binary
############################################################

add_executable(learn_alternating_binary src/learn_alternating_binary.cpp src/dbm_alternating_binary.hpp)

target_link_libraries(learn_alternating_binary PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_alternating_binary PUBLIC ${Q3C1_LIB})

############################################################
# learn_dbm_centered
############################################################

add_executable(learn_dbm_centered src/learn_dbm_centered.cpp src/dbm_centered.hpp)

target_link_libraries(learn_dbm_centered PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_dbm_centered PUBLIC ${Q3C1_LIB})

############################################################
# learn_rbm_centered
############################################################

add_executable(learn_rbm_centered src/learn_rbm_centered.cpp src/rbm_centered.hpp)

target_link_libraries(learn_rbm_centered PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_rbm_centered PUBLIC ${Q3C1_LIB})

############################################################
# learn_rbm_obs
############################################################

add_executable(learn_rbm_obs src/learn_rbm_obs.cpp src/rbm_obs.hpp)

target_link_libraries(learn_rbm_obs PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_rbm_obs PUBLIC ${Q3C1_LIB})

############################################################
# learn_rbm
############################################################

add_executable(learn_rbm src/learn_rbm.cpp src/rbm.hpp)

target_link_libraries(learn_rbm PUBLIC ${DBLZ_LIB})
target_link_libraries(learn_rbm PUBLIC ${Q3C1_LIB})

############################################################
# rbm_sample_traj
############################################################

add_executable(rbm_sample_traj src/rbm_sample_traj.cpp src/rbm_centered.hpp)

target_link_libraries(rbm_sample_traj PUBLIC ${DBLZ_LIB})
target_link_libraries(rbm_sample_traj PUBLIC ${Q3C1_LIB})

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)
