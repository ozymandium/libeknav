project(jbrandmeyer_libeknav)
cmake_minimum_required(VERSION 2.8)

set(LIBEKNAV_OPTIONS "-DRANK_ONE_UPDATES")
set(MORE_OPTIONS "-DEIGEN_DONT_ALIGN -DEIGEN_DONT_VECTORIZE ${LIBEKNAV_OPTIONS}")
set(OPTIMIZE "-Os -ffast-math -g -ffunction-sections -fdata-sections")
set(WARNINGS "-Wall -Wextra -std=c++11")
set(eknav_CXX_FLAGS "-pipe -fshow-column -fPIC ${OPTIMIZE} ${WARNINGS} ${MORE_OPTIONS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${eknav_CXX_FLAGS}")

include_directories(${PROJECT_SOURCE_DIR} posix)

find_package(Boost REQUIRED COMPONENTS system thread)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
list(APPEND eknav_LINK_LIBS
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_THREAD_LIBRARY}
)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# # print all available variables
# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

# set(INS_QKF_SRCS
#   ins_qkf_predict_ned.cpp
#   ins_qkf_observe_gps_p.cpp
#   ins_qkf_observe_vector.cpp
#   ins_qkf_predict.cpp
#   diagnostics.cpp
#   basic_ins_qkf.cpp
#   pr_ins_qkf.hpp
#   pr_ins_qkf.cpp
# )

set(INS_QKF_NED_SRCS
  ins_qkf_observe_vector.cpp
  ins_qkf_predict_ned.cpp
  basic_ins_qkf.cpp
  diagnostics.cpp
  posix/timer.hpp
  posix/timer.cpp
)

set(INS_QKF_ECEF_SRCS
  ins_qkf_observe_gps_p.cpp
  ins_qkf_observe_vector.cpp
  ins_qkf_predict.cpp
  basic_ins_qkf.cpp
  diagnostics.cpp
  posix/timer.hpp
  posix/timer.cpp
  pr_ins_qkf.hpp
  pr_ins_qkf.cpp
)

add_library(eknav-ned ${INS_QKF_NED_SRCS})
target_link_libraries(eknav-ned ${eknav_LINK_LIBS})

add_library(eknav-ecef ${INS_QKF_ECEF_SRCS})
target_link_libraries(eknav-ecef ${eknav_LINK_LIBS})

add_executable(test_ins_qkf
  test_ins_qkf.cpp
  ins_qkf_observe_gps_pvt
  posix/random_seed.hpp
  posix/random_seed.cpp
)
target_link_libraries(test_ins_qkf
  eknav-ecef
  ${eknav_LINK_LIBS}
)

add_executable(monte_carlo_pr_ins_qkf
  monte_carlo_pr_ins_qkf.cpp
  posix/random_seed.hpp
  posix/random_seed.cpp
)
target_link_libraries(monte_carlo_pr_ins_qkf
  eknav-ecef
  ${eknav_LINK_LIBS}
)
