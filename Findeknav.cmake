find_path(eknav_INCLUDE_DIR eknav/ins_qkf.hpp ${CMAKE_INSTALL_PREFIX}/include "$ENV{NAMER_ROOT}")
find_library(eknav_LIBRARY eknav ${CMAKE_INSTALL_PREFIX}/lib/eknav "$ENV{NAMER_ROOT}")

find_package(Boost REQUIRED COMPONENTS system thread)
find_package(Eigen3 REQUIRED)
set(eknav_LIBRARIES
  ${eknav_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_THREAD_LIBRARY}
)
set(eknav_INCLUDE_DIRS
  ${eknav_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(eknav DEFAULT_MSG eknav_LIBRARY eknav_INCLUDE_DIR)

# message(STATUS "eknav_INCLUDE_DIR: ${eknav_INCLUDE_DIR}")
# message(STATUS "eknav_LIBRARY: ${eknav_LIBRARY}")
# message(STATUS "Found eknav: ${eknav_FOUND}")
