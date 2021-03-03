
#find_package(Eigen3 REQUIRED)
find_package(OpenMP)
find_package(Boost COMPONENTS program_options REQUIRED)

set (SOURCES
    single_rate_sine.cpp
    multi_rate_sine.cpp
    input.cpp
)

# TODO: Remove hard-coded paths
set(Boost_INCLUDE_DIRS "C:\\Lib\\Boost\\include\\boost-1_70\\")
#include_directories(${Boost_INCLUDE_DIRS})

#message($ENV{YAML_CPP_INCLUDE_DIR})
message("Yaml path: $ENV{YAML_CPP_LIB_DIR}")

# 
include_directories(${YAML_CPP_INCLUDE_DIRS})

add_executable(friction_dynamics main.cpp ${SOURCES})
target_include_directories(friction_dynamics PUBLIC $ENV{EIGEN3_INCLUDE_DIR} $ENV{YAML_CPP_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
target_link_libraries(friction_dynamics PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES})
target_link_libraries(friction_dynamics PRIVATE yaml-cpp)
#target_link_directories(friction_dynamics PUBLIC ${Boost_INCLUDE_DIR})
target_link_directories(friction_dynamics PUBLIC "C:\\Lib\\Boost\\lib")