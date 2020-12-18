add_executable(test main.cpp)
target_link_libraries(test cg)
include_directories(${EIGEN_INCLUDE_DIRS})
