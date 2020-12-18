set(CG_OBJS cg.cpp)

add_library(cg SHARED ${CG_OBJS})
target_include_directories(cg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_dependencies(cg Eigen)
include_directories(${EIGEN_INCLUDE_DIRS})

if(ENABLE_OPENMP)
  target_link_libraries(cg PUBLIC OpenMP::OpenMP_CXX)
endif()
