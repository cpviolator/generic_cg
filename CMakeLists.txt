cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the project name and version
project(cg VERSION 1.0.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(ExternalProject)

# Link OpempMP CMake automatically adds 
# _OPENMP to the CXX flags
option(ENABLE_OPENMP "enable OpenMP" ON)
if(ENABLE_OPENMP)
  find_package(OpenMP)
endif()

# Add Eigen for small, dense algebra
set(EIGEN_VERSION 3.3.7)
set(EIGEN_DOWNLOAD_LOCATION ${CMAKE_SOURCE_DIR}/externals/eigen/${EIGEN_VERSION}.tar.bz2)
set(EIGEN_URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.bz2)
set(EIGEN_SHA 685adf14bd8e9c015b78097c1dc22f2f01343756f196acdc76a678e1ae352e11)
if(NOT EXISTS ${EIGEN_DOWNLOAD_LOCATION})
  message(STATUS "Checking for Eigen tarball and downloading if necessary.")
endif()
file(
  DOWNLOAD ${EIGEN_URL} ${EIGEN_DOWNLOAD_LOCATION}
  EXPECTED_HASH SHA256=${EIGEN_SHA}
  STATUS EIGEN_DOWNLOADED
  )
list(GET EIGEN_DOWNLOADED 0 EIGEN_DOWNLOADED_CODE)
list(GET EIGEN_DOWNLOADED 1 EIGEN_DOWNLOADED_MSG)
if(${EIGEN_DOWNLOADED_CODE})
  message( SEND_ERROR "Could not download Eigen automatically (${EIGEN_DOWNLOADED_MSG}). "
    "Please download eigen from ${EIGEN_URL} and save it to ${EIGEN_DOWNLOAD_LOCATION} "
    "and try running cmake again.")
endif()

ExternalProject_Add(
  Eigen
  URL ${CMAKE_SOURCE_DIR}/externals/eigen/${EIGEN_VERSION}.tar.bz2
  URL_HASH SHA256=${EIGEN_SHA}
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/externals/eigen/
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_Property(Eigen source_dir)
set(EIGEN_INCLUDE_DIRS ${source_dir})

add_subdirectory(src)
add_subdirectory(test)

add_subdirectory(pybind11)
pybind11_add_module(cg_python_binding SHARED src/cg_python_binding.cpp src/cg.cpp)
target_include_directories(cg_python_binding PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
