cmake_minimum_required(VERSION 2.8)

project(REPPnP)

IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF ()

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

# Better support for external libs (-fPIC flag)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3 ")

# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  add_definitions(-DCOMPILEDWITHC11)
  message(STATUS "Using flag -std=c++11.")
elseif (COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  add_definitions(-DCOMPILEDWITHC0X)
  message(STATUS "Using flag -std=c++0x.")
else ()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif ()


LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(OpenCV 3.0 QUIET)
if (NOT OpenCV_FOUND)
  find_package(OpenCV 2.4.3 QUIET)
  if (NOT OpenCV_FOUND)
    message(FATAL_ERROR "OpenCV > 2.4.3 not found.")
  endif ()
endif ()


#Find External Packages to refer to them later.
find_package(Eigen3 3.1.0 REQUIRED)
find_package(GTest REQUIRED)

#Add include Header files
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/transform/include
    ${EIGEN3_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}/unsupported
    ${GTEST_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/transform/test/eigen_checks/include
    )
message(${EIGEN3_INCLUDE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

#Build Libraries of this Project
add_library(${PROJECT_NAME} SHARED
    include/reppnp.h src/reppnp.cpp
    transform/include/transform_interface.h
    transform/include/transform.h transform/src/transform.cpp
    transform/include/scaled_transform.h transform/src/scaled_transform.cpp)

#Link External Libs to Project
target_link_libraries(
    ${PROJECT_NAME}
    ${EIGEN3_LIBS}
    glog
    pthread
    )

# Testing
#Build executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

#Build Library for testing
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)


add_executable(TestExecutable_REPPnP test/test_REPPnP.cpp)
target_link_libraries(TestExecutable_REPPnP ${GTEST_LIBRARIES} ${PROJECT_NAME} )