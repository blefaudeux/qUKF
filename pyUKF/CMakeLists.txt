# Define the project's name
project(pyUKF)
cmake_minimum_required(VERSION 2.6)

# Get to know where all the dependencies are..
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../libUKF/headers)
message("INFO: Additional lookup directory:\n ${CMAKE_MODULE_PATH}")

# Add all the files we're interested in (mostly for QTCreator)
file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Use Boost.python
set(Boost_USE_STATIC_LIBS  OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED)

# Find Python 2.7 (Boost <1.54)
find_path(PYTHON_INCLUDE_DIR python2.7/Python.h HINTS "/usr/include")
set(PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR}/python2.7)
find_library(PYTHON_LIBRARIES NAMES python2.7 libpython2.7 HINTS "/usr/lib/x86_64-linux-gnu")

message("\nINFO: Python libs : \n ${PYTHON_INCLUDE_DIR}\n ${PYTHON_LIBRARIES}\n")

find_package(eigen3)
if(NOT eigen3_FOUND)
  set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")

  if(NOT EIGEN3_INCLUDE_DIR )
    # Try to manually define Eigen location
    set(ENV{EIGEN3_INCLUDE_DIR} "/usr/include/eigen3/")
    set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")

    # Raise error if that's not enough
    if(NOT EIGEN3_INCLUDE_DIR )
      message(FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
    endif()

    include_directories( "${EIGEN3_INCLUDE_DIR}" )
  endif()
endif()

# Get the UKF lib and includes:

## Find system-wide
#find_path(UKF_INCLUDE_DIR qUKF/unscented_KF.h
#          HINTS "/usr/local/include"
#          PATH_SUFFIXES qUKF)

#find_library(UKF_LIBRARY NAMES UKFd libUKFd
#             HINTS "/usr/local/lib")

# Only use the local versions
find_path(UKF_INCLUDE_DIR unscented_KF.h
          HINTS "${PROJECT_SOURCE_DIR}/../libUKF"
          PATH_SUFFIXES include)

find_library(UKF_LIBRARY NAMES UKFd libUKFd
             HINTS ${UKF_INCLUDE_DIR}/../lib)

# Define the included directories
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${UKF_INCLUDE_DIR})

message("\nINFO: UKF lib found in :\n ${UKF_LIBRARY}\n ${UKF_INCLUDE_DIR}\n")

find_package(Boost COMPONENTS python thread)
message("\nINFO: Boost lib found in :\n ${Boost_LIBRARIES}\n")

# Set debug or release configurations
if( CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]" )
    message( "INFO: DEBUG BUILD" )
    # for debug type builds, turn on verbose makefiles
    SET(CMAKE_VERBOSE_MAKEFILE ON)

    # Tell other CMake files that we're doing a debug build
    SET( DEBUG_BUILD 1 )

    # Tell C/C++ that we're doing a debug build
    ADD_DEFINITIONS( -DDEBUG )
endif()

if( CMAKE_BUILD_TYPE MATCHES "[Rr][Ee][Ll][Ee][Aa][Ss][Ee]" )
    message( "INFO: RELEASE BUILD" )
endif()

# GCC flags
SET( GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -Wall -Wextra -pedantic")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

# Create static and dynamic libraries
add_library(pyUKF SHARED ${SRC})

# Fix the static compilation issue
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(pyUKF PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

# Set the link libraries :
set(LINK_LIBS  ${UKF_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
target_link_libraries(pyUKF ${LINK_LIBS} )

### Define the installation procedure
install(TARGETS pyUKF DESTINATION ${PROJECT_SOURCE_DIR}/Python)
