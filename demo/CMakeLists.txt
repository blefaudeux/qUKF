# Define the project's name
project(Demo)

cmake_minimum_required(VERSION 2.6)

# Make sure the compiler can find include files from our UKF library. 
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/headers)
include_directories(${PROJECT_SOURCE_DIR}/../libUKF/headers)
include_directories(${PROJECT_SOURCE_DIR}/../libUKF/lib)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# GCC flags
SET( GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -Wall -Wextra -pedantic")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

# Try to find the needed packages
find_package(Boost REQUIRED)

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

find_package( OpenCV REQUIRED )

# Set debug/release flags
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

# Add all the files we're interested in (mostly for QTCreator)
file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/headers/*.h)

# Make sure the linker can find the libUKF library once it is built. 
set(LIBUKF_BINARY_DIR ${PROJECT_SOURCE_DIR}/../libUKF/lib)

link_directories (${LIBUKF_BINARY_DIR}) 
MESSAGE("LibUKF Path :\n ${LIBUKF_BINARY_DIR}\n")

# Add executable "demo"
add_executable(demo ${SRC} ${HEADERS})

# Link the executable to the UKF library.
target_link_libraries (demo UKFd ${OpenCV_LIBS})
