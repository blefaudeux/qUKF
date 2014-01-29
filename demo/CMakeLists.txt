# Define the project's name
project(Demo)

cmake_minimum_required(VERSION 2.6)

# Make sure the compiler can find include files from our UKF library. 
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/headers)
include_directories(${PROJECT_SOURCE_DIR}/../libUKF/headers)

set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")

if(NOT EIGEN3_INCLUDE_DIR )
  # Try to manually define Eigen location
  set(ENV{EIGEN3_INCLUDE_DIR} "/usr/include/eigen3/")
  set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")

  # Raise error if that's not enough
  if(NOT EIGEN3_INCLUDE_DIR )
    message(FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
  endif()
endif()

include_directories( "${EIGEN3_INCLUDE_DIR}" )

# Try to find the needed packages
find_package(Boost REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add all the files we're interested in (mostly for QTCreator)
file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/headers/*.h)

MESSAGE("Cpp files : \n ${SRC}\n")
MESSAGE("H files : \n ${HEADERS}\n")

# Add the needed packages
find_package( OpenCV REQUIRED )

# Make sure the linker can find the libUKF library once it is built. 

set(LIBUKF_BINARY_DIR ${PROJECT_SOURCE_DIR}/../libUKF-build)

link_directories (${LIBUKF_BINARY_DIR}) 
MESSAGE("LibUKF Path :\n ${LIBUKF_BINARY_DIR}\n")

# Add executable "demo"
add_executable(demo ${SRC} ${HEADERS})

# Link the executable to the UKF library.
target_link_libraries (demo libUKF_dynamic ${OpenCV_LIBS})
