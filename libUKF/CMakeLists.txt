# Define the project's name
project(LIBUKF)

cmake_minimum_required(VERSION 2.6)

# Add all the files we're interested in (mostly for QTCreator)
file(GLOB ${LIBUKF_SOURCE_DIR}/src *.cpp)
file(GLOB ${LIBUKF_SOURCE_DIR}/headers *.h)

# Define the included directories
include_directories(${LIBUKF_SOURCE_DIR}/src)
include_directories(${LIBUKF_SOURCE_DIR}/headers)

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

# Create static and dynamic libraries
add_library (libUKF_static STATIC src/unscented_KF.cpp) 
add_library (libUKF_dynamic SHARED src/unscented_KF.cpp) 

# Set the link libraries :
TARGET_LINK_LIBRARIES(libUKF_static ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(libUKF_dynamic ${Boost_LIBRARIES} )

# Define the installation procedure
install(TARGETS libUKF_static DESTINATION lib)
install(FILES unscented_KF.h DESTINATION include)
