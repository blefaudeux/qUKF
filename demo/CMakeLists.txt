# Define the project's name
project(UKF_DEMO)

cmake_minimum_required(VERSION 2.6)

# Make sure the compiler can find include files from our UKF library. 
include_directories(${UKF_DEMO_SOURCE_DIR}/src)
include_directories(${UKF_DEMO_SOURCE_DIR}/headers)
include_directories(${UKF_DEMO_SOURCE_DIR}/../libUKF)

# Add all the files we're interested in (mostly for QTCreator)
file(GLOB ${UKF_DEMO_SOURCE_DIR}/src *.cpp)
file(GLOB ${UKF_DEMO_SOURCE_DIR}/headers *.h)

# Add the needed packages
find_package( OpenCV REQUIRED )

# Make sure the linker can find the libUKF library once it is built. 
link_directories (${LIBUKF_BINARY_DIR}) 

# Add executable "demo"
add_executable (demo src/test_unit_qUKF.cpp) 

# Link the executable to the Hello library. 
target_link_libraries (demo libUKF)  
