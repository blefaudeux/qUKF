# Make sure the compiler can find include files from our Hello library. 
include_directories (${QUKF_SOURCE_DIR}/libUKF) 

# Make sure the linker can find the libUKF library once it is built. 
link_directories (${QUKF_BINARY_DIR}/libUKF) 

# Add executable called "helloDemo" that is built from the source files 
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found. 
add_executable (demo test_unit_qUKF.cxx) 

# Link the executable to the Hello library. 
target_link_libraries (demo libUKF)  
