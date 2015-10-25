cmake_minimum_required (VERSION 2.6) 
project (qUKF) 

# Recurse into the "qUKF" and "demo" subdirectories. This does not actually 
# cause another cmake executable to run. The same process will walk through 
# the project's entire directory structure.
add_subdirectory (libUKF) 
add_subdirectory (pyUKF)

add_dependencies( pyUKF UKFd)

if(ADD_DEMO)
	add_subdirectory (demo) 
	add_dependencies( demo UKFd)
endif()


