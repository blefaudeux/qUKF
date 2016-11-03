cmake_minimum_required (VERSION 2.6) 
project (qUKF) 

SET(DEMO            FALSE     CACHE  BOOL    "Build the demo app")

add_subdirectory (libUKF) 


if(DEMO)
	add_subdirectory (demo) 
	add_dependencies( demo UKFd)
endif()

