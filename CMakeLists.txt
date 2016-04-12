cmake_minimum_required (VERSION 2.6) 
project (qUKF) 

SET(DEMO            FALSE     CACHE  BOOL    "Build the demo app")
SET(PYTHON_BINDINGS FALSE     CACHE  BOOL    "ADd python bindings -broken as of now-")


add_subdirectory (libUKF) 


if(DEMO)
	add_subdirectory (demo) 
	add_dependencies( demo UKFd)
endif()

if(PYTHON_BINDINGS)
	add_subdirectory (pyUKF)
	add_dependencies( pyUKF UKFd)
endif()

