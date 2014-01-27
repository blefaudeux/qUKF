#-------------------------------------------------
#
# Project created by QtCreator 2014-01-27T16:36:08
#
#-------------------------------------------------

TARGET = libUKF
TEMPLATE = lib

DEFINES += LIBUKF_LIBRARY

INCLUDEPATH += /usr/include \
        /usr/include/eigen3 \
        headers

SOURCES += src/eigen_tools.cpp \
    src/sigma_point.cpp \
    src/sigma_q_point.cpp \
	  src/statistic_tools.cpp \
	  src/unscented_KF.cpp

HEADERS += headers/eigen_tools.h \
    headers/sigma_point.h \
    headers/sigma_q_point.h \
    headers/statistic_tools.h \
    headers/unscented_KF.h \
	  headers/def.h

