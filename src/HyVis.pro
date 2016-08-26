#-------------------------------------------------
#
# Project created by QtCreator 2015-06-18T15:11:43
#
#-------------------------------------------------

QT       += core gui opengl network

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = HyVis
TEMPLATE = app

LIBS += -L/opt/local/lib -lhdf5

INCLUDEPATH += /opt/local/include /opt/local/ /opt/local/include/mpich-gcc49 


SOURCES += main.cpp\
        controlwindow.cpp \
    viewer.cpp \
    geometry.cpp \
    filedata.cpp \
    querywindow.cpp \
    qcustomplot.cpp \
    plot2dviewer.cpp \
    config.cpp

HEADERS  += controlwindow.h \
    viewer.h \
    geometry.h \
    filedata.h \
    querywindow.h \
    qcustomplot.h \
    plot2dviewer.h \
    colormaps.h \
    config.h

FORMS    += controlwindow.ui \
    querywindow.ui \
    plot2dviewer.ui \
    config.ui

RESOURCES += \
    Shaders.qrc

OTHER_FILES += \
    FWhite.fsh


# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc )
QMAKE_LFLAGS += $$system(mpicxx )
QMAKE_CXXFLAGS += $$system(mpicxx ) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx ) -DMPICH_IGNORE_CXX_SEEK -std=c++11

