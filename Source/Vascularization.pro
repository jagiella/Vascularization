 INCLUDEPATH += ../../NiX
 DEFINES += USE_GUI
 unix:LIBS += ~/lib/libNiX.so
 macx-g++:LIBS -= ~/lib/libNiX.so
 macx-g++:LIBS += /usr/local/lib/libNiX.dylib
 LFLAGS   += -lm
# CFLAGS   += -fopenmp -O3 -Wall -Wextra
 CXXFLAGS += -Wall -Wextra -ansi
# !!ATTENTION!! OPENMP INTERFERRES WITH QT AND QTHREAD
 QMAKE_LFLAGS   += -lm
# QMAKE_CFLAGS   += -fopenmp -g -Wall -Wextra
 QMAKE_CXXFLAGS += -Wall -Wextra -ansi
 QT          += opengl

# SOURCE SETTINGS

 HEADERS     = glwidget.hpp \
               helper.hpp \
               widget.hpp \
               window.hpp \
               Perfusion.hpp \
               VesselGraph.hpp

 SOURCES     = glwidget.cpp \
               helper.cpp \
               main.cpp \
               widget.cpp \
               window.cpp \
               VesselGraph.cpp \
               Perfusion.cpp
               
win32:RC_FILE    = Vascularization.rc
macx-g++:RC_FILE = Vascularization.icns



# OUTPUT SETTINGS

# object file directory
OBJECTS_DIR = $$PWD/../GUI/obj

# executable directory
DESTDIR = $$PWD/../GUI/bin

# moc file directory
MOC_DIR = $$PWD/../GUI/moc

 # install
 target.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
 sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS 2dpainting.pro
 sources.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
 INSTALLS += target sources

 symbian: include($$QT_SOURCE_TREE/examples/symbianpkgrules.pri)
