TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    channelindexpair.cpp \
    symblock.cpp \
    ccd.cpp \
    qstate.cpp \
    qdotspbasis.cpp \
    Coulomb_Functions.cpp \
    generalspclass.cpp \
    channelset.cpp \
    channel.cpp \
    qdotchannelset.cpp \
    qdotHFbasis.cpp

HEADERS += \
    generalspclass.h \
    channel.h \
    channelindexpair.h \
    symblock.h \
    channelset.h \
    ccd.h \
    qstate.h \
    qdotspbasis.h \
    Coulomb_Functions.hpp \
    qdotchannelset.h \
    qdotHFbasis.h

