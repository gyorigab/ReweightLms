#-------------------------------------------------
#
# Project created by QtCreator 2012-12-21T21:26:50
#
#-------------------------------------------------

QT       -= core

QT       -= gui

TARGET = lms 
CONFIG   += console
CONFIG   -= app_bundle

#CONFIG += c++11

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS   += -fopenmp

TEMPLATE = app

LIBS += /usr/local/pgsql/lib/libpq.a

INCLUDEPATH += /usr/local/pgsql/include/
INCLUDEPATH += ../src/matvec/


HEADERS +=  matvec/array.h\
            matvec/bandmat.h\
            matvec/choldec.h\
            matvec/covmat.h\
            matvec/gso.h\
            matvec/hilbert.h\
            matvec/inderr.h\
            matvec/jacobian.h\
            matvec/matbase.h\
            matvec/mat.h\
            matvec/matvecbase.h\
            matvec/matvec.h\
            matvec/memrep.h\
            matvec/pinv.h\
            matvec/sez\
            matvec/sortvec.h\
            matvec/svd.h\
            matvec/symmat.h\
            matvec/transmat.h\
            matvec/transvec.h\
            matvec/vecbase.h\
            matvec/vec.h\
            smatrix.h \
            conjungate_gradients.h \
            operators.h \
            conjungate_gradients_nr.h \
            median.h \
            data.h \
            sample_model.h \
            general_cmp.h \
            settings.h \
            logger.h \
            db_results.h \
            db_settings.h \
            lms_settings.h \
            lms_exact.h

SOURCES +=   main.cpp \
             operators.cpp
