

#
# If pkg-config isn't installed on your system, comment the following lines and
# set the fasterac flags as indicated with your own paths:
#
FASTERAC_CFLAGS = -I/usr/include
FASTERAC_LIBS   = -L/usr/lib -lfasterac
#
#FASTERAC_CFLAGS = $(shell pkg-config --cflags libfasterac)
#FASTERAC_LIBS   = $(shell pkg-config --libs   libfasterac)
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS   = $(shell root-config --glibs)
#


CC     = g++
CFLAGS = ${ROOT_CFLAGS} ${FASTERAC_CFLAGS} -Wall -lSpectrum
LIBS   = ${FASTERAC_LIBS} ${ROOT_LIBS} 
SRCEXE = $(shell ls *.C)
EXE    = $(SRCEXE:.C=)


all : $(EXE)

$(EXE): $(SRCEXE)
	${CC} $@.C -o $@ ${INCFLAGS} ${CFLAGS} ${LIBS} -w

clean :
	rm -f *.o
	rm -f $(EXE)


