PARALLEL:=1

OBJ = obj
SRC := src
MOD := mod

FFLAGS :=
F90FLAGS :=  
INCLUDES :=
LINKS :=

CPP_DIRECTIVES :=
CPP_DIRECTIVES += -DDEBUG_LEVEL=0
CPP_DIRECTIVES += -DPARALLEL=${PARALLEL}

ThisMachine := $(shell uname -n)
include Makefile.$(ThisMachine)

EXEC = xp2f.${ThisMachine}

OBJECTS = $(wildcard ${OBJ}/*.o)

${EXEC}: ${SRC}/p2f.F90  
	${F90} ${F90FLAGS} ${SRC}/p2f.F90 -o ${EXEC} ${INCLUDES} ${LINKS} ${OBJECTS} ${CPP_DIRECTIVES}

${OBJ}/%.o: ${SRC}/%.f
	${F77} -c ${FFLAGS} ${INCLUDES} $< -o $@ 

${OBJ}/%.o: ${SRC}/%.f90
	${F90} -c ${F90FLAGS} ${INCLUDES} $< -o $@ 

${OBJ}/%.o: ${SRC}/%.F90
	${F90} -c ${F90FLAGS} ${INCLUDES} $< -o $@ ${CPP_DIRECTIVES}

include Makefile.deps

clean:
	rm ${OBJ}/*.o ${MOD}/*.mod ${EXEC}  
