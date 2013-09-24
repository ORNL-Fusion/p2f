OBJ = obj
SRC := src
MOD := mod

FFLAGS :=
F90FLAGS :=
INCLUDES :=
LINKS :=

ThisMachine := $(shell uname -n)
include Makefile.$(ThisMachine)

EXEC = xp2f.${ThisMachine}

OBJECTS = $(wildcard ${OBJ}/*.o)

${EXEC}: ${SRC}/p2f.f90  
	${F90} ${F90FLAGS} ${SRC}/p2f.f90 -o ${EXEC} ${INCLUDES} ${LINKS} ${OBJECTS} 

${OBJ}/%.o: ${SRC}/%.f
	${F77} -c ${FFLAGS} ${INCLUDES} $< -o $@ 

${OBJ}/%.o: ${SRC}/%.f90
	${F90} -c ${F90FLAGS} ${INCLUDES} $< -o $@ 

${OBJ}/%.o: ${SRC}/%.F90
	${F90} -c ${F90FLAGS} ${INCLUDES} $< -o $@ 

include Makefile.deps

clean:
	rm ${OBJ}/*.o ${MOD}/*.mod ${EXEC}  
