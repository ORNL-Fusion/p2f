OBJ = obj
SRC = src
MOD = mod

ifeq ($(MACHINE),dlghp)
	HOME = /home/dg6/code
	F90	= ${HOME}/openmpi/gnu_64/bin/mpif90 -J${MOD}/
	WARN = -Wall -march=core2 -O3 -fbounds-check
	DISLIN = -I ${HOME}/dislin/dislin64/gf -L ${HOME}/dislin/dislin64 -ldislin 
	#PGPLOT = -L ${HOME}/pgplot -lpgplot -lX11
	NETCDF = -I ${HOME}/netcdf/netcdf_gnu64/include -L ${HOME}/netcdf/netcdf_gnu64/lib -lnetcdf
endif
ifeq ($(MACHINE),franklin)
	F90 = ftn -module ${MOD}/
	WARN = 
endif
ifeq ($(MACHINE),jaguar)
	F90 = ftn -module ${MOD}/
	WARN = 
	NETCDF = ${NETCDF_FLIB}
endif

ifeq ($(HOSTNAME),benten.gat.com)

        #F90 = ${HOME}/openmpi/openmpi.gcc/bin/mpif90 -J${MOD}/
        #WARN = 
        #NETCDF = -I ${HOME}/netcdf/netcdf.gcc/include -L ${HOME}/netcdf/netcdf.gcc/lib -lnetcdf

        F90 = ${HOME}/openmpi/openmpi.pgi/bin/mpif90 -module ${MOD}/
        WARN = 
        NETCDF = -I ${HOME}/netcdf/netcdf.pgi/include -L ${HOME}/netcdf/netcdf.pgi/lib -lnetcdf

        MACHINE = ${HOSTNAME}

endif

EXEC = xp2f.${MACHINE}

OBJECTS = ${OBJ}/eqdsk.o ${OBJ}/dlg.o ${OBJ}/fitpack.o ${OBJ}/gc_terms.o ${OBJ}/interp.o ${OBJ}/gc_integrate.o ${OBJ}/read_particle_list.o ${OBJ}/rzvv_grid.o ${OBJ}/init_mpi.o ${OBJ}/write_f_rzvv.o ${OBJ}/read_namelist.o ${OBJ}/beselI.o ${OBJ}/constants.o

ifeq (${MACHINE},dlghp)
p2f: ${SRC}/p2f.f90 ${OBJECTS}
	${F90} ${SRC}/p2f.f90 -o ${EXEC} ${OBJECTS} ${WARN} ${NETCDF} ${DISLIN}
else
p2f: ${SRC}/p2f.f90 ${OBJECTS}
	${F90} ${SRC}/p2f.f90 -o ${EXEC} ${OBJECTS} ${WARN} ${NETCDF}
endif

${MOD}/write_f_rzvv.mod: ${SRC}/write_f_rzvv.f90 ${OBJ}/write_f_rzvv.o
${OBJ}/write_f_rzvv.o: ${SRC}/write_f_rzvv.f90 
	${F90} -c ${SRC}/write_f_rzvv.f90 -o ${OBJ}/write_f_rzvv.o ${WARN} ${NETCDF}

${MOD}/gc_integrate.mod: ${SRC}/gc_integrate.F90 ${OBJ}/gc_integrate.o 

ifeq ($(MACHINE),dlghp)
${OBJ}/gc_integrate.o: ${SRC}/gc_integrate.F90 ${MOD}/interp.mod ${MOD}/gc_terms.mod ${MOD}/rzvv_grid.mod ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/gc_integrate.F90 -DUSE_DISLIN -o ${OBJ}/gc_integrate.o ${WARN} ${DISLIN}
else
${OBJ}/gc_integrate.o: ${SRC}/gc_integrate.F90 ${MOD}/interp.mod ${MOD}/gc_terms.mod ${MOD}/rzvv_grid.mod ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/gc_integrate.F90 -o ${OBJ}/gc_integrate.o ${WARN} 
endif

${MOD}/rzvv_grid.mod: ${SRC}/rzvv_grid.f90 ${OBJ}/rzvv_grid.o ${MOD}/read_namelist.mod
${OBJ}/rzvv_grid.o: ${SRC}/rzvv_grid.f90 ${MOD}/read_particle_list.mod
	${F90} -c ${SRC}/rzvv_grid.f90 -o ${OBJ}/rzvv_grid.o ${WARN}

${MOD}/read_particle_list.mod: ${SRC}/read_particle_list.f90 ${OBJ}/read_particle_list.o
${OBJ}/read_particle_list.o: ${SRC}/read_particle_list.f90 ${MOD}/init_mpi.mod ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/read_particle_list.f90 -o ${OBJ}/read_particle_list.o ${WARN} ${NETCDF}

${MOD}/init_mpi.mod: ${SRC}/init_mpi.f90 ${OBJ}/init_mpi.o
${OBJ}/init_mpi.o: ${SRC}/init_mpi.f90 ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/init_mpi.f90 -o ${OBJ}/init_mpi.o ${WARN}

${MOD}/interp.mod: ${SRC}/interp.f90 ${OBJ}/interp.o 
${OBJ}/interp.o: ${SRC}/interp.f90 ${MOD}/gc_terms.mod ${MOD}/eqdsk.mod
	${F90} -c ${SRC}/interp.f90 -o ${OBJ}/interp.o ${WARN}

${MOD}/gc_terms.mod: ${SRC}/gc_terms.f90 ${OBJ}/gc_terms.o  
${OBJ}/gc_terms.o: ${SRC}/gc_terms.f90 ${MOD}/eqdsk.mod ${MOD}/constants.mod
	${F90} -c ${SRC}/gc_terms.f90 -o ${OBJ}/gc_terms.o ${WARN}

${MOD}/eqdsk.mod: ${SRC}/eqdsk.f90 ${OBJ}/eqdsk.o  
${OBJ}/eqdsk.o: ${SRC}/eqdsk.f90 ${MOD}/dlg.mod 
	${F90} -c ${SRC}/eqdsk.f90 -o ${OBJ}/eqdsk.o ${WARN}

${MOD}/constants.mod: ${SRC}/constants.f90 ${OBJ}/constants.o
${OBJ}/constants.o: ${SRC}/constants.f90 ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/constants.f90 -o ${OBJ}/constants.o ${WARN}

${MOD}/read_namelist.mod: ${SRC}/read_namelist.f90 ${OBJ}/read_namelist.o
${OBJ}/read_namelist.o: ${SRC}/read_namelist.f90
	${F90} -c ${SRC}/read_namelist.f90 -o ${OBJ}/read_namelist.o ${WARN}

${MOD}/dlg.mod: ${SRC}/dlg.f90 ${OBJ}/dlg.o
${OBJ}/dlg.o: ${SRC}/dlg.f90
	${F90} -c ${SRC}/dlg.f90 -o ${OBJ}/dlg.o ${WARN} ${NETCDF}

${OBJ}/beselI.o: ${SRC}/bessel/beselI.f90
	${F90} -c ${SRC}/bessel/beselI.f90 -o ${OBJ}/beselI.o ${WARN}

${OBJ}/fitpack.o: ${SRC}/fitpack.f
	${F90} -c ${SRC}/fitpack.f -o ${OBJ}/fitpack.o

clean:
	rm ${OBJ}/*.o ${MOD}/*.mod ${EXEC}  
