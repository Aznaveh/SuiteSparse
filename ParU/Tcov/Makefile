default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CLIB = $(LDFLAGS) -L../../lib -lumfpack -lsuitesparseconfig -lcholmod -lamd \
        -lcolamd $(LIB_WITH_PARTITION)  $(LDLIBS)
FLIB = $(LAPACK) $(BLAS)


I = -I../../include  -I ../Include -I../../UMFPACK/Source -I../../AMD/Include


C = $(CXX) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) $(CONFIG_GPU) $(I) \
	$(CHOLMOD_CONFIG)

LIBS = $(CLIB) $(FLIB) $(TBB) $(GPULIB) 
#-lasan

FLAG = -Wno-write-strings 
#-O0 -fsanitize=address -g


print-%  : ; @echo $* = $($*)

# this test requires gcc
CC = gcc
CXX= g++
#BLAS = -lrefblas
LAPACK = -llapack

# for statement coverage (with gcov; see go) and picky compiler warnings
CF = -pg \
	-O0 -g -fprofile-arcs -ftest-coverage \
	-Wall -W -Wshadow -Winline -Wno-unused-parameter \
	-Wredundant-decls -Wdisabled-optimization \
	-fexceptions -fopenmp

# Using an optimized BLAS can cause problems in Valgrind
# FLIB = -llapack_plain -lblas_plain -lgfortran
# FLIB = -llapack_plain -lblas_plain -lgfortran -lg2c
  FLIB = $(LAPACK) $(BLAS)


all: parutest
	


#-------------------------------------------------------------------------------
INC = ../Include/Parallel_LU.hpp

OBJ = \
    paru_mem.o \
    paru_analyze.o \
	paru_init_rowFronts.o\
	paru_tuples.o\
	paru_front.o\
	paru_factorize.o\
	paru_create_element.o\
	paru_numeric_assemble.o\
	paru_trsm.o\
	paru_dgemm.o\
	paru_print.o\
	paru_init_rel.o\
	paru_update_rel_ind.o\
	paru_update_rowDeg.o\
	paru_finalize.o\
	paru_write.o\
	paru_cumsum.o

$(OBJ): $(INC)

I = -I ../../UMFPACK/Source -I../../UMFPACK/Include -I ../../CHOLMOD/Include\
	-I ../../AMD/Include  -I ../../SuiteSparse_config -I../Include
FLAG = -Wno-write-strings 
C = $(CXX) $(CF)  $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I) $(FLAG)
#-------------------------------------------------------------------------------
parutest: libraries $(INC) $(OBJ) 
	$(C)  umf_test.cpp -o umfout $(OBJ) $(LIBS)
#-------------------------------------------------------------------------------

paru_mem.o: ../Source/paru_mem.cpp
	$(C) -c $<


paru_analyze.o: ../Source/paru_analyze.cpp
	$(C) -c $<
#	- ln -s ../Source/paru_analyze.cpp
#	$(C) -c paru_analyze.cpp

paru_init_rowFronts.o: ../Source/paru_init_rowFronts.cpp
	$(C) -c $<

paru_tuples.o: ../Source/paru_tuples.cpp
	$(C) -c $<

#paru_assemble.o: ../Source/paru_assemble.cpp
#	$(C) -c $<

paru_front.o: ../Source/paru_front.cpp
	$(C) -c $<

paru_factorize.o: ../Source/paru_factorize.cpp
	$(C) -c $<

paru_create_element.o: ../Source/paru_create_element.cpp
	$(C) -c $<

paru_numeric_assemble.o: ../Source/paru_numeric_assemble.cpp
	$(C) -c $<

paru_trsm.o: ../Source/paru_trsm.cpp
	$(C) -c $<

paru_dgemm.o: ../Source/paru_dgemm.cpp
	$(C) -c $<

paru_print.o: ../Source/paru_print.cpp
	$(C) -c $<

paru_init_rel.o: ../Source/paru_init_rel.cpp
	$(C) -c $<

paru_update_rel_ind.o: ../Source/paru_update_rel_ind.cpp
	$(C) -c $<

paru_write.o: ../Source/paru_write.cpp
	$(C) -c $<

paru_update_rowDeg.o: ../Source/paru_update_rowDeg.cpp
	$(C) -c $<

paru_finalize.o: ../Source/paru_finalize.cpp
	$(C) -c $<

paru_cumsum.o: ../Source/paru_cumsum.cpp
	$(C) -c $<

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# libraries compiled without test coverage
#-------------------------------------------------------------------------------

libraries:
	( cd ../.. && $(MAKE) metis )
	( cd ../../SuiteSparse_config ; $(MAKE) library TCOV=yes )
	( cd ../../AMD && $(MAKE) library TCOV=yes )
	( cd ../../COLAMD && $(MAKE) library TCOV=yes )
	( cd ../../CCOLAMD && $(MAKE) library TCOV=yes )
	( cd ../../CAMD && $(MAKE) library TCOV=yes )
	( cd ../../CHOLMOD && $(MAKE) library TCOV=yes )
	( cd ../../UMFPACK && $(MAKE) library TCOV=yes )

#-------------------------------------------------------------------------------
#


#-------------------------------------------------------------------------------
run:
	./umfout < ../Matrix/b1_ss.mtx  	#7x7
#	./umfout < ../Matrix/ParUTst/tmp.mtx

VGFLAGS =  --leak-check=yes  --log-file="valgrindlogfile.log"\
		   --track-origins=yes\
		   #--show-leak-kinds=all

			#comment leak-check=yes or full

val:

	valgrind	$(VGFLAGS)	./umfout < ../Matrix/ParUTst/tmp.mtx
#	valgrind	$(VGFLAGS)	./umfout < ../Matrix/west0067.mtx 
#	valgrind	$(VGFLAGS)	./umfout < ../Matrix/b1_ss.mtx

	
debug:
	$(C) $(FLAG) -g umf_test.cpp -o umfout $(LIBS)
	gdb ./umfout < ../Matrix/c32.mtx		#3x2 	 #OK

gprofile:
	$(C) $(FLAG) -no-pie -O0 -g -pg umf_test.cpp -o gprof.o $(LIBS)
	./gprof.o < ../Matrix/west0067.mtx 
	gprof gprof.o > analysis.txt


vprofile:

	$(C) $(FLAG) -g umf_test.cpp -o vprof $(OBJ) $(LIBS)
	valgrind	--tool=callgrind -v ./vprof< ../Matrix/west0067.mtx
#	valgrind	--tool=callgrind -v ./vprof< ../Matrix/ParUTst/tmp.mtx

purge: distclean

distclean: clean
	- $(RM) umfout 
	- $(RM) analysis.txt tcov gprof
	- $(RM) gprof.o vprof gmon.out callgrind.out.*
	- $(RM) paneltest
	- $(RM) *.dot pfile tfile
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

az:
	(cd ..; $(MAKE) )
	($(MAKE) )
	(cd Tcov; $(MAKE) run )
