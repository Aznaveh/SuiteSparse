#===============================================================================
# ParU/Tcov/Makefile
#===============================================================================


default: go

valgrind: vgo

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CLIB = $(LDFLAGS) -L../../lib -lumfpack -lsuitesparseconfig -lcholmod -lamd \
        -lcolamd $(LIB_WITH_PARTITION)  $(LDLIBS)
FLIB = $(LAPACK) $(BLAS)


I = -I../../include  -I ../Include -I../../UMFPACK/Source -I../../AMD/Include


C = $(CXX) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) $(CONFIG_GPU) $(I) \
	$(CHOLMOD_CONFIG)

LIBS = $(CLIB) $(FLIB) $(TBB) $(GPULIB) 
#-lasan


print-%  : ; @echo $* = $($*)

# for statement coverage (with gcov; see go) and picky compiler warnings
# Aznaveh: I took out -Wredundant-decls because I overloade new and delete
CF = -pg \
	-O0 -g -fprofile-arcs -ftest-coverage \
	-Wall -W -Wshadow -Winline -Wno-unused-parameter \
	 -Wdisabled-optimization \
	-fexceptions -fopenmp 

# Using an optimized BLAS can cause problems in Valgrind
# FLIB = -llapack_plain -lblas_plain -lgfortran
# FLIB = -llapack_plain -lblas_plain -lgfortran -lg2c
  FLIB = $(LAPACK) $(BLAS)


all: parutest
	


#-------------------------------------------------------------------------------
INC = ../Include/ParU.hpp 
C11 = -std=c++11
CC11 = -std=c11


OBJ = \
    paru_mem.o \
    paru_memset.o \
    paru_memcpy.o \
    paru_heap.o \
    paru_hash.o \
    paru_pivotal.o \
    paru_analyze.o \
	paru_init_rowFronts.o\
	paru_tuples.o\
	paru_front.o\
	paru_factorize.o\
	paru_fs_factorize.o\
	paru_create_element.o\
	paru_assemble_row2U.o\
	paru_trsm.o\
	paru_dgemm.o\
	paru_print.o\
	paru_init_rel.o\
	paru_update_rowDeg.o\
	paru_update_rel_ind.o\
	paru_write.o\
	paru_cumsum.o\
	paru_intersection.o\
	paru_bin_search.o\
	paru_assemble.o\
	paru_prior_assemble.o \
	paru_perm.o \
	paru_lsolve.o \
	paru_usolve.o \
	paru_solve.o \
	paru_residual.o \
	paru_backward.o \
	paru_gaxpy.o \
	paru_norms.o \
	paru_Diag_update.o \
	paru_full_summed.o \
	paru_tasked_dgemm.o \
	paru_tasked_trsm.o \


$(OBJ): $(INC)

I = -I ../../UMFPACK/Source -I../../UMFPACK/Include -I ../../CHOLMOD/Include\
	-I ../../AMD/Include  -I ../../SuiteSparse_config -I../Include

FLAG = -Wno-write-strings


#C = $(CXX) $(CF) $(SPQR_CONFIG) $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I)
C = $(CXX) $(CF)  $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I) $(FLAG) $(C11)
C2 = $(CC) $(CF)  $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I) $(FLAG) $(CC11)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
parutest: libraries paru_test.cpp $(INC) $(OBJ) 
	$(C)  paru_test.cpp -o x_test $(OBJ) $(LIBS) -lm

go: parutest
	- ./x_test < ../Matrix/west0067.mtx > parutest_out.txt
	- ./cov

vgo: parutest
	- ./x_test < ../Matrix/west0067.mtx > parutest_out.txt
	#- ./cov


#-------------------------------------------------------------------------------
paru_mem.o: ../Source/paru_mem.cpp
	$(C) -c $<

paru_memset.o: ../Source/paru_memset.cpp
	$(C) -c $<

paru_memcpy.o: ../Source/paru_memcpy.cpp
	$(C) -c $<

paru_heap.o: ../Source/paru_heap.cpp
	$(C) -c $<

paru_hash.o: ../Source/paru_hash.cpp
	$(C) -c $<
	
paru_pivotal.o: ../Source/paru_pivotal.cpp
	$(C) -c $<

paru_analyze.o: ../Source/paru_analyze.cpp
	$(C) -c $<

paru_init_rowFronts.o: ../Source/paru_init_rowFronts.cpp
	$(C) -c $<

paru_tuples.o: ../Source/paru_tuples.cpp
	$(C) -c $<

paru_front.o: ../Source/paru_front.cpp
	$(C) -c $<

paru_factorize.o: ../Source/paru_factorize.cpp
	$(C) -c $<

paru_fs_factorize.o: ../Source/paru_fs_factorize.cpp
	$(C) -c $<

paru_create_element.o: ../Source/paru_create_element.cpp
	$(C) -c $<

paru_assemble_row2U.o: ../Source/paru_assemble_row2U.cpp
	$(C) -c $<

paru_trsm.o: ../Source/paru_trsm.cpp
	$(C) -c $<

paru_dgemm.o: ../Source/paru_dgemm.cpp
	$(C) -c $<

paru_print.o: ../Source/paru_print.cpp
	$(C) -c $<

paru_init_rel.o: ../Source/paru_init_rel.cpp
	$(C) -c $<

paru_write.o: ../Source/paru_write.cpp
	$(C) -c $<

paru_update_rowDeg.o: ../Source/paru_update_rowDeg.cpp
	$(C) -c $<

paru_update_rel_ind.o: ../Source/paru_update_rel_ind.cpp
	$(C) -c $<

paru_cumsum.o: ../Source/paru_cumsum.cpp
	$(C) -c $<

paru_intersection.o: ../Source/paru_intersection.cpp
	$(C) -c $<

paru_bin_search.o: ../Source/paru_bin_search.cpp
	$(C) -c $<

paru_prior_assemble.o:  ../Source/paru_prior_assemble.cpp
	$(C) -c $<

paru_assemble.o: ../Source/paru_assemble.cpp
	$(C) -c $<

paru_perm.o: ../Source/paru_perm.cpp
	$(C) -c $<

paru_lsolve.o: ../Source/paru_lsolve.cpp
	$(C) -c $<

paru_usolve.o: ../Source/paru_usolve.cpp
	$(C) -c $<

paru_solve.o: ../Source/paru_solve.cpp
	$(C) -c $<

paru_residual.o: ../Source/paru_residual.cpp
	$(C) -c $<

paru_backward.o: ../Source/paru_backward.cpp
	$(C) -c $<

paru_gaxpy.o: ../Source/paru_gaxpy.cpp
	$(C) -c $<

paru_norms.o: ../Source/paru_norms.cpp
	$(C) -c $<

paru_Diag_update.o: ../Source/paru_Diag_update.cpp
	$(C) -c $<

paru_full_summed.o: ../Source/paru_full_summed.cpp
	$(C) -c $<

paru_tasked_dgemm.o: ../Source/paru_tasked_dgemm.cpp
	$(C) -c $<

paru_tasked_trsm.o: ../Source/paru_tasked_trsm.cpp
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

purge: distclean

distclean: clean
	- $(RM) x_test
	- $(RM) parutest_out.txt
	- $(RM) analysis.txt tcov gprof
	- $(RM) gprof.o vprof gmon.out callgrind.out.*
	- $(RM) paneltest
	- $(RM) *.dot pfile tfile
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)
