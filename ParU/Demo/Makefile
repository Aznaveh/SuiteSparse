default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CLIB = $(LDFLAGS) -L../../lib -lparu -lspqr -lsuitesparseconfig -lcholmod -lamd \
        -lcolamd $(LIB_WITH_PARTITION) $(LIB_WITH_SPQRGPU) $(LDLIBS)
FLIB = $(LAPACK) $(BLAS)

I_WITH_SPQRGPU = 

I = -I../../include $(I_WITH_SPQRGPU) -I ../Include


C = $(CXX) $(CF) $(SPQR_CONFIG) $(CONFIG_PARTITION) $(CONFIG_GPU) $(I) \
	$(CHOLMOD_CONFIG)

LIBS = $(CLIB) $(FLIB) $(TBB) $(GPULIB)

all:
	$(C) testazny.cpp -o testazny $(LIBS)

run:
#	./testazny < ../Matrix/GD01_b.mtx
#	./testazny < ../Matrix/LFAT5.mtx
#	./testazny < ../Matrix/a2.mtx #2x2
#	./testazny < ../Matrix/b1_ss.mtx
#	./testazny < ../Matrix/bcspwr01.mtx
#	./testazny < ../Matrix/bfwa62.mtx
#	./testazny < ../Matrix/GD06_theory.mtx
#	./testazny < ../Matrix/lfat5b.mtx
#	./testazny < ../Matrix/pwr01b.mtx
#	./testazny < ../Matrix/cage5.mtx
#
###  The following have some sort of structural problem
#	./testazny < ../Matrix/Tina_AskCal.mtx
#	./testazny < ../Matrix/problem.mtx
#	./testazny < ../Matrix/Tina_AskCal_perm.mtx
#	./testazny < ../Matrix/Ragusa16.mtx #Structuraly rank defecient
#	./testazny < ../Matrix/n3c4-b4.mtx
#	./testazny < ../Matrix/r2.mtx
#	./testazny < ../Matrix/s32.mtx
#	./testazny < ../Matrix/lpi_itest6.mtx
#	./testazny < ../Matrix/lpi_galenet.mtx
#	./testazny < ../Matrix/c32.mtx
#	./testazny < ../Matrix/Groebner_id2003_aug.mtx
#	./testazny < ../Matrix/c2.mtx #2x2 4nz
#	./testazny < ../Matrix/a1.mtx
#	./testazny < ../Matrix/a4.mtx
#	./testazny < ../Matrix/Ragusa16_pattern.mtx
#	./testazny < ../Matrix/GD98_a.mtx #Structuraly rank defecient
#	./testazny < ../Matrix/a04.mtx
#	./testazny < ../Matrix/a0.mtx #all zeros

##The following doesn't work yet
#	./testazny < ../Matrix/lp_e226.mtx #!
#	./testazny < ../Matrix/lp_e226_transposed.mtx #!
	./testazny < ../Matrix/lp_share1b.mtx #!
#	./testazny < ../Matrix/young1c.mtx #!
#	./testazny < ../Matrix/Franz6_id1959_aug.mtx #!
#	./testazny < ../Matrix/arrow.mtx #! some problem I don't know
#	./testazny < ../Matrix/ash219.mtx #! 219x85 #crash on this

VGFLAGS = --tool=memcheck --leak-check=yes --show-reachable=yes -v\
		  --num-callers=20 --track-fds=yes

val:
#	valgrind	$(VGFLAGS)	./testazny  < ../Matrix/GD01_b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/LFAT5.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a2.mtx #2x2
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/b1_ss.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/bcspwr01.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/bfwa62.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/GD06_theory.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lfat5b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/pwr01b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/west0067.mtx #a 67x67 matrix
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/cage5.mtx
#
#
#	valgrind	$(VGFLAGS)	./testazny <  ../Matrix/Tina_AskCal.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/problem.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/problem.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/Ragusa16.mtx 
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/n3c4-b4.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/r2.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lpi_itest6.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lpi_itest6.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lpi_galenet.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a1.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/c2.mtx #2x2 4nz
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a4.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/Ragusa16_pattern.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/GD98_a.mtx

#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/s32.mtx  #Done with some errors
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/c32.mtx  #Done with some errors
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/Groebner_id2003_aug.mtx  #Done with several errors

# Test these:
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a04.mtx

debug:
	g++ -g   -O3 -fexceptions -fPIC -fopenmp   -I../../include   aztest.cpp -o \
            dtest.out -L../../lib -lspqr -lsuitesparseconfig -lcholmod -lamd \
            -lcolamd -lccolamd -lcamd -lmetis  -lm -lrt \
            -Wl,-rpath=../../lib -llapack -lopenblas 


purge: distclean

distclean: clean
	- $(RM) testazny
	- $(RM) *.dot pfile tfile
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)
