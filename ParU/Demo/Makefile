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

FLAG = -Wno-write-strings

all:
	$(C) $(FLAG) testazny.cpp -o testazny $(LIBS)
	

run:
#	./testazny < ../Matrix/a0.mtx  		#OK #0x0
#	./testazny < ../Matrix/a04.mtx 		#OK #0x4
#	./testazny < ../Matrix/a1.mtx 		#OK #Empty row
#	./testazny < ../Matrix/a2.mtx 		#OK
#	./testazny < ../Matrix/a4.mtx 		#OK #Empty row
#	./testazny < ../Matrix/arrow.mtx  #Not suitalbe for MATLAB  #Huge
#	./testazny < ../Matrix/ash219.mtx   #NO
#	./testazny < ../Matrix/b1_ss.mtx  	#OK #not positive definite
#	./testazny < ../Matrix/bfwa62.mtx   #NO m<n in the middle
#	./testazny < ../Matrix/c2.mtx		#OK #2x2 4nz
#	./testazny < ../Matrix/c32.mtx		#OK #3x2 
#	./testazny < ../Matrix/cage5.mtx    #NO   # m >=n in one of the fronts
#	./testazny < ../Matrix/Franz6_id1959_aug.mtx #NO SEG FAULT
#	./testazny < ../Matrix/GD98_a.mtx   #NO		 #Structuraly rank defecient
#	./testazny < ../Matrix/Groebner_id2003_aug.mtx #HUGE not good for MATLAB
#	./testazny < ../Matrix/LFAT5.mtx	#OK for err=1e-10
#	./testazny < ../Matrix/lfat5b.mtx   #OK for err=2e-16
#	./testazny < ../Matrix/lp_e226_transposed.mtx   #NO SEG FAULT
#	./testazny < ../Matrix/lp_e226.mtx #NO  # m >=n in one of the fronts
#	./testazny < ../Matrix/lp_share1b.mtx #NO  # m >=n assertion failure
#	./testazny < ../Matrix/lpi_galenet.mtx #NO # m >=n assertion failure
#	./testazny < ../Matrix/lpi_itest6.mtx  #NO # m >=n assertion failure
#	./testazny < ../Matrix/n3c4-b4.mtx    #NO  # m >=n assertion failure
#	./testazny < ../Matrix/permuted_b1_ss.mtx  #OK
#	./testazny < ../Matrix/problem.mtx    #NO  # m >=n assertion failure
#	./testazny < ../Matrix/pwr01b.mtx	  #NO  # m >=n assertion failure
#	./testazny < ../Matrix/r2.mtx         #NO  # m >=n assertion failure
#	./testazny < ../Matrix/Ragusa16.mtx   #OK #EMPTY row
#	./testazny < ../Matrix/s32.mtx		  #OK
#	./testazny < ../Matrix/Tina_AskCal.mtx       #OK #EMPTY col
#	./testazny < ../Matrix/Tina_AskCal_perm.mtx  #OK #EMPTY col
#	./testazny < ../Matrix/west0067.mtx 	#NO  # m >=n assertion failure
#	./testazny < ../Matrix/young1c.mtx     	#NO  # m >=n assertion failure
#	
#	./testazny < ../Matrix/ParUTst/az11_11.mtx  		#OK
	./testazny < ../Matrix/ParUTst/az88_2.mtx	  		#OK
#	./testazny < ../Matrix/ParUTst/bcsstk01.mtx  #NO  # m >=n assertion failure
#	./testazny < ../Matrix/ParUTst/bcsstk03.mtx  #NO  # m >=n assertion failure
#	./testazny < ../Matrix/ParUTst/bcsstk22.mtx  #NO  # m >=n assertion failure
#	./testazny < ../Matrix/ParUTst/bcsstm01.mtx	  		#OK
#	./testazny < ../Matrix/ParUTst/bcsstm02.mtx	  		#OK
#	./testazny < ../Matrix/ParUTst/bcsstm05.mtx	  		#OK  		
#	./testazny < ../Matrix/ParUTst/bcsstm06.mtx	  		#OK  
#	./testazny < ../Matrix/ParUTst/bcsstm19.mtx		 	#OK  
#	./testazny < ../Matrix/ParUTst/bcsstm20.mtx		 	#OK  
#	./testazny < ../Matrix/ParUTst/ex5.mtx  	#NO  # m >=n assertion failure
#	./testazny < ../Matrix/ParUTst/mesh1e1.mtx  #NO SEG FAULT + ALL
#	./testazny < ../Matrix/ParUTst/LF10/LF10.mtx #NO  # m >=n assertion failure
#	./testazny < ../Matrix/ParUTst/mesh1em1/mesh1em1.mtx
#


VGFLAGS = --tool=memcheck --leak-check=full --show-reachable=yes -v\
		  --track-origins=yes --num-callers=20 \
		  --track-fds=yes --log-file="valgrindlogfile.log"

val:
	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lp_share1b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/LFAT5.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a2.mtx 
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/b1_ss.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lfat5b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/west0067.mtx 
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
#
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/s32.mtx
#
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/c32.mtx 
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/Groebner_id2003_aug.mtx  
#
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a04.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lp_e226.mtx

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

az:
	(cd ..; $(MAKE) )
	(cd Demo; $(MAKE) run )
