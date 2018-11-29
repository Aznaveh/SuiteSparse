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
#	./testazny < ../Matrix/arrow.mtx    #Huge
#	./testazny < ../Matrix/ash219.mtx   #OK
#	./testazny < ../Matrix/b1_ss.mtx  	#OK 
#	./testazny < ../Matrix/bfwa62.mtx   #OK 
#	./testazny < ../Matrix/c2.mtx		#OK #2x2 4nz
#	./testazny < ../Matrix/c32.mtx		#OK #3x2 
#	./testazny < ../Matrix/cage3.mtx    #OK
#	./testazny < ../Matrix/cage5.mtx    #OK   #Good sample
#	./testazny < ../Matrix/Franz6_id1959_aug.mtx #Huge
#	./testazny < ../Matrix/GD98_a.mtx   #OK		 #Empty Row
#	./testazny < ../Matrix/Groebner_id2003_aug.mtx #HUGE not good for MATLAB
#	./testazny < ../Matrix/LFAT5.mtx	#OK
#	./testazny < ../Matrix/lfat5b.mtx   #OK 
#
#	./testazny < ../Matrix/lp_e226_transposed.mtx   #SEG FAULT
#
#	./testazny < ../Matrix/lp_e226.mtx 		#OK  	# m >=n 
#	./testazny < ../Matrix/lp_share1b.mtx 	#OK  	# m >=n 
#	./testazny < ../Matrix/lpi_galenet.mtx 	#OK 	# m >=n 
#	./testazny < ../Matrix/lpi_itest6.mtx  	#OK  	# m >=n 
#	./testazny < ../Matrix/n3c4-b4.mtx    	#OK   	# m >=n 
#	./testazny < ../Matrix/permuted_b1_ss.mtx  #OK
#	./testazny < ../Matrix/problem.mtx    #OK  # m >=n
#	./testazny < ../Matrix/pwr01b.mtx	  #OK  #Good sample
#	./testazny < ../Matrix/r2.mtx         #OK  # m >=n failure
#	./testazny < ../Matrix/Ragusa16.mtx   #OK #EMPTY row
#	./testazny < ../Matrix/s32.mtx		  #OK
#	./testazny < ../Matrix/Tina_AskCal.mtx       #OK #EMPTY col
#	./testazny < ../Matrix/Tina_AskCal_perm.mtx  #OK #EMPTY col
#	./testazny < ../Matrix/west0067.mtx 	#OK #Goood sample
#	./testazny < ../Matrix/young1c.mtx     	#OK  
#	./testazny < ../Matrix/ParUTst/az11_11.mtx   #OK
#	./testazny < ../Matrix/ParUTst/az88_2.mtx	 #OK
#	./testazny < ../Matrix/ParUTst/bcsstk01.mtx  #OK
#	./testazny < ../Matrix/ParUTst/bcsstk03.mtx  #NO  #Memory problem
#	./testazny < ../Matrix/ParUTst/bcsstk22.mtx  #OK
#	./testazny < ../Matrix/ParUTst/bcsstm01.mtx	  		#OK
#	./testazny < ../Matrix/ParUTst/bcsstm02.mtx	  		#OK
#	./testazny < ../Matrix/ParUTst/bcsstm05.mtx	  		#OK  		
#	./testazny < ../Matrix/ParUTst/bcsstm06.mtx	  		#OK  
#	./testazny < ../Matrix/ParUTst/bcsstm19.mtx		 	#OK  
#	./testazny < ../Matrix/ParUTst/bcsstm20.mtx		 	#OK  
#	./testazny < ../Matrix/ParUTst/ex5.mtx  	#OK  #Good sample
#	./testazny < ../Matrix/ParUTst/mesh1e1.mtx  #OK
#	./testazny < ../Matrix/ParUTst/LF10/LF10.mtx #OK  #Good sample
#	./testazny < ../Matrix/ParUTst/mesh1em1/mesh1em1.mtx #OK
#	./testazny < ../Matrix/ParUTst/Trefethen_20/Trefethen_20.mtx
#	./testazny < ../Matrix/ParUTst/Trefethen_20b/Trefethen_20b.mtx
#	./testazny < ../Matrix/ParUTst/curtis54/curtis54.mtx
	./testazny < ../Matrix/ParUTst/nos1/nos1.mtx #NO  #Memory problem
#	./testazny < ../Matrix/ParUTst/nos4/nos4.mtx #OK
#	./testazny < ../Matrix/ParUTst/494_bus/494_bus.mtx
#	./testazny < ../Matrix/ParUTst/mesh1em6/mesh1em6.mtx



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
