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
#-lasan

FLAG = -Wno-write-strings 
#-O0 -fsanitize=address -g

all:
	$(C) $(FLAG) testazny.cpp -o testazny $(LIBS)
	

run:
#	./testazny < ../Matrix/a0.mtx  		 #0x0 #OK
#	./testazny < ../Matrix/a04.mtx 		 #0x4 #OK
#	./testazny < ../Matrix/a1.mtx 		 #Empty row #OK
#	./testazny < ../Matrix/a2.mtx 		#OK
#	./testazny < ../Matrix/a4.mtx 		 #Empty row #OK
#	./testazny < ../Matrix/arrow.mtx    #Huge
#	./testazny < ../Matrix/ash219.mtx   #219x85
#	./testazny < ../Matrix/b1_ss.mtx  	#7x7
#	./testazny < ../Matrix/bfwa62.mtx   #62x62
#	./testazny < ../Matrix/c2.mtx		#2x2 4nz #OK
#	./testazny < ../Matrix/c32.mtx		#3x2 	 #OK
#	./testazny < ../Matrix/cage3.mtx    #5x5
#	./testazny < ../Matrix/cage5.mtx    #37x37
#	./testazny < ../Matrix/Franz6_id1959_aug.mtx #Huge
#	./testazny < ../Matrix/GD98_a.mtx   #Empty Row OK
#	./testazny < ../Matrix/Groebner_id2003_aug.mtx #HUGE not good for MATLAB
#	./testazny < ../Matrix/LFAT5.mtx	#14x14
#	./testazny < ../Matrix/lfat5b.mtx   #14x14
#
#	./testazny < ../Matrix/ParUTst/tmp.mtx
#
#	./testazny < ../Matrix/lp_e226_transposed.mtx   #427x223
#	./testazny < ../Matrix/lp_e226.mtx 		# m >=n 
#	./testazny < ../Matrix/lp_share1b.mtx 	# m >=n 
#	./testazny < ../Matrix/lpi_galenet.mtx 	# m >=n 
#	./testazny < ../Matrix/lpi_itest6.mtx  	# m >=n 
#
#	./testazny < ../Matrix/ParUTst/lpi_galenet_transposed.mtx
#	./testazny < ../Matrix/ParUTst/lpk_b2_transposed.mtx
#	./testazny < ../Matrix/ParUTst/lp_kb2/lp_kb2.mtx
#	./testazny < ../Matrix/ParUTst/lpi_woodinfe/lpi_woodinfe.mtx
#	./testazny < ../Matrix/ParUTst/lpi_ex72a/lpi_ex72a.mtx
#
#	./testazny < ../Matrix/n3c4-b4.mtx    	# m >=n 
#	./testazny < ../Matrix/permuted_b1_ss.mtx  #7x7
#	./testazny < ../Matrix/problem.mtx      # m >=n
#	./testazny < ../Matrix/pwr01b.mtx	    #sample 39x39
#	./testazny < ../Matrix/r2.mtx           # m >=n failure
#	./testazny < ../Matrix/Ragusa16.mtx    #EMPTY row
#	./testazny < ../Matrix/s32.mtx		   #3x2
#	./testazny < ../Matrix/Tina_AskCal.mtx       #EMPTY col
#	./testazny < ../Matrix/Tina_AskCal_perm.mtx  #EMPTY col
#
	./testazny < ../Matrix/west0067.mtx 	#67x67!
	./testazny < ../Matrix/ParUTst/west0132/west0132.mtx
	./testazny < ../Matrix/ParUTst/west0156/west0156.mtx
	./testazny < ../Matrix/ParUTst/west0167/west0167.mtx
	./testazny < ../Matrix/ParUTst/west0381/west0381.mtx
	./testazny < ../Matrix/ParUTst/west0479/west0479.mtx
	./testazny < ../Matrix/ParUTst/west0497/west0497.mtx
	./testazny < ../Matrix/ParUTst/west0655/west0655.mtx
	./testazny < ../Matrix/ParUTst/west0989/west0989.mtx #!


	./testazny < ../Matrix/young1c.mtx     	#841x841
	./testazny < ../Matrix/ParUTst/az11_11.mtx   #14x14
	./testazny < ../Matrix/ParUTst/az88_2.mtx	 #4x4
	./testazny < ../Matrix/ParUTst/bcsstk01.mtx #NO CBs 48x48 
	./testazny < ../Matrix/ParUTst/bcsstk03.mtx  
	./testazny < ../Matrix/ParUTst/bcsstk22.mtx 
	./testazny < ../Matrix/ParUTst/bcsstm01.mtx	
	./testazny < ../Matrix/ParUTst/bcsstm02.mtx	
	./testazny < ../Matrix/ParUTst/bcsstm05.mtx	  		
	./testazny < ../Matrix/ParUTst/bcsstm06.mtx	  
	./testazny < ../Matrix/ParUTst/bcsstm19.mtx	  
	./testazny < ../Matrix/ParUTst/bcsstm20.mtx	  
	./testazny < ../Matrix/ParUTst/ex5.mtx  	#27x27
	./testazny < ../Matrix/ParUTst/mesh1e1.mtx  #48x48
	./testazny < ../Matrix/ParUTst/LF10/LF10.mtx #18x18
	./testazny < ../Matrix/ParUTst/mesh1em1/mesh1em1.mtx #48x48
	./testazny < ../Matrix/ParUTst/Trefethen_20/Trefethen_20.mtx #20x20
	./testazny < ../Matrix/ParUTst/Trefethen_20b/Trefethen_20b.mtx #19x19
	./testazny < ../Matrix/ParUTst/curtis54/curtis54.mtx #54x54
	./testazny < ../Matrix/ParUTst/nos1/nos1.mtx #237x237
	./testazny < ../Matrix/ParUTst/nos4/nos4.mtx  #100x100
	./testazny < ../Matrix/ParUTst/494_bus/494_bus.mtx #494x494
	./testazny < ../Matrix/ParUTst/mesh1em6/mesh1em6.mtx #48x48



#VGFLAGS = --tool=memcheck --leak-check=yes --show-reachable=yes -v\
#		  --track-origins=yes --num-callers=20 \
#		  --track-fds=yes --log-file="valgrindlogfile.log"
VGFLAGS =  --leak-check=yes  --log-file="valgrindlogfile.log"




val:
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lp_share1b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/c32.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/LFAT5.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a2.mtx 
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/b1_ss.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/cage3.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lfat5b.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/west0067.mtx 
#	valgrind	$(VGFLAGS)  ./testazny < ../Matrix/ParUTst/west0132/west0132.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0156/west0156.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0167/west0167.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0381/west0381.mtx 
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0479/west0479.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0497/west0497.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0655/west0655.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/ParUTst/west0989/west0989.mtx
#
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/cage5.mtx
	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lp_e226_transposed.mtx   #427x223
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
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/c32.mtx 
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/Groebner_id2003_aug.mtx  
#
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/a04.mtx
#	valgrind	$(VGFLAGS)	./testazny < ../Matrix/lp_e226.mtx

debug:
	$(C) $(FLAG) -g testazny.cpp -o testazny $(LIBS)
	gdb ./testazny < ../Matrix/c32.mtx		#3x2 	 #OK

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
