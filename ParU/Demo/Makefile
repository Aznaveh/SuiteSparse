default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CLIB = $(LDFLAGS) -L../../lib -lparu -lumfpack -lsuitesparseconfig -lcholmod -lamd \
        -lcolamd $(LIB_WITH_PARTITION)  $(LDLIBS)
FLIB = $(LAPACK) $(BLAS)


I = -I../../include  -I ../Include -I../../UMFPACK/Source -I../../AMD/Include


C = $(CXX) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) $(CONFIG_GPU) $(I) \
	$(CHOLMOD_CONFIG)

LIBS = $(CLIB) $(FLIB) $(TBB) $(GPULIB) 
#-lasan

FLAG = -Wno-write-strings 
#-O0 -fsanitize=address -g

all:
	$(C) $(FLAG) umf_test.cpp -o umfout $(LIBS)


run:
	./umfout < ../Matrix/b1_ss.mtx  	#7x7
#	./umfout < ../Matrix/ParUTst/tmp.mtx

VGFLAGS =  --leak-check=yes  --log-file="valgrindlogfile.log"\
		   --track-origins=yes\
		   #--show-leak-kinds=all

			#comment leak-check=yes or full

val:

#	valgrind	$(VGFLAGS)	./umfout < ../Matrix/ParUTst/tmp.mtx
	valgrind	$(VGFLAGS)	./umfout < ../Matrix/west0067.mtx 
#	valgrind	$(VGFLAGS)	./umfout < ../Matrix/b1_ss.mtx

	
debug:
	$(C) $(FLAG) -g umf_test.cpp -o umfout $(LIBS)
	gdb ./umfout < ../Matrix/c32.mtx		#3x2 	 #OK

gprofile:
	$(C) $(FLAG) -no-pie -O0 -g -pg umf_test.cpp -o gprof.o $(LIBS)
	./gprof.o < ../Matrix/west0067.mtx 
	gprof gprof.o > analysis.txt


vprofile:

	$(C) $(FLAG) -g umf_test.cpp -o vprof $(LIBS)
#	valgrind	--tool=callgrind -v ./vprof< ../Matrix/west0067.mtx
	valgrind	--tool=callgrind -v ./vprof< ../Matrix/ParUTst/tmp.mtx

purge: distclean

distclean: clean
	- $(RM) umfout analysis.txt
	- $(RM) gprof.o vprof gmon.out analysis.txt callgrind.out.*
	- $(RM) paneltest
	- $(RM) *.dot pfile tfile
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

az:
	(cd ..; $(MAKE) )
	(cd Demo; $(MAKE) run )
