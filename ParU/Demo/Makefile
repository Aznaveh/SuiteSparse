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

#FLAG = -Wno-write-strings -std=c++11
FLAG = -Wno-write-strings 
#-O0 -fsanitize=address -g
C11 = -std=c++11
C17 = -std=c++17
C20 = -std=c++20

all:
	$(C) $(FLAG) $(C11) paru_demo.cpp -o x_paru_demo $(LIBS)


run:
	./paru_demo < ../Matrix/b1_ss.mtx  	#7x7

VGFLAGS = --leak-check=yes  --log-file="valgrindlogfile.log"\
		   --show-leak-kinds=all --track-origins=yes 

			#comment leak-check=yes or full

val:

	valgrind	$(VGFLAGS)	./x_paru_demo < ../Matrix/ParUTst/tmp.mtx
#	valgrind	$(VGFLAGS)	./x_paru_demo < ../Matrix/west0067.mtx 
#	valgrind	$(VGFLAGS)	./x_paru_demo < ../Matrix/b1_ss.mtx

	
	#compling everything with -g
debug:
	(cd ../Lib; $(MAKE) purge; $(MAKE) DEBUG=yes)
	(cd ../Demo;)
	$(C) $(FLAG) -g paru_demo.cpp -o x_paru_demo $(LIBS)

	#compling everything with -pg
gprof:
	(cd ../Lib; $(MAKE) purge; $(MAKE) PROF=yes)
	$(C) $(FLAG) -no-pie -O0  -pg paru_demo.cpp -o gprof.o $(LIBS)
	./gprof.o < ../Matrix/ParUTst/tmp.mtx
	gprof gprof.o gmon.out > analysis.txt


cchg: debug
	valgrind	--tool=cachegrind -v ./x_paru_demo < ../Matrix/ParUTst/tmp.mtx
	
cllg: debug
	valgrind	--tool=callgrind -v ./x_paru_demo < ../Matrix/ParUTst/tmp.mtx

mssf: debug
	valgrind	--tool=massif -v ./x_paru_demo --log-file="valgrindlogfile.log"\
		< ../Matrix/ParUTst/tmp.mtx


tcov:
	$(C) $(FLAG) -Wall --coverage paru_demo.cpp -o tcov $(LIBS)
	./tcov < ../Matrix/west0067.mtx 
	gcov ../Source/*.cpp paru_demo.cpp 

purge: distclean

distclean: clean
	- $(RM) umfout 
	- $(RM) x_paru_demo
	- $(RM) analysis.txt tcov gprof
	- $(RM) gprof.o vprof gmon.out callgrind.out.*
	- $(RM) paneltest
	- $(RM) *.dot pfile tfile
	- $(RM) -r $(PURGE)
	- $(RM) test* 
	- $(RM) gumfout 

clean:
	- $(RM) -r $(CLEAN)

az:
	(cd ..; $(MAKE) )
	(cd Demo; $(MAKE) run )
rev:
	$(foreach var,$(.VARIABLES),$(info $(var) = $($(var)))) #print all var
	#@echo "$(BLAS)" #print some val

.PHONY: test
test:
	- rm test
	$(C) $(FLAG) $(C17) -g blas_test.cpp -o test $(LIBS)
#	$(C) $(FLAG) $(C11) blas_NoLoopAlloc.cpp -o test0 $(LIBS)
#	$(C) $(FLAG) $(C11) blas_1LoopAlloc.cpp -o test1 $(LIBS)
#	$(C) $(FLAG) $(C11) blas_3LoopAlloc.cpp -o test3 $(LIBS)
valtest:
	valgrind --leak-check=yes  --log-file="valgrindlogfile.log"\
		   --show-leak-kinds=all --track-origins=yes ./test	
